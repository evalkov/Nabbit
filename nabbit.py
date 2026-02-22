#!/usr/bin/env python3
"""
Nabbit — Nanobody NGS Enrichment Analysis Pipeline
=================================================
HPC pipeline for analyzing nanobody phage/yeast display enrichment
from paired-end Illumina sequencing data (fastq.gz).

Features:
  - Paired-end read merging, translation, quality filtering
  - Enrichment analysis with fold-change, CPM, and log2FC metrics
  - CDR annotation (IMGT-based VHH motif anchoring)
  - Clone clustering by CDR3 identity/similarity
  - Cross-round clone tracking with chimera detection
  - Saturation/rarefaction analysis with diversity estimation
  - Productivity statistics per round
  - AbLang pseudo-perplexity scoring (optional, --ablang)
  - IgBLAST V/D/J gene annotation (optional, --igblast)
  - Cluster PCA and CDR3 phylogenetic tree (optional, biopython)

Dependencies:
  module load python/3.11
  pip install --user numpy pandas scipy matplotlib
  pip install --user ablang        # optional, for --ablang flag
  pip install --user biopython scikit-learn  # optional, for PCA/tree
  conda install -c bioconda igblast          # optional, for --igblast

Usage:
  python nabbit.py \\
      --fastq-dir /data/$USER/nanobody_sequencing \\
      --output-dir /data/$USER/results \\
      --threads 8 --ablang --igblast
"""

import argparse
import gzip
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import uuid
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import entropy

# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


# ===========================================================================
# Utility functions
# ===========================================================================
def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]

def translate_dna(dna: str, frame: int = 0) -> str:
    protein = []
    for i in range(frame, len(dna) - 2, 3):
        codon = dna[i:i+3].upper()
        protein.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(protein)

def mean_phred(qual_string: str) -> float:
    if not qual_string:
        return 0.0
    return np.mean([ord(c) - 33 for c in qual_string])

def read_fastq(filepath: str):
    opener = gzip.open if filepath.endswith('.gz') else open
    with opener(filepath, 'rt') as f:
        batch = []
        for line in f:
            batch.append(line)
            if len(batch) == 4:
                yield batch[0].strip()[1:], batch[1].strip(), batch[3].strip()
                batch = []

def find_best_overlap_merge(r1: str, r2: str, min_overlap: int = 20) -> Tuple[Optional[str], int]:
    r2_rc = reverse_complement(r2)
    max_possible = min(len(r1), len(r2_rc))
    for overlap in range(max_possible, min_overlap - 1, -1):
        if r1[-overlap:] == r2_rc[:overlap]:
            return r1 + r2_rc[overlap:], overlap
    return None, 0

def write_fasta(entries: List[Tuple[str, str]], filepath: str):
    with open(filepath, 'w') as f:
        for header, seq in entries:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


# ===========================================================================
# FASTQ discovery
# ===========================================================================
def discover_fastqs(fastq_dir: str, sample_map: Optional[str] = None,
                    rounds: Optional[List[int]] = None) -> Dict[int, Dict[str, str]]:
    if sample_map:
        log.info(f"Reading sample map: {sample_map}")
        mapping = {}
        with open(sample_map) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    mapping[int(parts[0])] = {'R1': parts[1], 'R2': parts[2]}
        if rounds:
            mapping = {k: v for k, v in mapping.items() if k in rounds}
        return mapping

    fastq_dir = Path(fastq_dir)
    pattern = re.compile(r'(?P<num>\d+).*?R(?P<read>[12])', re.IGNORECASE)
    files: Dict[int, Dict[str, str]] = {}
    for ext in ('*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq'):
        for fp in fastq_dir.rglob(ext):
            m = pattern.search(fp.name)
            if m:
                rnd = int(m.group('num'))
                files.setdefault(rnd, {})[f"R{m.group('read')}"] = str(fp.resolve())

    valid = {}
    for rnd, reads in sorted(files.items()):
        if 'R1' in reads and 'R2' in reads:
            if rounds is None or rnd in rounds:
                valid[rnd] = reads
            else:
                log.debug(f"Round {rnd}: skipped (not in --rounds)")
        else:
            log.warning(f"Round {rnd}: missing pair, skipping")

    # Fallback: try _1/_2 suffix pattern (SRA, ENA, etc.)
    if not valid:
        sra_pattern = re.compile(r'^(.+)[_.]([12])\.(?:fastq|fq)(?:\.gz)?$', re.IGNORECASE)
        sra_groups: Dict[str, Dict[str, str]] = {}
        for ext in ('*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq'):
            for fp in fastq_dir.rglob(ext):
                m = sra_pattern.match(fp.name)
                if m:
                    prefix = m.group(1)
                    read = m.group(2)
                    sra_groups.setdefault(prefix, {})[f"R{read}"] = str(fp.resolve())
        for rnd_idx, prefix in enumerate(sorted(sra_groups), start=1):
            reads = sra_groups[prefix]
            if 'R1' in reads and 'R2' in reads:
                if rounds is None or rnd_idx in rounds:
                    valid[rnd_idx] = reads
                    log.info(f"SRA fallback: {prefix} → Round {rnd_idx}")

    return valid


# ===========================================================================
# Read processing
# ===========================================================================
def process_fastq_pair(r1_file: str, r2_file: str, min_overlap: int = 20,
                       min_prot_len: int = 90, max_prot_len: int = 180,
                       min_avg_qual: float = 20.0) -> Tuple[Counter, dict]:
    proteins = []
    total = merged = qual_fail = stop_fail = len_fail = 0

    for (_, s1, q1), (_, s2, q2) in zip(read_fastq(r1_file), read_fastq(r2_file)):
        total += 1
        if min_avg_qual > 0 and (mean_phred(q1) + mean_phred(q2)) / 2 < min_avg_qual:
            qual_fail += 1
            continue
        seq, ovlp = find_best_overlap_merge(s1, s2, min_overlap)
        if seq is None:
            continue
        merged += 1
        protein = translate_dna(seq)
        if '*' in protein:
            stop_fail += 1
            continue
        if not (min_prot_len <= len(protein) <= max_prot_len):
            len_fail += 1
            continue
        proteins.append(protein)

    stats = dict(total_reads=total, merged=merged, qual_fail=qual_fail,
                 stop_codon_fail=stop_fail, length_fail=len_fail,
                 valid_proteins=len(proteins))
    return Counter(proteins), stats

def _process_round(args):
    rnd, r1, r2, mo, mnp, mxp, mq = args
    counter, stats = process_fastq_pair(r1, r2, mo, mnp, mxp, mq)
    return rnd, counter, stats


# ===========================================================================
# CDR Annotation
# ===========================================================================
class VHHAnnotator:
    """
    Rule-based VHH CDR annotator using conserved framework anchor motifs.
    
    Anchors:
      FR1→CDR1: Conserved Cys ~pos 22 + 4 residues
      CDR1→FR2: W[FY]RQ motif
      FR2→CDR2: ~15 residue FR2
      CDR2→FR3: [YF][YF]DSVKG motif
      FR3→CDR3: YYC motif
      CDR3→FR4: WG[QK]GT motif
    """

    @staticmethod
    def annotate(sequence: str) -> dict:
        result = dict(FR1='', CDR1='', FR2='', CDR2='', FR3='', CDR3='', FR4='',
                      CDR3_length=0, annotation_quality='good')
        seq = sequence.upper()

        # Find FR4 anchor (WG[QKR]G[TLS])
        fr4_match = None
        for m in re.finditer(r'WG[QKRE]G[TLS]', seq):
            fr4_match = m
        if fr4_match is None:
            for m in re.finditer(r'WG.G', seq):
                fr4_match = m
        if fr4_match is None:
            result['annotation_quality'] = 'failed_no_FR4'
            result['FR1'] = seq
            return result
        fr4_start = fr4_match.start()

        # Find CDR3 start (YYC anchor)
        cdr3_anchor = None
        for m in re.finditer(r'Y[YF]C', seq[:fr4_start]):
            cdr3_anchor = m
        if cdr3_anchor:
            cdr3_start = cdr3_anchor.end()
        else:
            for i in range(fr4_start - 1, max(fr4_start - 20, 0), -1):
                if seq[i] == 'C':
                    cdr3_start = i + 1
                    break
            else:
                cdr3_start = max(0, fr4_start - 15)

        n_term = seq[:cdr3_start]

        # Find FR1/CDR1 boundary (conserved Cys ~pos 22)
        cys_pos = None
        for i in range(15, min(35, len(n_term))):
            if n_term[i] == 'C':
                cys_pos = i
                break
        if cys_pos is None:
            result['annotation_quality'] = 'partial_no_Cys22'
            cys_pos = min(21, len(n_term) - 1)
        fr1_end = cys_pos + 4

        # Find CDR1/FR2 boundary (WFRQ anchor)
        fr2_anchor = None
        for m in re.finditer(r'W[FYWLV]R[QKR]', n_term[fr1_end:]):
            fr2_anchor = m
            break
        if fr2_anchor is None:
            for m in re.finditer(r'W[FYWLV]R', n_term[fr1_end:]):
                fr2_anchor = m
                break
        if fr2_anchor:
            cdr1_end = fr1_end + fr2_anchor.start()
        else:
            cdr1_end = min(fr1_end + 12, len(n_term))
            result['annotation_quality'] = 'partial_no_FR2_anchor'

        fr2_end = min(cdr1_end + 15, len(n_term))

        # Find CDR2/FR3 boundary
        fr3_anchor = None
        for pat in (r'[YF][YFA]D[SA]VKG', r'[YF]D[SA]VKG', r'D[SA][VI]KG'):
            for m in re.finditer(pat, n_term[fr2_end:]):
                fr3_anchor = m
                break
            if fr3_anchor:
                break
        if fr3_anchor:
            cdr2_end = fr2_end + fr3_anchor.start()
        else:
            cdr2_end = min(fr2_end + 10, len(n_term))
            if result['annotation_quality'] == 'good':
                result['annotation_quality'] = 'partial_no_FR3_anchor'

        result['FR1'] = n_term[:fr1_end]
        result['CDR1'] = n_term[fr1_end:cdr1_end]
        result['FR2'] = n_term[cdr1_end:fr2_end]
        result['CDR2'] = n_term[fr2_end:cdr2_end]
        result['FR3'] = n_term[cdr2_end:cdr3_start]
        result['CDR3'] = seq[cdr3_start:fr4_start]
        result['FR4'] = seq[fr4_start:]
        result['CDR3_length'] = fr4_start - cdr3_start
        return result


# ===========================================================================
# Clone clustering (CDR3-based)
# ===========================================================================
def _nw_identity(s1: str, s2: str) -> float:
    """Quick NW identity between two short sequences (for CDR3 clustering)."""
    m, n = len(s1), len(s2)
    if m == 0 or n == 0:
        return 0.0
    score = np.zeros((m+1, n+1), dtype=np.int16)
    for i in range(1, m+1):
        score[i, 0] = i * -2
    for j in range(1, n+1):
        score[0, j] = j * -2
    for i in range(1, m+1):
        for j in range(1, n+1):
            score[i, j] = max(
                score[i-1, j-1] + (2 if s1[i-1] == s2[j-1] else -1),
                score[i-1, j] - 2,
                score[i, j-1] - 2)
    # Traceback for identity count
    i, j = m, n
    matches = total = 0
    while i > 0 or j > 0:
        total += 1
        if i > 0 and j > 0:
            diag = score[i-1, j-1] + (2 if s1[i-1] == s2[j-1] else -1)
            if score[i, j] == diag:
                if s1[i-1] == s2[j-1]:
                    matches += 1
                i -= 1; j -= 1
                continue
        if i > 0 and score[i, j] == score[i-1, j] - 2:
            i -= 1
        else:
            j -= 1
    return matches / total * 100 if total > 0 else 0.0


def cluster_by_cdr3(annotations: pd.DataFrame, threshold: float = 0.8) -> pd.DataFrame:
    annotations = annotations.copy()

    # Exact CDR3 → clone_id
    cdr3_to_clone = {}
    counter = 0
    clone_ids = []
    for cdr3 in annotations['CDR3']:
        if cdr3 not in cdr3_to_clone:
            counter += 1
            cdr3_to_clone[cdr3] = f"clone_{counter:04d}"
        clone_ids.append(cdr3_to_clone[cdr3])
    annotations['clone_id'] = clone_ids

    # Similarity-based families (single-linkage)
    unique_cdr3s = list(cdr3_to_clone.keys())
    n = len(unique_cdr3s)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for i, j in combinations(range(n), 2):
        c1, c2 = unique_cdr3s[i], unique_cdr3s[j]
        if not c1 or not c2:
            continue
        if abs(len(c1) - len(c2)) > max(len(c1), len(c2)) * (1 - threshold):
            continue
        if _nw_identity(c1, c2) >= threshold * 100:
            union(i, j)

    fam_map = {}
    fam_counter = 0
    cdr3_to_fam = {}
    for i in range(n):
        root = find(i)
        if root not in fam_map:
            fam_counter += 1
            fam_map[root] = f"family_{fam_counter:03d}"
        cdr3_to_fam[unique_cdr3s[i]] = fam_map[root]

    annotations['clone_family'] = [cdr3_to_fam.get(c, 'unassigned') for c in annotations['CDR3']]
    return annotations


# ===========================================================================
# Chimera detection
# ===========================================================================
def detect_chimeras(annotations: pd.DataFrame) -> pd.DataFrame:
    annotations = annotations.copy()
    fr_sig = (annotations['FR1'].fillna('') + '|' +
              annotations['FR2'].fillna('') + '|' +
              annotations['FR3'].fillna(''))

    cdr3_by_fr = defaultdict(set)
    fr_by_cdr3 = defaultdict(set)
    for sig, cdr3 in zip(fr_sig, annotations['CDR3']):
        cdr3_by_fr[sig].add(cdr3)
        fr_by_cdr3[cdr3].add(sig)

    flags, notes = [], []
    for sig, cdr3 in zip(fr_sig, annotations['CDR3']):
        n_cdr3 = len(cdr3_by_fr[sig])
        n_fr = len(fr_by_cdr3[cdr3])
        if n_cdr3 > 1 and n_fr > 1:
            flags.append(True)
            notes.append(f"FR shared with {n_cdr3} CDR3s; CDR3 shared with {n_fr} FRs")
        elif n_cdr3 > 1:
            flags.append(False)
            notes.append(f"Same FR, {n_cdr3} different CDR3s")
        elif n_fr > 1:
            flags.append(False)
            notes.append(f"Same CDR3, {n_fr} different FRs")
        else:
            flags.append(False)
            notes.append('')

    annotations['is_potential_chimera'] = flags
    annotations['chimera_note'] = notes
    return annotations


# ===========================================================================
# Cross-round clone tracking
# ===========================================================================
def track_clones(rounds_data: Dict[str, Counter], annotations: pd.DataFrame,
                 top_n: int = 50) -> pd.DataFrame:
    round_names = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
    annotator = VHHAnnotator()

    seq_to_family = dict(zip(annotations['Sequence'], annotations['clone_family']))
    seq_to_cdr3 = dict(zip(annotations['Sequence'], annotations['CDR3']))

    family_counts = defaultdict(lambda: defaultdict(int))
    family_cdr3 = {}

    for rnd in round_names:
        total = sum(rounds_data[rnd].values())
        if total == 0:
            continue
        for seq, count in rounds_data[rnd].most_common(top_n * 10):
            if seq not in seq_to_family:
                ann = annotator.annotate(seq)
                cdr3 = ann['CDR3']
                matched = False
                for existing_cdr3, fam in dict(zip(annotations['CDR3'], annotations['clone_family'])).items():
                    if cdr3 == existing_cdr3:
                        seq_to_family[seq] = fam
                        seq_to_cdr3[seq] = cdr3
                        matched = True
                        break
                if not matched:
                    seq_to_family[seq] = f"novel_{hash(seq) % 10000:04d}"
                    seq_to_cdr3[seq] = cdr3

            fam = seq_to_family[seq]
            family_counts[fam][rnd] += count
            if fam not in family_cdr3:
                family_cdr3[fam] = seq_to_cdr3.get(seq, '')

    rows = []
    for fam in family_counts:
        row = {'clone_family': fam, 'representative_CDR3': family_cdr3.get(fam, '')}
        for rnd in round_names:
            total = sum(rounds_data[rnd].values())
            c = family_counts[fam].get(rnd, 0)
            freq = c / total if total > 0 else 0
            row[f'{rnd}_count'] = c
            row[f'{rnd}_freq'] = freq
            row[f'{rnd}_cpm'] = freq * 1_000_000
        f1 = row.get(f'{round_names[0]}_freq', 0)
        fl = row.get(f'{round_names[-1]}_freq', 0)
        row['enrichment'] = fl / f1 if f1 > 0 else (float('inf') if fl > 0 else 0)
        first_cpm = row.get(f'{round_names[0]}_cpm', 0)
        last_cpm = row.get(f'{round_names[-1]}_cpm', 0)
        row['log2fc'] = np.log2((last_cpm + 10) / (first_cpm + 10))
        rows.append(row)

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(f'{round_names[-1]}_freq', ascending=False).reset_index(drop=True)
    return df


# ===========================================================================
# AbLang scoring
# ===========================================================================
def score_with_ablang(sequences: List[str], device: str = 'cpu') -> pd.DataFrame:
    """
    Score sequences with AbLang heavy-chain model.
    Returns pseudo-perplexity (lower = more antibody-like) and
    per-position likelihood statistics.
    """
    try:
        import ablang
    except ImportError:
        log.error("AbLang not installed. Install with: pip install --user ablang")
        log.error("Skipping AbLang scoring.")
        return pd.DataFrame()

    log.info(f"Loading AbLang heavy-chain model (device={device})...")
    model = ablang.pretrained("heavy")
    model.freeze()

    log.info(f"Scoring {len(sequences)} sequences...")

    # Get res-likelihoods (logits, shape: batch x seq_len+2 x 20)
    likelihoods = model(sequences, mode='likelihood')

    # Build AA-to-column mapping: vocab indices 1-20 are the 20 amino acids
    vocab = model.tokenizer.vocab_to_aa
    aa_order = [vocab[i] for i in range(1, 21)]
    aa_to_col = {aa: col for col, aa in enumerate(aa_order)}

    # Convert logits to probabilities
    from scipy.special import softmax
    probs_all = softmax(likelihoods, axis=-1)

    rows = []
    for i, seq in enumerate(sequences):
        probs = probs_all[i]  # shape: (seq_len+2, 20)

        log_probs = []
        for pos, aa in enumerate(seq):
            col = aa_to_col.get(aa)
            if col is not None and (pos + 1) < len(probs):
                prob = float(probs[pos + 1][col])  # +1 to skip start token
                if prob > 0:
                    log_probs.append(np.log(prob))

        if log_probs:
            mean_nll = -np.mean(log_probs)
            pseudo_ppl = np.exp(mean_nll)
            mean_prob = np.exp(np.mean(log_probs))
            min_prob = np.exp(min(log_probs))
        else:
            pseudo_ppl = float('inf')
            mean_prob = 0
            min_prob = 0

        rows.append({
            'Sequence': seq,
            'ablang_pseudo_ppl': pseudo_ppl,
            'ablang_mean_prob': mean_prob,
            'ablang_min_prob': min_prob,
            'ablang_n_positions': len(log_probs),
        })

        if (i + 1) % 10 == 0:
            log.info(f"  Scored {i+1}/{len(sequences)}")

    return pd.DataFrame(rows)


# ===========================================================================
# Enrichment analysis
# ===========================================================================
def compute_enrichment(rounds_data: Dict[str, Counter], top_n: int = 50):
    round_names = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
    all_seqs = set()
    for rd in rounds_data.values():
        all_seqs.update(rd.keys())
    log.info(f"Total unique sequences across all rounds: {len(all_seqs):,}")

    totals = {r: sum(rounds_data[r].values()) for r in round_names}
    rows = []
    for seq in all_seqs:
        row = {'Sequence': seq, 'Length': len(seq)}
        for r in round_names:
            c = rounds_data[r][seq]
            row[f'{r}_count'] = c
            row[f'{r}_freq'] = c / totals[r] if totals[r] > 0 else 0
            row[f'{r}_cpm'] = (c / totals[r]) * 1_000_000 if totals[r] > 0 else 0
        rows.append(row)

    df = pd.DataFrame(rows)
    last_round = None
    for r in reversed(round_names):
        if totals.get(r, 0) > 0:
            last_round = r
            break
    if last_round is None:
        return df, df.head(0)

    df = df.sort_values(f'{last_round}_freq', ascending=False).reset_index(drop=True)

    # Fold-changes (linear)
    first = round_names[0]
    ps = 1e-8
    df[f'{first}_to_{last_round}_fold'] = (df[f'{last_round}_freq'] + ps) / (df[f'{first}_freq'] + ps)
    for i in range(len(round_names) - 1):
        rp, rn = round_names[i], round_names[i+1]
        df[f'{rp}_to_{rn}_fold'] = (df[f'{rn}_freq'] + ps) / (df[f'{rp}_freq'] + ps)

    # Log2 fold-changes (CPM-based, alpseq-style with pseudocount of 10)
    df[f'{first}_to_{last_round}_log2fc'] = np.log2(
        (df[f'{last_round}_cpm'] + 10) / (df[f'{first}_cpm'] + 10))
    for i in range(len(round_names) - 1):
        rp, rn = round_names[i], round_names[i+1]
        df[f'{rp}_to_{rn}_log2fc'] = np.log2(
            (df[f'{rn}_cpm'] + 10) / (df[f'{rp}_cpm'] + 10))

    return df, df.head(top_n).copy()


# ===========================================================================
# Diversity / convergence
# ===========================================================================
def compute_diversity(rounds_data):
    round_names = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
    rows = []
    for rnd in round_names:
        c = rounds_data[rnd]
        total = sum(c.values())
        if total == 0:
            continue
        freqs = np.array([v / total for v in c.values()])
        rows.append({
            'Round': rnd, 'Total_reads': total, 'Unique_sequences': len(c),
            'Shannon_diversity': entropy(freqs, base=2),
            'Simpson_index': float(np.sum(freqs**2)),
            'Top1_dominance_%': c.most_common(1)[0][1] / total * 100,
            'Top10_dominance_%': sum(v for _, v in c.most_common(10)) / total * 100,
        })
    return pd.DataFrame(rows)

def compute_convergence(rounds_data):
    round_names = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
    rows = []
    for i, j in combinations(range(len(round_names)), 2):
        r1, r2 = round_names[i], round_names[j]
        s1, s2 = set(rounds_data[r1].keys()), set(rounds_data[r2].keys())
        shared = len(s1 & s2)
        rows.append({
            'Round_A': r1, 'Round_B': r2,
            'Unique_A': len(s1), 'Unique_B': len(s2),
            'Shared': shared,
            'Jaccard': shared / len(s1 | s2) if s1 | s2 else 0,
        })
    return pd.DataFrame(rows)


# ===========================================================================
# Saturation / rarefaction analysis (alpseq-style)
# ===========================================================================
def compute_saturation(rounds_data: Dict[str, Counter], annotations: pd.DataFrame,
                       n_points: int = 100) -> pd.DataFrame:
    """Rarefaction analysis: subsample reads and count unique sequences at each depth."""
    round_names = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
    rows = []
    for rnd in round_names:
        counter = rounds_data[rnd]
        # Use unique protein sequences (weighted by read count), excluding singletons
        seq_counts = {seq: count for seq, count in counter.items() if count > 1}
        if not seq_counts:
            continue
        # Build weighted vector of sequence labels
        seq_vec = []
        for seq, count in seq_counts.items():
            seq_vec.extend([seq] * count)
        seq_arr = np.array(seq_vec)
        total = len(seq_arr)
        sample_sizes = np.linspace(1, total, n_points, dtype=int)
        sample_sizes = np.unique(sample_sizes)
        prev_unique = 0
        prev_size = 0
        for ss in sample_sizes:
            indices = np.random.choice(total, size=ss, replace=False)
            n_unique = len(set(seq_arr[indices]))
            deriv = (n_unique - prev_unique) / (ss - prev_size) if ss > prev_size else 0
            rows.append({
                'Round': rnd,
                'sample_size': int(ss),
                'number_unique': n_unique,
                'derivative': deriv,
            })
            prev_unique = n_unique
            prev_size = ss
    return pd.DataFrame(rows)


def estimate_diversity(saturation_df: pd.DataFrame) -> pd.DataFrame:
    """Estimate library diversity from saturation curves (asymptote at min positive derivative)."""
    rows = []
    for rnd, grp in saturation_df.groupby('Round'):
        pos_deriv = grp[grp['derivative'] > 0]
        if pos_deriv.empty:
            continue
        min_deriv_row = pos_deriv.loc[pos_deriv['derivative'].idxmin()]
        est = int(min_deriv_row['number_unique'])
        # Round appropriately
        if est > 1_000_000:
            est = round(est / 100_000) * 100_000
        elif est > 100_000:
            est = round(est / 10_000) * 10_000
        elif est > 10_000:
            est = round(est / 1_000) * 1_000
        elif est > 1_000:
            est = round(est / 100) * 100
        rows.append({'Round': rnd, 'approx_diversity': est})
    return pd.DataFrame(rows)


# ===========================================================================
# Diverse sequence selection (alpseq-style)
# ===========================================================================
def select_diverse_sequences(top_df: pd.DataFrame, annotations: pd.DataFrame,
                             n: int = 50) -> List[Tuple[str, str]]:
    """Select diverse top sequences: one representative per clone family, fill remaining."""
    fam_map = dict(zip(annotations['Sequence'], annotations['clone_family']))
    cdr3_map = dict(zip(annotations['Sequence'], annotations['CDR3']))
    clone_map = dict(zip(annotations['Sequence'], annotations['clone_id']))

    selected = []
    seen_families = set()

    # First pass: one per family
    for _, row in top_df.iterrows():
        if len(selected) >= n:
            break
        seq = row['Sequence']
        fam = fam_map.get(seq, 'unknown')
        if fam not in seen_families:
            seen_families.add(fam)
            cid = clone_map.get(seq, 'unk')
            cdr3 = cdr3_map.get(seq, '')
            header = f"rank{int(row.name)+1}_{cid}_fam{fam}_CDR3len{len(cdr3)}"
            selected.append((header, seq))

    # Second pass: fill remaining slots with next-highest unselected
    for _, row in top_df.iterrows():
        if len(selected) >= n:
            break
        seq = row['Sequence']
        if any(s[1] == seq for s in selected):
            continue
        cid = clone_map.get(seq, 'unk')
        cdr3 = cdr3_map.get(seq, '')
        header = f"rank{int(row.name)+1}_{cid}_CDR3len{len(cdr3)}"
        selected.append((header, seq))

    return selected


# ===========================================================================
# Productivity stats (alpseq-style)
# ===========================================================================
def compute_productivity_stats(all_stats: Dict[str, dict]) -> pd.DataFrame:
    """Reshape processing stats into per-round productivity breakdown."""
    round_names = sorted(all_stats.keys(), key=lambda x: int(x.replace('Round', '')))
    rows = []
    for rnd in round_names:
        s = all_stats[rnd]
        total = s['total_reads']
        rows.append({
            'Round': rnd,
            'total_reads': total,
            'merged': s['merged'],
            'merge_rate': s['merged'] / total * 100 if total > 0 else 0,
            'qual_fail': s['qual_fail'],
            'stop_codon_fail': s['stop_codon_fail'],
            'length_fail': s['length_fail'],
            'valid_proteins': s['valid_proteins'],
            'pct_productive': s['valid_proteins'] / total * 100 if total > 0 else 0,
        })
    return pd.DataFrame(rows)


# ===========================================================================
# IgBLAST annotation (optional)
# ===========================================================================
def _igblast_chunk(sequences: List[Tuple[str, str]], igblastp: str,
                   germline_db_v: str, igdata: str,
                   threads: int = 1) -> List[dict]:
    """
    Run igblastp on a chunk of protein sequences and parse the output.
    Returns a list of annotation dicts (one per query sequence).
    """
    tmp_fasta = None
    tmp_out = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_in:
            for header, seq in sequences:
                tmp_in.write(f">{header}\n{seq}\n")
            tmp_fasta = tmp_in.name

        tmp_out = tmp_fasta + '.igblast.txt'

        cmd = [
            igblastp,
            '-query', tmp_fasta,
            '-germline_db_V', germline_db_v,
            '-organism', 'alpaca',
            '-domain_system', 'imgt',
            '-num_threads', str(threads),
            '-outfmt', '7 qseqid sseqid pident evalue',
            '-out', tmp_out,
        ]
        env = os.environ.copy()
        env['IGDATA'] = igdata

        result = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)
        if result.returncode != 0:
            log.warning(f"IgBLAST returned non-zero exit code: {result.returncode}")
            log.warning(f"stderr: {result.stderr[:500]}")
            return []

        if not os.path.exists(tmp_out):
            log.warning("IgBLAST produced no output file")
            return []

        ann_rows = []
        current_query = None
        domain_info = {}
        v_hits = []

        with open(tmp_out) as f:
            for line in f:
                line = line.strip()
                if line.startswith('# Query:'):
                    if current_query is not None:
                        best_v = v_hits[0] if v_hits else ('', 0.0)
                        ann_rows.append({
                            'sequence_id': current_query,
                            'v_call': best_v[0],
                            'v_identity': best_v[1],
                            **domain_info,
                        })
                    current_query = line.split('# Query:')[1].strip()
                    domain_info = {}
                    v_hits = []
                elif line.startswith('# Alignment summary'):
                    pass
                elif not line.startswith('#') and not line.startswith('Total') and '\t' in line:
                    parts = line.split('\t')
                    if parts[0].endswith('-IMGT') or parts[0].startswith('CDR3'):
                        region = parts[0].split('-')[0].replace(' (germline)', '')
                        if len(parts) >= 8:
                            domain_info[f'{region}_identity'] = float(parts[7]) if parts[7] != 'N/A' else 0
                    elif parts[0] == 'V' and len(parts) >= 5:
                        v_hits.append((parts[2], float(parts[3])))

            if current_query is not None:
                best_v = v_hits[0] if v_hits else ('', 0.0)
                ann_rows.append({
                    'sequence_id': current_query,
                    'v_call': best_v[0],
                    'v_identity': best_v[1],
                    **domain_info,
                })

        return ann_rows

    except subprocess.TimeoutExpired:
        log.warning("IgBLAST chunk timed out after 300s")
        return []
    except Exception as e:
        log.warning(f"IgBLAST chunk failed: {e}")
        return []
    finally:
        for path in (tmp_fasta, tmp_out):
            if path and os.path.exists(path):
                os.unlink(path)


def run_igblast(sequences: List[Tuple[str, str]], igblast_db: str,
                threads: int = 1,
                chunk_size: int = 0) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run igblastp on protein sequences and parse the domain/hit output.
    sequences: list of (header, seq_aa) tuples.
    chunk_size: if >0, split sequences into chunks and run in parallel.
    Returns: (annotations_df, gene_usage_df)
    """
    igblastp = shutil.which('igblastp')
    if igblastp is None:
        # Check common locations outside PATH
        for candidate in [
            Path.home() / 'bin' / 'igblastp',
            Path('/private/tmp/ncbi-igblast-1.22.0/bin/igblastp'),
            Path('/usr/local/bin/igblastp'),
        ]:
            if candidate.is_file() and os.access(candidate, os.X_OK):
                igblastp = str(candidate)
                break
    if igblastp is None:
        log.warning("igblastp not found in PATH. Skipping IgBLAST annotation.")
        return pd.DataFrame(), pd.DataFrame()

    db_path = Path(igblast_db)
    germline_db_v_aa = str(db_path / 'databases' / 'imgt_alpaca_ighv_aa')
    germline_db_v_nt = str(db_path / 'databases' / 'imgt_alpaca_ighv')
    germline_db_v = germline_db_v_aa if os.path.exists(germline_db_v_aa + '.phr') else germline_db_v_nt
    igdata = str(db_path / 'igdata')

    # Decide whether to chunk
    if chunk_size > 0 and len(sequences) > chunk_size:
        chunks = [sequences[i:i + chunk_size] for i in range(0, len(sequences), chunk_size)]
        n_parallel = min(len(chunks), threads)
        threads_per_chunk = max(1, threads // n_parallel)
        log.info(f"  IgBLAST: splitting {len(sequences)} sequences into {len(chunks)} chunks "
                 f"({n_parallel} parallel, {threads_per_chunk} threads each)")

        ann_rows = []
        with ThreadPoolExecutor(max_workers=n_parallel) as executor:
            futures = {
                executor.submit(_igblast_chunk, chunk, igblastp, germline_db_v, igdata, threads_per_chunk): idx
                for idx, chunk in enumerate(chunks)
            }
            for future in as_completed(futures):
                idx = futures[future]
                try:
                    rows = future.result()
                    ann_rows.extend(rows)
                    log.info(f"  IgBLAST chunk {idx + 1}/{len(chunks)}: {len(rows)} annotations")
                except Exception as e:
                    log.warning(f"  IgBLAST chunk {idx + 1}/{len(chunks)} failed: {e}")
    else:
        ann_rows = _igblast_chunk(sequences, igblastp, germline_db_v, igdata, threads)

    if not ann_rows:
        return pd.DataFrame(), pd.DataFrame()

    annotations_df = pd.DataFrame(ann_rows)

    # Compute V gene usage percentages
    usage_rows = []
    v_calls = annotations_df['v_call'].dropna()
    v_calls = v_calls[v_calls != '']
    total = len(v_calls)
    if total > 0:
        for gene, count in v_calls.value_counts().items():
            usage_rows.append({
                'gene_type': 'V',
                'gene': gene,
                'count': count,
                'percentage': count / total * 100,
            })
    gene_usage_df = pd.DataFrame(usage_rows)

    return annotations_df, gene_usage_df


# ===========================================================================
# Cluster PCA (requires biopython + scikit-learn)
# ===========================================================================
def compute_cluster_pca(top_df: pd.DataFrame,
                        annotations: pd.DataFrame) -> Optional[pd.DataFrame]:
    """PCA of enriched clusters using pairwise alignment scores."""
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
        from sklearn.decomposition import PCA
    except ImportError:
        log.warning("biopython/scikit-learn not available, skipping cluster PCA")
        return None

    seqs = top_df['Sequence'].tolist()
    n = len(seqs)
    if n < 3:
        return None

    # Cap at 100 sequences for performance
    seqs = seqs[:100]
    n = len(seqs)

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    # Build pairwise alignment score matrix
    score_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            s = aligner.score(seqs[i], seqs[j])
            score_mat[i, j] = s
            score_mat[j, i] = s

    # Convert similarity to distance
    max_score = score_mat.max()
    dist_mat = max_score - score_mat

    pca = PCA(n_components=2)
    coords = pca.fit_transform(dist_mat)

    fam_map = dict(zip(annotations['Sequence'], annotations['clone_family']))
    round_names = sorted(top_df.columns[top_df.columns.str.endswith('_cpm')].tolist())
    last_cpm_col = round_names[-1] if round_names else None
    first_round = sorted(
        [c for c in top_df.columns if c.endswith('_cpm')],
        key=lambda x: int(x.replace('_cpm', '').replace('Round', ''))
    )
    log2fc_col = [c for c in top_df.columns if c.endswith('_log2fc')]
    last_log2fc = log2fc_col[0] if log2fc_col else None

    cdr3_map = dict(zip(annotations['Sequence'], annotations['CDR3'])) if 'CDR3' in annotations.columns else {}
    rows = []
    for i in range(n):
        seq = seqs[i]
        row_data = top_df[top_df['Sequence'] == seq].iloc[0] if seq in top_df['Sequence'].values else None
        rows.append({
            'Sequence': seq,
            'rank': i + 1,
            'CDR3': cdr3_map.get(seq, ''),
            'PC1': coords[i, 0],
            'PC2': coords[i, 1],
            'clone_family': fam_map.get(seq, 'unknown'),
            'enrichment': float(row_data[last_log2fc]) if row_data is not None and last_log2fc else 0,
            'cpm': float(row_data[last_cpm_col]) if row_data is not None and last_cpm_col else 0,
        })

    return pd.DataFrame(rows)


# ===========================================================================
# CDR3 phylogenetic tree (requires biopython)
# ===========================================================================
def build_cdr3_tree(annotations: pd.DataFrame, output_dir: str, max_seqs: int = 50):
    """Build neighbor-joining tree of top CDR3 sequences."""
    try:
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Align import PairwiseAligner
        import Bio.Phylo as Phylo
    except ImportError:
        log.warning("biopython not available, skipping CDR3 tree")
        return

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        log.warning("matplotlib not available, skipping CDR3 tree plot")
        return

    cdr3s = annotations[annotations['CDR3'].str.len() > 0].head(max_seqs)
    if len(cdr3s) < 3:
        return

    seqs = cdr3s['CDR3'].tolist()
    labels = [f"{row['clone_id']}_{row['clone_family']}" for _, row in cdr3s.iterrows()]

    # Pad sequences to same length for MSA
    max_len = max(len(s) for s in seqs)
    records = []
    seen_ids = set()
    for i, (seq, label) in enumerate(zip(seqs, labels)):
        padded = seq + '-' * (max_len - len(seq))
        # Sanitize label and ensure uniqueness
        safe_label = re.sub(r'[^a-zA-Z0-9_]', '_', label)[:25]
        base_label = safe_label
        suffix = 1
        while safe_label in seen_ids:
            safe_label = f"{base_label}_{suffix}"
            suffix += 1
        seen_ids.add(safe_label)
        records.append(SeqRecord(Seq(padded), id=safe_label, description=''))

    aln = MultipleSeqAlignment(records)

    # Distance matrix and NJ tree
    calculator = DistanceCalculator('blosum62')
    try:
        dm = calculator.get_distance(aln)
    except Exception:
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    fig, ax = plt.subplots(figsize=(12, max(8, len(seqs) * 0.3)))
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.set_title('CDR3 Phylogenetic Tree (Neighbor-Joining)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cdr3_phylo_tree.png'), dpi=300, bbox_inches='tight')
    plt.close()
    log.info("CDR3 phylogenetic tree saved")


# ===========================================================================
# Plots
# ===========================================================================
def generate_plots(top_df, diversity_df, rounds_data, annotations, output_dir,
                    ablang_df=None, clone_tracking=None, all_stats=None,
                    saturation_df=None, diversity_est=None, igblast_df=None,
                    gene_usage_df=None, cluster_pca=None):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        log.warning("matplotlib not available, skipping plots")
        return

    round_names = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
    rds = [r for r in round_names if sum(rounds_data[r].values()) > 0]
    labels = [r.replace('Round', 'R') for r in rds]

    # Enrichment trajectories
    fig, ax = plt.subplots(figsize=(12, 7))
    for idx, row in top_df.head(10).iterrows():
        freqs = [row.get(f'{r}_freq', 0) * 100 for r in rds]
        ax.plot(labels, freqs, marker='o', linewidth=2, label=f"Rank {idx+1} ({int(row['Length'])}aa)")
    ax.set_xlabel('Selection Round', fontsize=13, fontweight='bold')
    ax.set_ylabel('Relative Frequency (%)', fontsize=13, fontweight='bold')
    ax.set_title('Enrichment Trajectories — Top 10', fontsize=15, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'enrichment_trajectories.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Top sequence detail
    if len(top_df) > 0:
        tr = top_df.iloc[0]
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5))
        a1.plot(labels, [tr.get(f'{r}_count', 0) for r in rds], marker='o', linewidth=2, markersize=10, color='#2E86AB')
        a1.set_xlabel('Round'); a1.set_ylabel('Read Count'); a1.set_title('Top Enriched: Counts'); a1.grid(True, alpha=0.3)
        a2.plot(labels, [tr.get(f'{r}_freq', 0)*100 for r in rds], marker='s', linewidth=2, markersize=10, color='#A23B72')
        a2.set_xlabel('Round'); a2.set_ylabel('Frequency (%)'); a2.set_title('Top Enriched: Frequency'); a2.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'top_sequence_detail.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # Diversity
    if not diversity_df.empty:
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5))
        rl = diversity_df['Round'].str.replace('Round', 'R')
        a1.bar(rl, diversity_df['Shannon_diversity'], color='#4ECDC4', edgecolor='black')
        a1.set_ylabel('Shannon Diversity (bits)'); a1.set_title('Library Diversity'); a1.grid(True, alpha=0.3, axis='y')
        a2.bar(rl, diversity_df['Top1_dominance_%'], color='#FF6B6B', edgecolor='black', label='Top 1')
        a2.bar(rl, diversity_df['Top10_dominance_%'], color='#FFE66D', edgecolor='black', alpha=0.5, label='Top 10')
        a2.set_ylabel('Dominance (%)'); a2.set_title('Sequence Dominance'); a2.legend(); a2.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'diversity_across_rounds.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # CDR3 lengths
    if not annotations.empty:
        cdr3_lens = annotations[annotations['CDR3_length'] > 0]['CDR3_length']
        if len(cdr3_lens) > 0:
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.hist(cdr3_lens, bins=range(0, int(cdr3_lens.max()) + 2), color='#45B7D1', edgecolor='black', alpha=0.8)
            ax.set_xlabel('CDR3 Length (aa)'); ax.set_ylabel('Count'); ax.set_title('CDR3 Length Distribution')
            ax.grid(True, alpha=0.3, axis='y')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'cdr3_length_distribution.png'), dpi=300, bbox_inches='tight')
            plt.close()

    # Plot 5: Pseudo-PPL vs Enrichment
    if ablang_df is not None and not ablang_df.empty and len(top_df) > 0:
        round_names_sorted = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
        first_r, last_r = round_names_sorted[0], round_names_sorted[-1]
        fold_col = f'{first_r}_to_{last_r}_fold'
        if fold_col in top_df.columns:
            merged = top_df.merge(ablang_df, on='Sequence', how='inner')
            if not merged.empty:
                # Add clone family info from annotations
                fam_map = dict(zip(annotations['Sequence'], annotations['clone_family']))
                merged['clone_family'] = merged['Sequence'].map(fam_map).fillna('unknown')

                fig, ax = plt.subplots(figsize=(10, 7))
                families = merged['clone_family'].unique()
                cmap = plt.colormaps['tab10'].resampled(max(len(families), 1))
                for idx_f, fam in enumerate(families):
                    sub = merged[merged['clone_family'] == fam]
                    ax.scatter(sub[fold_col], sub['ablang_pseudo_ppl'],
                               c=[cmap(idx_f % 10)], label=fam if len(families) <= 10 else None,
                               s=60, alpha=0.7, edgecolors='black', linewidths=0.5)
                ax.set_xscale('log')
                ax.set_xlabel('Enrichment Fold-Change (log scale)', fontsize=12, fontweight='bold')
                ax.set_ylabel('AbLang Pseudo-PPL (lower = more natural)', fontsize=12, fontweight='bold')
                ax.set_title('Pseudo-Perplexity vs Enrichment', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)
                # Label top 5 points by rank
                for rank_i in range(min(5, len(merged))):
                    row = merged.iloc[rank_i]
                    ax.annotate(f"Rank {rank_i+1}", (row[fold_col], row['ablang_pseudo_ppl']),
                                fontsize=8, fontweight='bold', ha='left',
                                xytext=(5, 5), textcoords='offset points')
                if len(families) <= 10:
                    ax.legend(fontsize=8, loc='best')
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, 'ppl_vs_enrichment.png'), dpi=300, bbox_inches='tight')
                plt.close()

    # Plot 6: Dominant Clone Variant Heatmap
    if not annotations.empty and len(top_df) > 0:
        # Find the dominant clone family (family of rank-1 sequence)
        top_seq = top_df.iloc[0]['Sequence']
        top_ann = annotations[annotations['Sequence'] == top_seq]
        if not top_ann.empty:
            dom_family = top_ann.iloc[0]['clone_family']
            dom_cdr3 = top_ann.iloc[0]['CDR3']
            # Get variants: same CDR3, different full sequences in top_df
            family_seqs = annotations[annotations['clone_family'] == dom_family]['Sequence'].tolist()
            variants = top_df[top_df['Sequence'].isin(family_seqs)].head(20)
            if len(variants) >= 2:
                seqs = variants['Sequence'].tolist()
                ref = seqs[0]
                max_len = max(len(s) for s in seqs)
                # Build mutation matrix: 0 = match ref, 1 = mismatch, 2 = gap/extension
                mat = np.zeros((len(seqs), max_len), dtype=int)
                for i, seq in enumerate(seqs):
                    for j in range(max_len):
                        if j < len(seq) and j < len(ref):
                            mat[i, j] = 0 if seq[j] == ref[j] else 1
                        else:
                            mat[i, j] = 2
                # Find variable positions (any mismatch across variants)
                var_cols = [j for j in range(max_len) if np.any(mat[:, j] != 0)]
                if var_cols:
                    sub_mat = mat[:, var_cols]
                    # Build text labels for the heatmap cells
                    cell_text = []
                    for i, seq in enumerate(seqs):
                        row_text = []
                        for j in var_cols:
                            row_text.append(seq[j] if j < len(seq) else '-')
                        cell_text.append(row_text)

                    fig, ax = plt.subplots(figsize=(max(6, len(var_cols) * 0.5 + 2), max(4, len(seqs) * 0.5 + 2)))
                    from matplotlib.colors import ListedColormap
                    cmap_hm = ListedColormap(['#E8F5E9', '#FF8A80', '#BDBDBD'])
                    ax.imshow(sub_mat, cmap=cmap_hm, aspect='auto', interpolation='nearest')
                    # Annotate cells with AA letters
                    for i in range(len(seqs)):
                        for j in range(len(var_cols)):
                            ax.text(j, i, cell_text[i][j], ha='center', va='center', fontsize=8, fontweight='bold')
                    ax.set_xticks(range(len(var_cols)))
                    ax.set_xticklabels([str(c + 1) for c in var_cols], fontsize=7, rotation=90)
                    ax.set_yticks(range(len(seqs)))
                    ylabels = [f"Rank {int(variants.iloc[i].name) + 1}" for i in range(len(seqs))]
                    ax.set_yticklabels(ylabels, fontsize=8)
                    ax.set_xlabel('Position', fontsize=11, fontweight='bold')
                    ax.set_title(f'Variant Mutations vs Rank 1 — {dom_family} (CDR3: {dom_cdr3[:15]}...)',
                                 fontsize=12, fontweight='bold')
                    plt.tight_layout()
                    plt.savefig(os.path.join(output_dir, 'dominant_clone_variants.png'), dpi=300, bbox_inches='tight')
                    plt.close()

    # Plot 7: Clone Family Stacked Area Chart
    if clone_tracking is not None and not clone_tracking.empty:
        round_names_sorted = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
        freq_cols = [f'{r}_freq' for r in round_names_sorted if f'{r}_freq' in clone_tracking.columns]
        if freq_cols:
            # Top 8 families by last-round frequency
            top_families = clone_tracking.head(8)
            other_freqs = []
            for col in freq_cols:
                total_top = top_families[col].sum()
                other_freqs.append(max(0, 1.0 - total_top) if total_top < 1.0 else 0)

            fig, ax = plt.subplots(figsize=(12, 7))
            x_labels = [r.replace('Round', 'R') for r in round_names_sorted if f'{r}_freq' in clone_tracking.columns]
            x = range(len(x_labels))
            colors = plt.colormaps['Set2'].resampled(9)
            stack_data = [top_families.iloc[i][freq_cols].astype(float).values * 100 for i in range(len(top_families))]
            if any(v > 0 for v in other_freqs):
                stack_data.append(np.array(other_freqs) * 100)
            labels = top_families['clone_family'].tolist()
            if any(v > 0 for v in other_freqs):
                labels.append('Other')
            ax.stackplot(x, *stack_data,
                         labels=labels,
                         colors=[colors(i) for i in range(len(stack_data))],
                         alpha=0.85)
            ax.set_xticks(list(x))
            ax.set_xticklabels(x_labels)
            ax.set_xlabel('Selection Round', fontsize=12, fontweight='bold')
            ax.set_ylabel('Frequency (%)', fontsize=12, fontweight='bold')
            ax.set_title('Clone Family Dynamics Across Rounds', fontsize=14, fontweight='bold')
            ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
            ax.grid(True, alpha=0.3, axis='y')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'clone_family_dynamics.png'), dpi=300, bbox_inches='tight')
            plt.close()

    # Plot 8: CDR3 Length vs Enrichment
    if not annotations.empty and len(top_df) > 0:
        round_names_sorted = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
        first_r, last_r = round_names_sorted[0], round_names_sorted[-1]
        fold_col = f'{first_r}_to_{last_r}_fold'
        freq_col = f'{last_r}_freq'
        if fold_col in top_df.columns and freq_col in top_df.columns:
            merged_cdr3 = top_df.merge(
                annotations[['Sequence', 'CDR3_length']],
                on='Sequence', how='inner'
            )
            merged_cdr3 = merged_cdr3[merged_cdr3['CDR3_length'] > 0]
            if not merged_cdr3.empty:
                fig, ax = plt.subplots(figsize=(10, 7))
                sizes = merged_cdr3[freq_col] / merged_cdr3[freq_col].max() * 300 + 20
                sc = ax.scatter(merged_cdr3['CDR3_length'], merged_cdr3[fold_col],
                                s=sizes, c='#2E86AB', alpha=0.6, edgecolors='black', linewidths=0.5)
                ax.set_yscale('log')
                ax.set_xlabel('CDR3 Length (aa)', fontsize=12, fontweight='bold')
                ax.set_ylabel('Enrichment Fold-Change (log scale)', fontsize=12, fontweight='bold')
                ax.set_title('CDR3 Length vs Enrichment', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)
                # Add size legend
                for pct, label in [(0.2, 'Low freq'), (0.6, 'Mid freq'), (1.0, 'High freq')]:
                    ax.scatter([], [], s=pct * 300 + 20, c='#2E86AB', alpha=0.6,
                               edgecolors='black', linewidths=0.5, label=label)
                ax.legend(title=f'{last_r} Frequency', loc='best', fontsize=9)
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, 'cdr3_length_vs_enrichment.png'), dpi=300, bbox_inches='tight')
                plt.close()

    # Plot 9: Read Count Bar Chart with Thresholds
    if all_stats:
        round_names_sorted = sorted(all_stats.keys(), key=lambda x: int(x.replace('Round', '')))
        fig, ax = plt.subplots(figsize=(10, 6))
        rl = [r.replace('Round', 'R') for r in round_names_sorted]
        valid_reads = [all_stats[r]['valid_proteins'] for r in round_names_sorted]
        bars = ax.bar(rl, valid_reads, color='#2E86AB', edgecolor='black', alpha=0.85)
        # Reference thresholds from alpseq
        for thresh, color, label in [(2_000_000, '#e83e8c', '2M ideal'),
                                     (500_000, '#fd7e14', '500K'),
                                     (100_000, '#ffc107', '100K min')]:
            if max(valid_reads) > thresh * 0.3:
                ax.axhline(y=thresh, color=color, linestyle='--', alpha=0.7, label=label)
        ax.set_xlabel('Selection Round', fontsize=12, fontweight='bold')
        ax.set_ylabel('Valid Reads', fontsize=12, fontweight='bold')
        ax.set_title('Read Counts Per Round', fontsize=14, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis='y')
        # Add value labels on bars
        for bar, val in zip(bars, valid_reads):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                    f'{val:,.0f}', ha='center', va='bottom', fontsize=9)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'read_counts_per_round.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # Plot 10: Productivity Funnel
    if all_stats:
        round_names_sorted = sorted(all_stats.keys(), key=lambda x: int(x.replace('Round', '')))
        fig, ax = plt.subplots(figsize=(12, 7))
        x = np.arange(len(round_names_sorted))
        width = 0.15
        stages = [
            ('total_reads', 'Total Reads', '#264653'),
            ('merged', 'Merged', '#2A9D8F'),
        ]
        # Compute pass-through stages
        stage_data = {}
        for rnd in round_names_sorted:
            s = all_stats[rnd]
            stage_data[rnd] = {
                'total_reads': s['total_reads'],
                'merged': s['merged'],
                'qual_pass': s['merged'] - s.get('qual_fail', 0) if s['merged'] > s.get('qual_fail', 0) else s['merged'],
                'no_stop': s['valid_proteins'] + s['length_fail'],
                'valid': s['valid_proteins'],
            }
        labels_stages = [
            ('total_reads', 'Total', '#264653'),
            ('merged', 'Merged', '#2A9D8F'),
            ('qual_pass', 'Qual Pass', '#E9C46A'),
            ('no_stop', 'No Stop', '#F4A261'),
            ('valid', 'Valid', '#E76F51'),
        ]
        for i, (key, label, color) in enumerate(labels_stages):
            vals = [stage_data[r][key] for r in round_names_sorted]
            ax.bar(x + i * width, vals, width, label=label, color=color, edgecolor='black', alpha=0.85)
        ax.set_xticks(x + width * 2)
        ax.set_xticklabels([r.replace('Round', 'R') for r in round_names_sorted])
        ax.set_xlabel('Selection Round', fontsize=12, fontweight='bold')
        ax.set_ylabel('Read Count', fontsize=12, fontweight='bold')
        ax.set_title('Productivity Funnel', fontsize=14, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'productivity_funnel.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # Plot 11: Saturation Curves
    if saturation_df is not None and not saturation_df.empty:
        rounds_in_sat = saturation_df['Round'].unique()
        n_rounds = len(rounds_in_sat)
        ncols = min(2, n_rounds)
        nrows = math.ceil(n_rounds / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5 * nrows), squeeze=False)
        for idx, rnd in enumerate(sorted(rounds_in_sat, key=lambda x: int(x.replace('Round', '')))):
            row_i, col_i = divmod(idx, ncols)
            ax = axes[row_i][col_i]
            rnd_data = saturation_df[saturation_df['Round'] == rnd]
            ax.plot(rnd_data['sample_size'], rnd_data['number_unique'],
                    color='#e83e8c', linewidth=2)
            ax.scatter(rnd_data['sample_size'], rnd_data['number_unique'],
                       color='#e83e8c', s=15, alpha=0.6)
            ax.set_xlabel('Sample Size', fontsize=11)
            ax.set_ylabel('Unique CDR3s', fontsize=11)
            ax.set_title(rnd.replace('Round', 'Round '), fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3)
        # Hide unused subplots
        for idx in range(n_rounds, nrows * ncols):
            row_i, col_i = divmod(idx, ncols)
            axes[row_i][col_i].set_visible(False)
        fig.suptitle('Saturation Curves (Rarefaction)', fontsize=15, fontweight='bold', y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'saturation_curves.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # Plot 12: Log2 FC vs Abundance
    if len(top_df) > 0:
        round_names_sorted = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
        last_r = round_names_sorted[-1]
        first_r = round_names_sorted[0]
        cpm_col = f'{last_r}_cpm'
        log2fc_col = f'{first_r}_to_{last_r}_log2fc'
        if cpm_col in top_df.columns and log2fc_col in top_df.columns:
            plot_df = top_df[top_df[cpm_col] > 0].copy()
            if not plot_df.empty:
                fig, ax = plt.subplots(figsize=(10, 7))
                ax.scatter(plot_df[cpm_col], plot_df[log2fc_col],
                           c='#2E86AB', s=50, alpha=0.6, edgecolors='black', linewidths=0.5)
                ax.axhline(y=1.5, color='#e83e8c', linestyle='--', alpha=0.7, label='log2FC = 1.5')
                ax.set_xscale('log')
                ax.set_xlabel(f'{last_r} CPM (log10 scale)', fontsize=12, fontweight='bold')
                ax.set_ylabel('Log2 Fold-Change', fontsize=12, fontweight='bold')
                ax.set_title('Log2 FC vs Abundance', fontsize=14, fontweight='bold')
                ax.legend(fontsize=10)
                ax.grid(True, alpha=0.3)
                # Label top 5
                for i in range(min(5, len(plot_df))):
                    r = plot_df.iloc[i]
                    ax.annotate(f"Rank {i+1}", (r[cpm_col], r[log2fc_col]),
                                fontsize=8, fontweight='bold', xytext=(5, 5),
                                textcoords='offset points')
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, 'log2fc_vs_abundance.png'), dpi=300, bbox_inches='tight')
                plt.close()

    # Plot 13: Convergence Heatmap
    convergence_df = compute_convergence(rounds_data)
    if not convergence_df.empty:
        round_names_sorted = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
        n_rnd = len(round_names_sorted)
        heatmap = np.zeros((n_rnd, n_rnd))
        rnd_to_idx = {r: i for i, r in enumerate(round_names_sorted)}
        for _, row in convergence_df.iterrows():
            i, j = rnd_to_idx[row['Round_A']], rnd_to_idx[row['Round_B']]
            heatmap[i, j] = row['Jaccard']
            heatmap[j, i] = row['Jaccard']
        np.fill_diagonal(heatmap, 1.0)

        fig, ax = plt.subplots(figsize=(8, 7))
        im = ax.imshow(heatmap, cmap='YlOrRd', vmin=0, vmax=1, aspect='auto')
        rl = [r.replace('Round', 'R') for r in round_names_sorted]
        ax.set_xticks(range(n_rnd))
        ax.set_xticklabels(rl, fontsize=10)
        ax.set_yticks(range(n_rnd))
        ax.set_yticklabels(rl, fontsize=10)
        # Annotate cells
        for i in range(n_rnd):
            for j in range(n_rnd):
                ax.text(j, i, f'{heatmap[i,j]:.3f}', ha='center', va='center',
                        fontsize=9, color='black' if heatmap[i,j] < 0.5 else 'white')
        plt.colorbar(im, ax=ax, label='Jaccard Similarity')
        ax.set_title('Convergence Heatmap (Jaccard Similarity)', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'convergence_heatmap.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # Plot 14: V/D/J Gene Usage Stacked Bars (if IgBLAST)
    if gene_usage_df is not None and not gene_usage_df.empty:
        gene_types = [gt for gt in ['V', 'D', 'J'] if gt in gene_usage_df['gene_type'].values]
        if gene_types:
            fig, axes = plt.subplots(1, len(gene_types), figsize=(6 * len(gene_types), 6), squeeze=False)
            for idx, gt in enumerate(gene_types):
                ax = axes[0][idx]
                gt_data = gene_usage_df[gene_usage_df['gene_type'] == gt].copy()
                # Filter V genes to IGHV3 family with >=5% usage (per alpseq)
                if gt == 'V':
                    gt_data = gt_data[(gt_data['gene'].str.contains('IGHV3', na=False)) |
                                      (gt_data['percentage'] >= 5)]
                gt_data = gt_data.sort_values('percentage', ascending=True)
                colors = plt.colormaps['Set2'].resampled(max(len(gt_data), 1))
                bars = ax.barh(gt_data['gene'], gt_data['percentage'],
                               color=[colors(i) for i in range(len(gt_data))],
                               edgecolor='black', alpha=0.85)
                ax.set_xlabel('Usage (%)', fontsize=11, fontweight='bold')
                ax.set_title(f'{gt} Gene Usage', fontsize=13, fontweight='bold')
                ax.grid(True, alpha=0.3, axis='x')
            fig.suptitle('V/D/J Gene Usage', fontsize=15, fontweight='bold')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'vdj_gene_usage.png'), dpi=300, bbox_inches='tight')
            plt.close()

    # Plot 15: Enriched Cluster PCA
    if cluster_pca is not None and not cluster_pca.empty:
        fig, ax = plt.subplots(figsize=(10, 8))
        sizes = cluster_pca['cpm'] / cluster_pca['cpm'].max() * 300 + 20 if cluster_pca['cpm'].max() > 0 else 50
        sc = ax.scatter(cluster_pca['PC1'], cluster_pca['PC2'],
                        s=sizes, c=cluster_pca['enrichment'], cmap='RdYlBu_r',
                        alpha=0.7, edgecolors='black', linewidths=0.5)
        plt.colorbar(sc, ax=ax, label='Enrichment (log2FC)')
        ax.set_xlabel('PC1', fontsize=12, fontweight='bold')
        ax.set_ylabel('PC2', fontsize=12, fontweight='bold')
        ax.set_title('Enriched Cluster PCA', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        # Label top 10
        for i in range(min(10, len(cluster_pca))):
            r = cluster_pca.iloc[i]
            ax.annotate(f"{r['clone_family']}", (r['PC1'], r['PC2']),
                        fontsize=7, fontweight='bold', xytext=(5, 5),
                        textcoords='offset points')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'cluster_pca.png'), dpi=300, bbox_inches='tight')
        plt.close()

    log.info(f"Plots saved to {output_dir}")


# ===========================================================================
# Report
# ===========================================================================
def write_report(top_df, annotations, clone_tracking, diversity_df, rounds_data, output_dir,
                 ablang_df=None, productivity_df=None, saturation_df=None,
                 diversity_est=None, igblast_df=None):
    round_names = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
    rds = [r for r in round_names if sum(rounds_data[r].values()) > 0]
    last = rds[-1] if rds else 'Round1'

    path = os.path.join(output_dir, 'top_sequences_report.txt')
    with open(path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"NANOBODY NGS ENRICHMENT REPORT — {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")

        for rnd in rds:
            f.write(f"  {rnd}: {sum(rounds_data[rnd].values()):>10,} reads, "
                    f"{len(rounds_data[rnd]):>8,} unique\n")

        if not diversity_df.empty:
            f.write("\nDIVERSITY\n" + "-" * 40 + "\n")
            for _, row in diversity_df.iterrows():
                f.write(f"  {row['Round']}: Shannon={row['Shannon_diversity']:.2f}, "
                        f"Top1={row['Top1_dominance_%']:.1f}%, Top10={row['Top10_dominance_%']:.1f}%\n")

        if productivity_df is not None and not productivity_df.empty:
            f.write(f"\nPRODUCTIVITY\n" + "-" * 40 + "\n")
            for _, row in productivity_df.iterrows():
                f.write(f"  {row['Round']}: {int(row['total_reads']):,} total → "
                        f"{int(row['merged']):,} merged ({row['merge_rate']:.1f}%) → "
                        f"{int(row['valid_proteins']):,} valid ({row['pct_productive']:.1f}%)\n")

        if diversity_est is not None and not diversity_est.empty:
            f.write(f"\nESTIMATED LIBRARY DIVERSITY\n" + "-" * 40 + "\n")
            for _, row in diversity_est.iterrows():
                f.write(f"  {row['Round']}: ~{int(row['approx_diversity']):,} unique CDR3s\n")

        f.write(f"\n{'=' * 80}\nTOP ENRICHED SEQUENCES IN {last.upper()}\n{'=' * 80}\n")
        for i in range(min(20, len(top_df))):
            row = top_df.iloc[i]
            seq = row['Sequence']
            f.write(f"\n{'=' * 80}\nRank {i+1} — {int(row['Length'])} aa\n{'=' * 80}\n")
            f.write(f"Sequence: {seq}\n\n")

            for rnd in rds:
                c = int(row.get(f'{rnd}_count', 0))
                fr = row.get(f'{rnd}_freq', 0)
                f.write(f"  {rnd}: {c:>8,} ({fr*100:7.3f}%)\n")

            fold_col = f'{rds[0]}_to_{last}_fold'
            if fold_col in row.index:
                v = row[fold_col]
                f.write(f"\n  Enrichment: {'de novo' if np.isinf(v) else f'{v:.1f}x'}\n")

            # CPM and log2FC
            cpm_col = f'{last}_cpm'
            log2fc_col = f'{rds[0]}_to_{last}_log2fc'
            if cpm_col in row.index:
                f.write(f"  CPM ({last}): {row[cpm_col]:,.1f}\n")
            if log2fc_col in row.index:
                f.write(f"  Log2FC: {row[log2fc_col]:.2f}\n")

            ann = annotations[annotations['Sequence'] == seq]
            if not ann.empty:
                a = ann.iloc[0]
                f.write(f"\n  CDR Annotation ({a['annotation_quality']}):\n")
                for region in ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4'):
                    f.write(f"    {region:5s}: {a[region]}\n")
                f.write(f"    CDR3 length: {a['CDR3_length']} aa\n")
                f.write(f"    Clone: {a['clone_id']} | Family: {a['clone_family']}\n")
                if a.get('is_potential_chimera', False):
                    f.write(f"    ⚠ CHIMERA: {a['chimera_note']}\n")
                elif a.get('chimera_note', ''):
                    f.write(f"    Note: {a['chimera_note']}\n")

            if ablang_df is not None and not ablang_df.empty:
                ab = ablang_df[ablang_df['Sequence'] == seq]
                if not ab.empty:
                    abr = ab.iloc[0]
                    f.write(f"\n  AbLang: pseudo-PPL={abr['ablang_pseudo_ppl']:.3f}, "
                            f"mean_prob={abr['ablang_mean_prob']:.4f}, "
                            f"min_prob={abr['ablang_min_prob']:.4f}\n")

            if igblast_df is not None and not igblast_df.empty and 'sequence_id' in igblast_df.columns:
                ann_row = annotations[annotations['Sequence'] == seq]
                if not ann_row.empty:
                    clone_id = ann_row.iloc[0].get('clone_id', '')
                    ig_match = igblast_df[igblast_df['sequence_id'].str.contains(clone_id, na=False)]
                    if not ig_match.empty:
                        igr = ig_match.iloc[0]
                        f.write(f"\n  IgBLAST: V={igr.get('v_call', 'N/A')}, "
                                f"D={igr.get('d_call', 'N/A')}, J={igr.get('j_call', 'N/A')}\n")

        if not clone_tracking.empty:
            f.write(f"\n\n{'=' * 80}\nCLONE FAMILY TRACKING\n{'=' * 80}\n")
            for i in range(min(20, len(clone_tracking))):
                row = clone_tracking.iloc[i]
                f.write(f"\n{row['clone_family']} — CDR3: {row['representative_CDR3']}\n")
                for rnd in rds:
                    c = int(row.get(f'{rnd}_count', 0))
                    fr = row.get(f'{rnd}_freq', 0)
                    f.write(f"  {rnd}: {c:>8,} ({fr*100:7.3f}%)\n")
                e = row.get('enrichment', 0)
                f.write(f"  Enrichment: {'de novo' if np.isinf(e) else f'{e:.1f}x'}\n")

    log.info(f"Report: {path}")


# ===========================================================================
# Interactive HTML Dashboard (requires plotly)
# ===========================================================================

DASHBOARD_CSS = """
:root{
  --bg:#0c0e13;--surface:#14171e;--surface2:#1a1e28;--border:#252a36;
  --text:#e2e4e9;--text-dim:#7a8094;--accent:#a78bfa;
  --teal:#06b6d4;--orange:#f97316;--green:#34d399;--red:#f87171;
  --mono:'IBM Plex Mono',monospace;
}
[data-theme=light]{
  --bg:#f5f5f7;--surface:#ffffff;--surface2:#eeeef0;--border:#d4d4d8;
  --text:#1a1a2e;--text-dim:#6b7280;--accent:#7c3aed;
  --teal:#0891b2;--orange:#ea580c;--green:#059669;--red:#dc2626;
}
*{margin:0;padding:0;box-sizing:border-box}
*::-webkit-scrollbar{width:8px;height:8px}
*::-webkit-scrollbar-track{background:var(--surface);border-radius:4px}
*::-webkit-scrollbar-thumb{background:var(--border);border-radius:4px}
*::-webkit-scrollbar-thumb:hover{background:var(--text-dim)}
*::-webkit-scrollbar-corner{background:var(--surface)}
*{scrollbar-width:thin;scrollbar-color:var(--border) var(--surface)}
input[type=checkbox]{-webkit-appearance:none;appearance:none;width:16px;height:16px;border:1.5px solid var(--border);
  border-radius:4px;background:var(--surface2);cursor:pointer;position:relative;vertical-align:middle;
  flex-shrink:0;transition:all .15s}
input[type=checkbox]:hover{border-color:var(--accent)}
input[type=checkbox]:checked{background:var(--accent);border-color:var(--accent)}
input[type=checkbox]:checked::after{content:'';position:absolute;left:4.5px;top:1.5px;width:4px;height:8px;
  border:solid #fff;border-width:0 2px 2px 0;transform:rotate(45deg)}
select{-webkit-appearance:none;appearance:none;background:var(--surface2);border:1px solid var(--border);
  color:var(--text);padding:6px 28px 6px 10px;border-radius:6px;font-size:12px;font-family:var(--mono);
  cursor:pointer;transition:all .15s;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%237a8094'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:right 8px center}
select:hover{border-color:var(--accent)}
select:focus{outline:none;border-color:var(--accent)}
select option{background:var(--surface2);color:var(--text)}
body{background:var(--bg);color:var(--text);font-family:'DM Sans',sans-serif;line-height:1.5;min-height:100vh}
.container{max-width:1280px;margin:0 auto;padding:32px 24px}
.header{margin-bottom:32px;display:flex;align-items:center;justify-content:space-between}
.header-left{display:flex;align-items:flex-start;flex-direction:column;gap:4px}
.logo{font-size:28px;font-weight:700;letter-spacing:-0.5px;color:var(--text);display:flex;align-items:center}
.logo span{color:var(--accent)}
.subtitle{font-size:14px;color:var(--text-dim);font-family:var(--mono)}
.header-right{display:flex;align-items:center;gap:16px}
.timestamp{font-size:11px;color:var(--text-dim);font-family:var(--mono)}
.theme-btn{background:var(--surface);border:1px solid var(--border);color:var(--text-dim);
  padding:8px 16px;border-radius:8px;cursor:pointer;font-size:12px;font-family:var(--mono);
  font-weight:500;transition:all .2s}
.theme-btn:hover{border-color:var(--accent);color:var(--text)}
.tabs{display:flex;gap:4px;margin-bottom:32px;background:var(--surface);border-radius:10px;
  padding:4px;width:fit-content}
.tab{padding:10px 24px;border-radius:8px;font-size:14px;font-weight:600;cursor:pointer;
  border:none;background:transparent;color:var(--text-dim);font-family:var(--mono);transition:all .2s}
.tab:hover{color:var(--text)}
.tab.active{background:var(--surface2);color:var(--text);box-shadow:0 1px 3px rgba(0,0,0,0.3)}
.tab-content{display:none}
.tab-content.active{display:block}
.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:16px;margin-bottom:28px}
.metric-card{background:var(--surface);border:1px solid var(--border);border-radius:10px;padding:20px}
.metric-label{font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:1.2px;
  color:var(--text-dim);margin-bottom:8px;font-family:var(--mono)}
.metric-value{font-size:24px;font-weight:700;font-family:var(--mono);color:var(--accent)}
.metric-value.teal{color:var(--teal)}
.metric-value.orange{color:var(--orange)}
.chart-card{background:var(--surface);border:1px solid var(--border);border-radius:10px;
  padding:20px;margin-bottom:16px}
.chart-title{font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:1.2px;
  color:var(--text-dim);margin-bottom:16px;font-family:var(--mono)}
.chart-legend{margin:-8px 0 12px}
.chart-legend summary{font-size:11px;color:var(--accent);cursor:pointer;font-family:var(--mono);
  user-select:none;list-style:none;display:inline-flex;align-items:center;gap:4px}
.chart-legend summary::before{content:'▸';transition:transform 0.15s}
.chart-legend[open] summary::before{transform:rotate(90deg)}
.chart-legend summary::-webkit-details-marker{display:none}
.chart-legend .legend-body{font-size:12px;color:var(--text-dim);line-height:1.6;
  padding:8px 0 4px;font-family:var(--mono)}
.chart-legend .legend-body b{color:var(--text);font-weight:600}
.chart-row{display:grid;grid-template-columns:1fr 1fr;gap:16px}
.placeholder{background:var(--surface);border-radius:10px;padding:40px;text-align:center;
  border:1px dashed var(--border);margin-bottom:16px}
.placeholder-icon{font-size:28px;margin-bottom:8px;opacity:.4}
.placeholder-text{color:var(--text-dim);font-size:12px;font-family:var(--mono)}
table.enrich-table{width:100%;border-collapse:collapse;font-size:12px;font-family:var(--mono)}
table.enrich-table th{background:var(--surface2);color:var(--text-dim);padding:10px 8px;
  text-align:left;cursor:pointer;user-select:none;position:sticky;top:0;
  border-bottom:1px solid var(--border);white-space:nowrap;font-size:11px;font-weight:600;
  text-transform:uppercase;letter-spacing:0.5px}
table.enrich-table th:hover{color:var(--text)}
table.enrich-table td{padding:8px;border-bottom:1px solid var(--border);white-space:nowrap;color:var(--text-dim)}
table.enrich-table tr:hover td{background:var(--surface2);color:var(--text)}
.table-controls{display:flex;gap:12px;align-items:center;margin-bottom:12px}
.table-search,.phylo-search{background:var(--surface2);border:1px solid var(--border);color:var(--text);
  padding:8px 14px;border-radius:8px;font-size:12px;font-family:var(--mono);width:280px}
.table-search:focus,.phylo-search:focus{outline:none;border-color:var(--accent)}
.sort-arrow{font-size:10px;margin-left:4px;opacity:.3}
th.sorted .sort-arrow{opacity:1;color:var(--accent)}
.download-btn{display:inline-flex;align-items:center;gap:8px;background:var(--surface);
  border:1px solid var(--border);color:var(--text-dim);padding:10px 20px;border-radius:8px;
  cursor:pointer;font-size:12px;font-weight:500;font-family:var(--mono);
  text-decoration:none;transition:all .2s;margin-right:8px}
.download-btn:hover{border-color:var(--accent);color:var(--accent)}
.download-btn svg{width:14px;height:14px;fill:currentColor}
.download-bar{display:flex;align-items:center;gap:8px;margin-bottom:20px;flex-wrap:wrap}
.btn-export{background:var(--surface2);border:1px solid var(--border);color:var(--text-dim);
  padding:8px 14px;border-radius:8px;cursor:pointer;font-size:12px;font-family:var(--mono);transition:all .2s}
.btn-export:hover{border-color:var(--accent);color:var(--accent)}
.filter-bar{display:flex;align-items:center;gap:16px;margin-bottom:12px;flex-wrap:wrap;
  padding:10px 16px;background:var(--surface);border-radius:8px;border:1px solid var(--border)}
.filter-bar label{font-size:12px;color:var(--text-dim);cursor:pointer;display:flex;
  align-items:center;gap:4px;font-family:var(--mono)}
.filter-bar input[type=checkbox]{margin-right:2px}
.family-group td{background:var(--accent)!important;color:#fff!important;font-weight:600;
  font-size:11px;letter-spacing:.5px;padding:6px 8px!important}
.pca-linked-row td{background:rgba(167,139,250,0.1)!important}
table.enrich-table tr.pca-highlight td{background:rgba(167,139,250,0.15)!important;color:var(--text)!important}
table.enrich-table tr.cluster-highlight td{background:rgba(244,165,138,0.15)!important;color:var(--text)!important}
.mermaid-wrap{background:transparent;border:none;border-radius:0;
  padding:24px;margin-bottom:16px;text-align:center}
.mermaid-wrap pre.mermaid{background:transparent;border:none;text-align:center}
.methods-text{color:var(--text-dim);font-size:14px;line-height:1.7;max-width:800px;margin:0 auto 24px}
.methods-text h3{font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:1.2px;
  color:var(--accent);margin-bottom:12px;font-family:var(--mono)}
.methods-text p{color:var(--text-dim)}
@media(max-width:900px){.chart-row{grid-template-columns:1fr}.metrics{grid-template-columns:1fr}}
.phylo-legend{display:flex;align-items:center;justify-content:center;gap:14px;padding:10px 0 4px;
  font-size:12px;color:var(--text-dim);font-family:var(--mono)}
.phylo-legend-label{font-weight:600;color:var(--text);margin-right:4px}
.phylo-legend-item{display:inline-flex;align-items:center;gap:3px}
.phylo-legend svg circle{fill:#f4a58a;stroke:#e07050;stroke-width:1}
.phylo-paginated{margin-top:16px}
.phylo-entries-bar{display:flex;align-items:center;gap:8px;margin-bottom:8px;font-size:12px;color:var(--text-dim);font-family:var(--mono)}
.phylo-page-size{font-size:11px;padding:4px 24px 4px 8px}
.custom-select{position:relative;display:inline-block}
.custom-select-trigger{background:var(--surface2);border:1px solid var(--border);color:var(--text);
  padding:6px 28px 6px 10px;border-radius:6px;font-size:12px;font-family:var(--mono);cursor:pointer;
  transition:all .15s;white-space:nowrap;user-select:none;position:relative}
.custom-select-trigger::after{content:'';position:absolute;right:9px;top:50%;transform:translateY(-50%);
  border-left:4px solid transparent;border-right:4px solid transparent;border-top:5px solid var(--text-dim);
  transition:transform .15s}
.custom-select.open .custom-select-trigger{border-color:var(--accent)}
.custom-select.open .custom-select-trigger::after{transform:translateY(-50%) rotate(180deg)}
.custom-select-trigger:hover{border-color:var(--accent)}
.custom-select-options{display:none;position:absolute;bottom:100%;left:0;margin-bottom:4px;
  background:var(--surface2);border:1px solid var(--border);border-radius:6px;min-width:100%;
  z-index:999;overflow:hidden;box-shadow:0 -4px 16px rgba(0,0,0,0.4)}
.custom-select.open .custom-select-options{display:block}
.custom-select.drop-down .custom-select-options{bottom:auto;top:100%;margin-bottom:0;margin-top:4px;
  box-shadow:0 4px 16px rgba(0,0,0,0.4)}
.custom-select-option{padding:6px 12px;font-size:12px;font-family:var(--mono);color:var(--text-dim);
  cursor:pointer;transition:all .1s;white-space:nowrap}
.custom-select-option:hover{background:var(--accent);color:#fff}
.custom-select-option.selected{color:var(--accent)}
.custom-select-option.selected:hover{color:#fff}
.phylo-page-controls{display:flex;align-items:center;justify-content:space-between;margin-top:10px;flex-wrap:wrap;gap:8px}
.phylo-page-info{font-size:12px;color:var(--text-dim);font-family:var(--mono)}
.phylo-page-nav{display:flex;align-items:center;gap:4px}
.phylo-page-btn{background:var(--surface2);border:1px solid var(--border);color:var(--accent);
  padding:4px 10px;border-radius:6px;cursor:pointer;font-size:11px;font-family:var(--mono);transition:all .2s}
.phylo-page-btn:hover{border-color:var(--accent);background:var(--surface)}
.phylo-page-btn.active{background:var(--accent);color:#fff;border-color:var(--accent)}
.phylo-page-btn[disabled]{opacity:.4;cursor:default}
.phylo-page-ellipsis{color:var(--text-dim);font-size:11px;padding:0 4px}
.btn-export-xl{background:var(--surface2);border:1px solid var(--border);color:var(--text-dim);
  padding:8px 14px;border-radius:8px;cursor:pointer;font-size:12px;font-family:var(--mono);transition:all .2s}
.btn-export-xl:hover{border-color:#34d399;color:#34d399}
"""

DASHBOARD_JS = """
var mermaidRendered=false;
function initTabs(){
  document.querySelectorAll('.tab').forEach(t=>{
    t.addEventListener('click',()=>{
      document.querySelectorAll('.tab').forEach(x=>x.classList.remove('active'));
      document.querySelectorAll('.tab-content').forEach(x=>x.classList.remove('active'));
      t.classList.add('active');
      document.getElementById(t.dataset.tab).classList.add('active');
      window.dispatchEvent(new Event('resize'));
      if(t.dataset.tab==='tab-methods'&&!mermaidRendered)renderMermaid();
    });
  });
}
function initTheme(){
  const btn=document.getElementById('theme-toggle');
  let dark=true;
  btn.addEventListener('click',()=>{
    dark=!dark;
    document.documentElement.setAttribute('data-theme',dark?'dark':'light');
    btn.textContent=dark?'Light Mode':'Dark Mode';
    const tc=dark?'#e2e4e9':'#1a1a2e';
    const dc=dark?'#7a8094':'#6b7280';
    const gc=dark?'rgba(255,255,255,0.05)':'rgba(0,0,0,0.06)';
    document.querySelectorAll('.plotly-graph-div').forEach(d=>{
      Plotly.relayout(d,{'font.color':tc,'xaxis.gridcolor':gc,'yaxis.gridcolor':gc,
        'xaxis.color':dc,'yaxis.color':dc,'legend.font.color':dc});
    });
    if(mermaidRendered){mermaidRendered=false;renderMermaid();}
  });
}
function initTableSort(){
  document.querySelectorAll('table.enrich-table').forEach(table=>{
    const headers=table.querySelectorAll('th');
    headers.forEach((th,col)=>{
      th.addEventListener('click',()=>{
        const tbody=table.querySelector('tbody');
        const rows=Array.from(tbody.querySelectorAll('tr'));
        const asc=!th.classList.contains('sorted-asc');
        headers.forEach(h=>{h.classList.remove('sorted','sorted-asc','sorted-desc')});
        th.classList.add('sorted',asc?'sorted-asc':'sorted-desc');
        rows.sort((a,b)=>{
          let va=a.cells[col].getAttribute('data-val')||a.cells[col].textContent;
          let vb=b.cells[col].getAttribute('data-val')||b.cells[col].textContent;
          let na=parseFloat(va),nb=parseFloat(vb);
          if(!isNaN(na)&&!isNaN(nb))return asc?na-nb:nb-na;
          return asc?va.localeCompare(vb):vb.localeCompare(va);
        });
        rows.forEach(r=>tbody.appendChild(r));
      });
    });
  });
}
function initTableFilter(){
  document.querySelectorAll('.table-search').forEach(input=>{
    const tableId=input.dataset.table;
    const table=document.getElementById(tableId);
    if(!table)return;
    input.addEventListener('input',()=>{
      const q=input.value.toLowerCase();
      table.querySelectorAll('tbody tr').forEach(r=>{
        r.style.display=r.textContent.toLowerCase().includes(q)?'':'none';
      });
    });
  });
}
function initCsvExport(){
  document.querySelectorAll('.btn-export').forEach(btn=>{
    btn.addEventListener('click',()=>{
      const tid=btn.dataset.table;
      const table=document.getElementById(tid);
      if(!table)return;
      const rows=[];
      table.querySelectorAll('thead tr').forEach(tr=>{
        const cells=[];
        tr.querySelectorAll('th').forEach(th=>cells.push('"'+th.textContent.replace(/[▲▼]/g,'').trim()+'"'));
        rows.push(cells.join(','));
      });
      table.querySelectorAll('tbody tr').forEach(tr=>{
        if(tr.classList.contains('family-group'))return;
        if(tr.style.display==='none')return;
        const cells=[];
        tr.querySelectorAll('td').forEach(td=>cells.push('"'+(td.getAttribute('data-val')||td.textContent).trim()+'"'));
        rows.push(cells.join(','));
      });
      const blob=new Blob([rows.join('\\n')],{type:'text/csv'});
      const a=document.createElement('a');
      a.href=URL.createObjectURL(blob);
      a.download=tid+'.csv';
      a.click();
      URL.revokeObjectURL(a.href);
    });
  });
}
function initPcaCrosstalk(){
  const pcaDiv=document.getElementById('chart_pca');
  const pcaTable=document.getElementById('pca-table');
  if(!pcaDiv||!pcaTable)return;
  const defaultColor='#DC7F9B';const defaultEdge='#760a2a';const defaultSize=8;
  const accentColor=getComputedStyle(document.documentElement).getPropertyValue('--accent').trim()||'#a78bfa';
  let pinnedRank=null;
  let baseSizes=null;
  function storeBaseSizes(){
    if(!baseSizes&&pcaDiv.data){
      baseSizes=[];
      for(let i=0;i<pcaDiv.data.length;i++){
        const s=pcaDiv.data[i].marker?pcaDiv.data[i].marker.size:defaultSize;
        baseSizes.push(Array.isArray(s)?s.slice():[s]);
      }
    }
  }
  function clearAll(){
    pcaTable.querySelectorAll('tbody tr.pca-highlight').forEach(r=>r.classList.remove('pca-highlight'));
    storeBaseSizes();
    if(pcaDiv.data){
      for(let i=0;i<pcaDiv.data.length;i++){
        const n=(pcaDiv.data[i].x||[]).length;
        Plotly.restyle(pcaDiv,{'marker.color':[Array(n).fill(defaultColor)],
          'marker.size':[baseSizes?baseSizes[i]:Array(n).fill(defaultSize)],
          'marker.line.color':[Array(n).fill(defaultEdge)]},i);
      }
    }
  }
  function getRank(pt){
    const cd=pt.customdata;
    if(Array.isArray(cd))return cd[0];
    return cd;
  }
  function showHighlightByRank(rank){
    storeBaseSizes();
    pcaTable.querySelectorAll('tbody tr').forEach(r=>{
      if(parseInt(r.dataset.rank)===rank)r.classList.add('pca-highlight');
    });
    if(pcaDiv.data){
      for(let i=0;i<pcaDiv.data.length;i++){
        const cd=pcaDiv.data[i].customdata||[];
        const colors=cd.map(v=>{const r=Array.isArray(v)?v[0]:v;return r===rank?accentColor:defaultColor;});
        const bs=baseSizes?baseSizes[i]:cd.map(()=>defaultSize);
        const sizes=cd.map((v,j)=>{const r=Array.isArray(v)?v[0]:v;return r===rank?Math.max((bs[j]||defaultSize)+4,14):(bs[j]||defaultSize);});
        const ec=cd.map(v=>{const r=Array.isArray(v)?v[0]:v;return r===rank?accentColor:defaultEdge;});
        Plotly.restyle(pcaDiv,{'marker.color':[colors],'marker.size':[sizes],'marker.line.color':[ec]},i);
      }
    }
  }
  // Chart hover → table
  pcaDiv.on('plotly_hover',function(evtData){
    if(pinnedRank!==null)return;
    if(!evtData||!evtData.points||!evtData.points.length)return;
    const rank=getRank(evtData.points[0]);
    if(rank===undefined||rank===null)return;
    clearAll();showHighlightByRank(rank);
    const row=pcaTable.querySelector('tbody tr[data-rank="'+rank+'"]');
    if(row)row.scrollIntoView({behavior:'smooth',block:'nearest'});
  });
  pcaDiv.on('plotly_unhover',function(){
    if(pinnedRank!==null)return;
    clearAll();
  });
  // Chart click → pin
  pcaDiv.on('plotly_click',function(evtData){
    if(!evtData||!evtData.points||!evtData.points.length)return;
    const rank=getRank(evtData.points[0]);
    if(rank===undefined||rank===null)return;
    if(pinnedRank===rank){pinnedRank=null;clearAll();return;}
    pinnedRank=rank;clearAll();showHighlightByRank(rank);
    const row=pcaTable.querySelector('tbody tr[data-rank="'+rank+'"]');
    if(row)row.scrollIntoView({behavior:'smooth',block:'nearest'});
  });
  // Lasso select (keep existing behavior)
  pcaDiv.on('plotly_selected',function(evtData){
    clearAll();pinnedRank=null;
    if(!evtData||!evtData.points)return;
    const idxs=new Set(evtData.points.map(p=>getRank(p)));
    pcaTable.querySelectorAll('tbody tr').forEach(r=>{
      if(idxs.has(parseInt(r.dataset.rank)))r.classList.add('pca-highlight');
    });
  });
  pcaDiv.on('plotly_deselect',function(){clearAll();pinnedRank=null;});
  // Table hover → chart
  pcaTable.querySelector('tbody').addEventListener('mouseenter',function(e){
    const row=e.target.closest('tr');
    if(!row||row.dataset.rank===undefined||pinnedRank!==null)return;
    clearAll();showHighlightByRank(parseInt(row.dataset.rank));
  },true);
  pcaTable.querySelector('tbody').addEventListener('mouseleave',function(e){
    const row=e.target.closest('tr');
    if(!row||pinnedRank!==null)return;
    clearAll();
  },true);
  // Table click → pin
  pcaTable.querySelector('tbody').addEventListener('click',function(e){
    const row=e.target.closest('tr');
    if(!row||row.dataset.rank===undefined)return;
    const rank=parseInt(row.dataset.rank);
    if(pinnedRank===rank){pinnedRank=null;clearAll();return;}
    pinnedRank=rank;clearAll();showHighlightByRank(rank);
  });
}
function initClusterCrosstalk(){
  const treeDiv=document.getElementById('chart_cluster_phylo');
  const clTable=document.getElementById('cluster-enrich-table');
  if(!treeDiv||!clTable)return;
  const lineageData=window.__clusterLineage||{};
  let pinnedCidx=null;
  let baseSizes=null;
  let baseColors=null;
  const accentColor=getComputedStyle(document.documentElement).getPropertyValue('--accent').trim()||'#a78bfa';
  function storeBase(){
    if(treeDiv.data&&treeDiv.data.length>1){
      if(!baseSizes){const s=treeDiv.data[1].marker.size;baseSizes=Array.isArray(s)?s.slice():[];}
      if(!baseColors){const c=treeDiv.data[1].marker.color;baseColors=Array.isArray(c)?c.slice():[];}
    }
  }
  function clearTree(){
    if(treeDiv.data&&treeDiv.data.length>1){
      storeBase();
      const n=treeDiv.data[1].x.length;
      Plotly.restyle(treeDiv,{'marker.color':[baseColors||Array(n).fill('#f4a58a')],
        'marker.size':[baseSizes||treeDiv.data[1].marker.size],
        'marker.line.color':[Array(n).fill('rgba(0,0,0,0.15)')]},1);
    }
    while(treeDiv.data&&treeDiv.data.length>2){
      Plotly.deleteTraces(treeDiv,treeDiv.data.length-1);
    }
  }
  function clearTable(){
    clTable.querySelectorAll('tbody tr.cluster-highlight').forEach(r=>r.classList.remove('cluster-highlight'));
  }
  function showHighlight(cidx){
    storeBase();
    // Highlight table row
    clTable.querySelectorAll('tbody tr').forEach(r=>{
      if(r.dataset.cidx===String(cidx))r.classList.add('cluster-highlight');
    });
    // Highlight tree node
    if(treeDiv.data&&treeDiv.data.length>1){
      const cd=treeDiv.data[1].customdata||[];
      const bc=baseColors||[];
      const colors=cd.map((v,j)=>v===cidx?accentColor:(bc[j]||'#f4a58a'));
      const sizes=(baseSizes||[]).map((s,j)=>cd[j]===cidx?Math.max(s+6,14):s);
      const ec=cd.map(v=>v===cidx?accentColor:'rgba(0,0,0,0.15)');
      Plotly.restyle(treeDiv,{'marker.color':[colors],'marker.size':[sizes],'marker.line.color':[ec]},1);
    }
    // Draw lineage path (polyline segments)
    const segs=lineageData[cidx];
    if(segs&&segs.length){
      const lx=[],ly=[];
      segs.forEach(s=>{s.forEach(p=>{lx.push(p[0]);ly.push(p[1]);});lx.push(null);ly.push(null);});
      Plotly.addTraces(treeDiv,{x:lx,y:ly,mode:'lines',
        line:{color:accentColor,width:3},showlegend:false,hoverinfo:'skip'});
    }
  }
  // Table row hover → tree
  clTable.querySelector('tbody').addEventListener('mouseenter',function(e){
    const row=e.target.closest('tr');
    if(!row||row.dataset.cidx===undefined||pinnedCidx!==null)return;
    clearTree();clearTable();
    showHighlight(parseInt(row.dataset.cidx));
  },true);
  clTable.querySelector('tbody').addEventListener('mouseleave',function(e){
    const row=e.target.closest('tr');
    if(!row||pinnedCidx!==null)return;
    clearTree();clearTable();
  },true);
  // Table row click → pin/unpin
  clTable.querySelector('tbody').addEventListener('click',function(e){
    const row=e.target.closest('tr');
    if(!row||row.dataset.cidx===undefined)return;
    const cidx=parseInt(row.dataset.cidx);
    if(pinnedCidx===cidx){pinnedCidx=null;clearTree();clearTable();return;}
    pinnedCidx=cidx;
    clearTree();clearTable();
    showHighlight(cidx);
  });
  // Tree node hover → table
  treeDiv.on('plotly_hover',function(evtData){
    if(pinnedCidx!==null)return;
    if(!evtData||!evtData.points||!evtData.points.length)return;
    const pt=evtData.points[0];
    if(pt.curveNumber!==1)return;
    const cidx=pt.customdata;
    if(cidx===undefined||cidx===null)return;
    clearTree();clearTable();
    showHighlight(cidx);
    const row=clTable.querySelector('tbody tr[data-cidx="'+cidx+'"]');
    if(row)row.scrollIntoView({behavior:'smooth',block:'nearest'});
  });
  treeDiv.on('plotly_unhover',function(){
    if(pinnedCidx!==null)return;
    clearTree();clearTable();
  });
  // Tree node click → pin/unpin
  treeDiv.on('plotly_click',function(evtData){
    if(!evtData||!evtData.points||!evtData.points.length)return;
    const pt=evtData.points[0];
    if(pt.curveNumber!==1)return;
    const cidx=pt.customdata;
    if(cidx===undefined||cidx===null)return;
    if(pinnedCidx===cidx){pinnedCidx=null;clearTree();clearTable();return;}
    pinnedCidx=cidx;
    clearTree();clearTable();
    showHighlight(cidx);
    const row=clTable.querySelector('tbody tr[data-cidx="'+cidx+'"]');
    if(row)row.scrollIntoView({behavior:'smooth',block:'nearest'});
  });
}
function initRoundPhyloCrosstalk(){
  document.querySelectorAll('[id^="chart_phylo_"]').forEach(treeDiv=>{
    const rl=treeDiv.id.replace('chart_phylo_','');
    const clTable=document.getElementById('phylo-table-'+rl);
    if(!clTable)return;
    const lineageData=window['__phyloLineage_'+rl]||{};
    let pinnedIdx=null;
    let baseSizes=null;
    const defaultColor='#f4a58a';const defaultEdgeColor='#e07050';
    const accentColor=getComputedStyle(document.documentElement).getPropertyValue('--accent').trim()||'#a78bfa';
    function storeBaseSizes(){
      if(!baseSizes&&treeDiv.data&&treeDiv.data.length>1){
        const s=treeDiv.data[1].marker.size;
        baseSizes=Array.isArray(s)?s.slice():[];
      }
    }
    function clearTree(){
      if(treeDiv.data&&treeDiv.data.length>1){
        storeBaseSizes();
        const n=treeDiv.data[1].x.length;
        Plotly.restyle(treeDiv,{'marker.color':[Array(n).fill(defaultColor)],
          'marker.size':[baseSizes||treeDiv.data[1].marker.size],
          'marker.line.color':[Array(n).fill(defaultEdgeColor)]},1);
      }
      while(treeDiv.data&&treeDiv.data.length>2){
        Plotly.deleteTraces(treeDiv,treeDiv.data.length-1);
      }
    }
    function clearTable(){
      clTable.querySelectorAll('tbody tr.cluster-highlight').forEach(r=>r.classList.remove('cluster-highlight'));
    }
    function showHighlight(lidx){
      storeBaseSizes();
      clTable.querySelectorAll('tbody tr').forEach(r=>{
        if(r.dataset.lidx===String(lidx))r.classList.add('cluster-highlight');
      });
      if(treeDiv.data&&treeDiv.data.length>1){
        const cd=treeDiv.data[1].customdata||[];
        const colors=cd.map(v=>v===lidx?accentColor:defaultColor);
        const sizes=(baseSizes||[]).map((s,j)=>cd[j]===lidx?Math.max(s+6,14):s);
        const ec=cd.map(v=>v===lidx?accentColor:defaultEdgeColor);
        Plotly.restyle(treeDiv,{'marker.color':[colors],'marker.size':[sizes],'marker.line.color':[ec]},1);
      }
      const segs=lineageData[lidx];
      if(segs&&segs.length){
        const lx=[],ly=[];
        segs.forEach(s=>{s.forEach(p=>{lx.push(p[0]);ly.push(p[1]);});lx.push(null);ly.push(null);});
        Plotly.addTraces(treeDiv,{x:lx,y:ly,mode:'lines',
          line:{color:accentColor,width:3},showlegend:false,hoverinfo:'skip'});
      }
    }
    // Table hover → tree
    clTable.querySelector('tbody').addEventListener('mouseenter',function(e){
      const row=e.target.closest('tr');
      if(!row||row.dataset.lidx===undefined||pinnedIdx!==null)return;
      clearTree();clearTable();showHighlight(parseInt(row.dataset.lidx));
    },true);
    clTable.querySelector('tbody').addEventListener('mouseleave',function(e){
      const row=e.target.closest('tr');
      if(!row||pinnedIdx!==null)return;
      clearTree();clearTable();
    },true);
    // Table click → pin/unpin
    clTable.querySelector('tbody').addEventListener('click',function(e){
      const row=e.target.closest('tr');
      if(!row||row.dataset.lidx===undefined)return;
      const lidx=parseInt(row.dataset.lidx);
      if(pinnedIdx===lidx){pinnedIdx=null;clearTree();clearTable();return;}
      pinnedIdx=lidx;clearTree();clearTable();showHighlight(lidx);
    });
    // Tree hover → table
    treeDiv.on('plotly_hover',function(evtData){
      if(pinnedIdx!==null)return;
      if(!evtData||!evtData.points||!evtData.points.length)return;
      const pt=evtData.points[0];
      if(pt.curveNumber!==1)return;
      const lidx=pt.customdata;
      if(lidx===undefined||lidx===null)return;
      clearTree();clearTable();showHighlight(lidx);
      const row=clTable.querySelector('tbody tr[data-lidx="'+lidx+'"]');
      if(row)row.scrollIntoView({behavior:'smooth',block:'nearest'});
    });
    treeDiv.on('plotly_unhover',function(){
      if(pinnedIdx!==null)return;
      clearTree();clearTable();
    });
    // Tree click → pin/unpin
    treeDiv.on('plotly_click',function(evtData){
      if(!evtData||!evtData.points||!evtData.points.length)return;
      const pt=evtData.points[0];
      if(pt.curveNumber!==1)return;
      const lidx=pt.customdata;
      if(lidx===undefined||lidx===null)return;
      if(pinnedIdx===lidx){pinnedIdx=null;clearTree();clearTable();return;}
      pinnedIdx=lidx;clearTree();clearTable();showHighlight(lidx);
      const row=clTable.querySelector('tbody tr[data-lidx="'+lidx+'"]');
      if(row)row.scrollIntoView({behavior:'smooth',block:'nearest'});
    });
  });
}
function initPplCrosstalk(){
  const pplDiv=document.getElementById('chart_ppl');
  const pplTable=document.getElementById('ppl-table');
  if(!pplDiv||!pplTable)return;
  const accentColor=getComputedStyle(document.documentElement).getPropertyValue('--accent').trim()||'#a78bfa';
  let pinnedIdx=null;
  let baseColors=null;let baseSizes=null;let baseEdgeColors=null;
  function storeBase(){
    if(!baseColors&&pplDiv.data){
      baseColors=[];baseSizes=[];baseEdgeColors=[];
      for(let i=0;i<pplDiv.data.length;i++){
        const m=pplDiv.data[i].marker||{};
        const c=m.color;baseColors.push(Array.isArray(c)?c.slice():c);
        const s=m.size;baseSizes.push(Array.isArray(s)?s.slice():s);
        const ec=m.line?m.line.color:null;baseEdgeColors.push(Array.isArray(ec)?ec.slice():ec);
      }
    }
  }
  function clearAll(){
    pplTable.querySelectorAll('tbody tr.cluster-highlight').forEach(r=>r.classList.remove('cluster-highlight'));
    storeBase();
    if(pplDiv.data){
      for(let i=0;i<pplDiv.data.length;i++){
        const n=(pplDiv.data[i].x||[]).length;
        const bc=baseColors[i];const bs=baseSizes[i];const bec=baseEdgeColors[i];
        Plotly.restyle(pplDiv,{
          'marker.color':[Array.isArray(bc)?bc:Array(n).fill(bc)],
          'marker.size':[Array.isArray(bs)?bs:Array(n).fill(bs)],
          'marker.line.color':[Array.isArray(bec)?bec:Array(n).fill(bec||'black')]
        },i);
      }
    }
  }
  function showHighlight(pidx){
    storeBase();
    pplTable.querySelectorAll('tbody tr').forEach(r=>{
      if(r.dataset.pidx===String(pidx))r.classList.add('cluster-highlight');
    });
    if(pplDiv.data){
      for(let i=0;i<pplDiv.data.length;i++){
        const cd=pplDiv.data[i].customdata||[];
        const bc=baseColors[i];const bs=baseSizes[i];
        const defC=Array.isArray(bc)?bc[0]:bc;
        const defS=Array.isArray(bs)?bs[0]:(bs||8);
        const colors=cd.map((v,j)=>v===pidx?accentColor:(Array.isArray(bc)?bc[j]:defC));
        const sizes=cd.map((v,j)=>v===pidx?Math.max((Array.isArray(bs)?bs[j]:defS)+4,14):(Array.isArray(bs)?bs[j]:defS));
        const ec=cd.map(v=>v===pidx?accentColor:'black');
        Plotly.restyle(pplDiv,{'marker.color':[colors],'marker.size':[sizes],'marker.line.color':[ec]},i);
      }
    }
  }
  // Chart hover → table
  pplDiv.on('plotly_hover',function(evtData){
    if(pinnedIdx!==null)return;
    if(!evtData||!evtData.points||!evtData.points.length)return;
    const pidx=evtData.points[0].customdata;
    if(pidx===undefined||pidx===null)return;
    clearAll();showHighlight(pidx);
    const row=pplTable.querySelector('tbody tr[data-pidx="'+pidx+'"]');
    if(row)row.scrollIntoView({behavior:'smooth',block:'nearest'});
  });
  pplDiv.on('plotly_unhover',function(){
    if(pinnedIdx!==null)return;
    clearAll();
  });
  // Chart click → pin
  pplDiv.on('plotly_click',function(evtData){
    if(!evtData||!evtData.points||!evtData.points.length)return;
    const pidx=evtData.points[0].customdata;
    if(pidx===undefined||pidx===null)return;
    if(pinnedIdx===pidx){pinnedIdx=null;clearAll();return;}
    pinnedIdx=pidx;clearAll();showHighlight(pidx);
    const row=pplTable.querySelector('tbody tr[data-pidx="'+pidx+'"]');
    if(row)row.scrollIntoView({behavior:'smooth',block:'nearest'});
  });
  // Table hover → chart
  pplTable.querySelector('tbody').addEventListener('mouseenter',function(e){
    const row=e.target.closest('tr');
    if(!row||row.dataset.pidx===undefined||pinnedIdx!==null)return;
    clearAll();showHighlight(parseInt(row.dataset.pidx));
  },true);
  pplTable.querySelector('tbody').addEventListener('mouseleave',function(e){
    const row=e.target.closest('tr');
    if(!row||pinnedIdx!==null)return;
    clearAll();
  },true);
  // Table click → pin
  pplTable.querySelector('tbody').addEventListener('click',function(e){
    const row=e.target.closest('tr');
    if(!row||row.dataset.pidx===undefined)return;
    const pidx=parseInt(row.dataset.pidx);
    if(pinnedIdx===pidx){pinnedIdx=null;clearAll();return;}
    pinnedIdx=pidx;clearAll();showHighlight(pidx);
  });
}
function initPcaFilters(){
  const container=document.getElementById('pca-filters');
  const pcaDiv=document.getElementById('chart_pca');
  if(!container||!pcaDiv)return;
  container.querySelectorAll('input[type=checkbox]').forEach(cb=>{
    cb.addEventListener('change',()=>{
      const checked=new Set();
      container.querySelectorAll('input[type=checkbox]:checked').forEach(c=>checked.add(c.value));
      const showAll=checked.has('all');
      for(let i=0;i<pcaDiv.data.length;i++){
        const name=pcaDiv.data[i].name||'';
        const vis=showAll||checked.has(name)||checked.has('top10')&&pcaDiv.data[i]._isTop10||
                   checked.has('top25')&&pcaDiv.data[i]._isTop25;
        Plotly.restyle(pcaDiv,{visible:vis?true:'legendonly'},i);
      }
    });
  });
}
function initFamilyGroupToggle(){
  const cb=document.getElementById('group-by-family');
  const table=document.getElementById('enrichment-table');
  if(!cb||!table)return;
  cb.addEventListener('change',()=>{
    const tbody=table.querySelector('tbody');
    const rows=Array.from(tbody.querySelectorAll('tr:not(.family-group)'));
    tbody.querySelectorAll('tr.family-group').forEach(r=>r.remove());
    if(cb.checked){
      rows.sort((a,b)=>{
        const fa=a.dataset.family||'';const fb=b.dataset.family||'';
        if(fa!==fb)return fa.localeCompare(fb);
        return parseInt(a.dataset.rank)-parseInt(b.dataset.rank);
      });
      let lastFam='';
      const nCols=table.querySelectorAll('thead th').length;
      rows.forEach(r=>{
        const fam=r.dataset.family||'';
        if(fam!==lastFam){
          const gh=document.createElement('tr');
          gh.className='family-group';
          gh.innerHTML='<td colspan="'+nCols+'">'+fam+'</td>';
          tbody.appendChild(gh);
          lastFam=fam;
        }
        tbody.appendChild(r);
      });
    }else{
      rows.sort((a,b)=>parseInt(a.dataset.rank)-parseInt(b.dataset.rank));
      rows.forEach(r=>tbody.appendChild(r));
    }
  });
}
async function renderMermaid(){
  if(typeof mermaid==='undefined')return;
  const el=document.querySelector('.mermaid');
  if(!el)return;
  if(!el.dataset.src)el.dataset.src=el.textContent;
  const isDark=document.documentElement.getAttribute('data-theme')!=='light';
  mermaid.initialize({startOnLoad:false,theme:isDark?'dark':'default',
    themeVariables:isDark?{primaryColor:'#e83e8c',primaryTextColor:'#e0e0e0',
      lineColor:'#4ECDC4',secondaryColor:'#16213e'}:
      {primaryColor:'#e83e8c',primaryTextColor:'#2d2d2d',lineColor:'#4ECDC4',secondaryColor:'#fbf0ed'}});
  el.removeAttribute('data-processed');
  el.innerHTML=el.dataset.src;
  const {svg}=await mermaid.render('mermaid-svg',el.dataset.src);
  el.innerHTML=svg;
  mermaidRendered=true;
}
function initPhyloPagination(){
  document.querySelectorAll('.phylo-paginated').forEach(wrapper=>{
    const table=wrapper.querySelector('table');
    const info=wrapper.querySelector('.phylo-page-info');
    const nav=wrapper.querySelector('.phylo-page-nav');
    const sizeSelect=wrapper.querySelector('.phylo-page-size');
    const searchInput=wrapper.querySelector('.phylo-search');
    if(!table)return;
    let pageSize=parseInt(sizeSelect?sizeSelect.value:'5');
    let currentPage=1;
    function allRows(){return Array.from(table.querySelectorAll('tbody tr'));}
    function getVisible(){return allRows().filter(r=>r.dataset.filtered!=='true');}
    function render(){
      const rows=getVisible();const total=rows.length;
      const totalPages=Math.max(1,Math.ceil(total/pageSize));
      if(currentPage>totalPages)currentPage=totalPages;
      const start=(currentPage-1)*pageSize;const end=Math.min(start+pageSize,total);
      allRows().forEach(r=>r.style.display='none');
      rows.slice(start,end).forEach(r=>r.style.display='');
      if(info)info.textContent=total>0?'Showing '+(start+1)+' to '+end+' of '+total+' entries':'No entries';
      if(nav){
        let html='<button class="phylo-page-btn" data-page="prev"'+(currentPage<=1?' disabled':'')+'>Previous</button>';
        const pages=new Set([1]);
        for(let p=Math.max(1,currentPage-1);p<=Math.min(totalPages,currentPage+1);p++)pages.add(p);
        pages.add(totalPages);
        const sorted=[...pages].sort((a,b)=>a-b);let last=0;
        sorted.forEach(p=>{
          if(p-last>1)html+='<span class="phylo-page-ellipsis">&hellip;</span>';
          html+='<button class="phylo-page-btn'+(p===currentPage?' active':'')+'" data-page="'+p+'">'+p+'</button>';
          last=p;
        });
        html+='<button class="phylo-page-btn" data-page="next"'+(currentPage>=totalPages?' disabled':'')+'>Next</button>';
        nav.innerHTML=html;
        nav.querySelectorAll('.phylo-page-btn').forEach(btn=>{
          btn.addEventListener('click',()=>{
            const p=btn.dataset.page;
            if(p==='prev')currentPage=Math.max(1,currentPage-1);
            else if(p==='next')currentPage=Math.min(totalPages,currentPage+1);
            else currentPage=parseInt(p);
            render();
          });
        });
      }
    }
    if(searchInput){searchInput.addEventListener('input',()=>{
      const q=searchInput.value.toLowerCase();
      allRows().forEach(r=>{r.dataset.filtered=r.textContent.toLowerCase().includes(q)?'false':'true';});
      currentPage=1;render();
    });}
    if(sizeSelect){sizeSelect.addEventListener('change',()=>{
      pageSize=parseInt(sizeSelect.value);currentPage=1;render();
    });}
    render();
  });
}
function initXlExport(){
  document.querySelectorAll('.btn-export-xl').forEach(btn=>{
    btn.addEventListener('click',()=>{
      const tid=btn.dataset.table;const table=document.getElementById(tid);if(!table)return;
      const html='<html><head><meta charset="utf-8"></head><body>'+table.outerHTML+'</body></html>';
      const blob=new Blob([html],{type:'application/vnd.ms-excel'});
      const a=document.createElement('a');a.href=URL.createObjectURL(blob);
      a.download=tid+'.xls';a.click();URL.revokeObjectURL(a.href);
    });
  });
}
function initCustomSelects(){
  document.querySelectorAll('select').forEach(sel=>{
    if(sel.dataset.customized)return;
    sel.dataset.customized='1';
    sel.style.display='none';
    const wrap=document.createElement('div');
    wrap.className='custom-select';
    sel.parentNode.insertBefore(wrap,sel);
    wrap.appendChild(sel);
    const trigger=document.createElement('div');
    trigger.className='custom-select-trigger';
    const opts=document.createElement('div');
    opts.className='custom-select-options';
    function buildOpts(){
      opts.innerHTML='';
      Array.from(sel.options).forEach(o=>{
        const d=document.createElement('div');
        d.className='custom-select-option'+(o.selected?' selected':'');
        d.textContent=o.textContent;d.dataset.value=o.value;
        d.addEventListener('click',e=>{
          e.stopPropagation();
          sel.value=o.value;sel.dispatchEvent(new Event('change',{bubbles:true}));
          opts.querySelectorAll('.custom-select-option').forEach(x=>x.classList.remove('selected'));
          d.classList.add('selected');
          trigger.textContent=o.textContent;
          wrap.classList.remove('open');
        });
        opts.appendChild(d);
      });
      const cur=sel.options[sel.selectedIndex];
      trigger.textContent=cur?cur.textContent:'';
    }
    buildOpts();
    wrap.appendChild(trigger);wrap.appendChild(opts);
    trigger.addEventListener('click',e=>{
      e.stopPropagation();
      document.querySelectorAll('.custom-select.open').forEach(w=>{if(w!==wrap)w.classList.remove('open');});
      wrap.classList.toggle('open');
      if(wrap.classList.contains('open')){
        const r=wrap.getBoundingClientRect();
        const spaceAbove=r.top;const spaceBelow=window.innerHeight-r.bottom;
        wrap.classList.toggle('drop-down',spaceBelow>spaceAbove||spaceBelow>200);
      }
    });
    new MutationObserver(buildOpts).observe(sel,{childList:true,attributes:true,subtree:true});
  });
  document.addEventListener('click',()=>{document.querySelectorAll('.custom-select.open').forEach(w=>w.classList.remove('open'));});
}
document.addEventListener('DOMContentLoaded',()=>{initTabs();initTheme();initTableSort();initTableFilter();
  initCsvExport();initPcaCrosstalk();initPcaFilters();initFamilyGroupToggle();
  initCustomSelects();initPhyloPagination();initXlExport();initClusterCrosstalk();initRoundPhyloCrosstalk();initPplCrosstalk();});
"""

DASHBOARD_PALETTE = ['#a78bfa','#06b6d4','#f97316','#34d399','#f87171',
                      '#fbbf24','#ec4899','#8b5cf6','#14b8a6','#fb923c']


def generate_dashboard(top_df, diversity_df, rounds_data, annotations, output_dir,
                       all_stats, clone_tracking, saturation_df, diversity_est,
                       productivity_df, convergence_df, ablang_df=None,
                       igblast_df=None, gene_usage_df=None, cluster_pca=None,
                       fasta_enrichment='', fasta_diversity=''):
    """Generate an interactive HTML dashboard. Requires plotly."""
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        log.warning("plotly not installed — skipping dashboard (pip install plotly)")
        return None

    round_names = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))
    rds = [r for r in round_names if sum(rounds_data[r].values()) > 0]
    labels = [r.replace('Round', 'R') for r in rds]
    first_r, last_r = round_names[0], round_names[-1]
    fold_col = f'{first_r}_to_{last_r}_fold'
    log2fc_col = f'{first_r}_to_{last_r}_log2fc'
    last_cpm_col = f'{last_r}_cpm'

    def _base_layout(**kw):
        layout = dict(
            font=dict(family='DM Sans, sans-serif', color='#e2e4e9', size=12),
            paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
            margin=dict(l=50, r=30, t=40, b=40),
            xaxis=dict(gridcolor='rgba(255,255,255,0.05)', color='#7a8094',
                       tickfont=dict(family='IBM Plex Mono, monospace', size=11)),
            yaxis=dict(gridcolor='rgba(255,255,255,0.05)', color='#7a8094',
                       tickfont=dict(family='IBM Plex Mono, monospace', size=11)),
            legend=dict(font=dict(size=11, color='#7a8094', family='IBM Plex Mono, monospace')),
        )
        layout.update(kw)
        return layout

    def _to_div(fig, uid=None):
        div_id = uid or f"chart_{id(fig)}"
        return fig.to_html(full_html=False, include_plotlyjs=False, div_id=div_id)

    def _legend(text):
        return (f'<details class="chart-legend"><summary>Figure legend</summary>'
                f'<div class="legend-body">{text}</div></details>')

    def _placeholder(icon, text, hint=''):
        return (f'<div class="placeholder"><div class="placeholder-icon">{icon}</div>'
                f'<div class="placeholder-text">{text}</div>'
                f'{"<div class=\"placeholder-text\" style=\"margin-top:6px;font-size:11px\">" + hint + "</div>" if hint else ""}'
                f'</div>')

    def _radial_layout(tree):
        """Equal-angle radial layout for a phylogenetic tree.
        Returns (positions, edges, lineages) where positions maps clade key to (x, y),
        edges is list of [(x,y), ...] polyline point lists (arc + radial elbow edges),
        and lineages maps each leaf key to the list of edge indices forming its
        root-to-leaf path."""
        import math

        def count_terminals(clade):
            if clade.is_terminal():
                return 1
            return sum(count_terminals(c) for c in clade.clades)

        positions = {}
        polar = {}  # clade_key -> (angle, radius)
        edges = []
        # Track parent relationship for lineage reconstruction
        parent_map = {}  # clade_key -> parent_clade_key
        edge_lookup = {}  # (parent_key, child_key) -> edge_index

        def _key(clade):
            return clade.name if clade.name else f"_i{id(clade)}"

        def layout(clade, angle_start, angle_end, radius):
            angle_mid = (angle_start + angle_end) / 2
            x = math.cos(angle_mid) * radius
            y = math.sin(angle_mid) * radius
            positions[_key(clade)] = (x, y)
            polar[_key(clade)] = (angle_mid, radius)
            if clade.is_terminal():
                return
            n_leaves = count_terminals(clade)
            cur = angle_start
            for child in clade.clades:
                child_leaves = count_terminals(child)
                child_span = (angle_end - angle_start) * child_leaves / n_leaves
                bl = child.branch_length if child.branch_length and child.branch_length > 0 else 0.01
                child_radius = radius + bl
                layout(child, cur, cur + child_span, child_radius)
                # Build elbow edge: arc at parent radius, then radial to child
                p_angle = angle_mid
                c_angle, _ = polar[_key(child)]
                pts = []
                # Arc from parent angle to child angle at parent radius
                if radius > 1e-9:
                    n_arc = 20
                    for k in range(n_arc + 1):
                        t = k / n_arc
                        a = p_angle + t * (c_angle - p_angle)
                        pts.append((math.cos(a) * radius, math.sin(a) * radius))
                else:
                    pts.append((x, y))
                # Radial line to child position
                pts.append(positions[_key(child)])
                edge_idx = len(edges)
                edges.append(pts)
                parent_map[_key(child)] = _key(clade)
                edge_lookup[(_key(clade), _key(child))] = edge_idx
                cur += child_span

        layout(tree.root, 0, 2 * math.pi, 0)

        # Build lineages: for each leaf, collect edge indices from root to leaf
        lineages = {}
        root_key = _key(tree.root)
        for clade_key in positions:
            # Only include terminal nodes (leaves have names like "c0", "c1", etc.)
            if clade_key.startswith('_i'):
                continue  # internal node
            # Walk from leaf up to root, collecting edge indices
            path_edges = []
            cur_key = clade_key
            while cur_key in parent_map:
                p_key = parent_map[cur_key]
                eidx = edge_lookup.get((p_key, cur_key))
                if eidx is not None:
                    path_edges.append(eidx)
                cur_key = p_key
            path_edges.reverse()  # root-to-leaf order
            lineages[clade_key] = path_edges

        return positions, edges, lineages

    def _build_round_phylo(sequences, counts, total_rnd, seq_ann, round_label):
        """Build a radial NJ phylogenetic tree figure for one round.
        Returns (fig, lineage_json) or (None, None)."""
        try:
            from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
            from Bio.Align import MultipleSeqAlignment
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
        except ImportError:
            return None, None

        if len(sequences) < 3:
            return None, None

        max_len = max(len(s) for s in sequences)
        records = []
        id_to_seq = {}
        for i, seq in enumerate(sequences):
            padded = seq + '-' * (max_len - len(seq))
            sid = f"s{i}"
            id_to_seq[sid] = seq
            records.append(SeqRecord(Seq(padded), id=sid, description=''))

        aln = MultipleSeqAlignment(records)
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

        positions, edges, lineages = _radial_layout(tree)

        fig = go.Figure()

        # Single edge trace with None separators (elbow polylines)
        edge_x, edge_y = [], []
        for pts in edges:
            for (px, py) in pts:
                edge_x.append(px)
                edge_y.append(py)
            edge_x.append(None)
            edge_y.append(None)
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y, mode='lines',
            line=dict(color='#9ca3af', width=1),
            showlegend=False, hoverinfo='skip'))

        # Leaf nodes: size ∝ log10(CPM), salmon color
        leaf_x, leaf_y, leaf_sizes, leaf_hover, leaf_customdata = [], [], [], [], []
        for i, (sid, seq) in enumerate(id_to_seq.items()):
            if sid not in positions:
                continue
            x, y = positions[sid]
            leaf_x.append(x)
            leaf_y.append(y)
            leaf_customdata.append(i)
            cnt = counts[seq]
            cpm = cnt / total_rnd * 1e6
            ann = seq_ann.get(seq, {})
            cdr3 = (ann.get('CDR3', '') if isinstance(ann, dict)
                    else (ann['CDR3'] if 'CDR3' in ann.index else ''))
            family = (ann.get('clone_family', '') if isinstance(ann, dict)
                      else (ann['clone_family'] if 'clone_family' in ann.index else ''))
            log_cpm = np.log10(max(cpm, 1))
            size = max(4, min(35, (log_cpm - 2.5) * 12))
            leaf_sizes.append(size)
            leaf_hover.append(
                f'CDR3: {cdr3}<br>Family: {family}<br>'
                f'CPM: {cpm:,.0f}<br>Count: {cnt:,}')

        fig.add_trace(go.Scatter(
            x=leaf_x, y=leaf_y, mode='markers',
            marker=dict(size=leaf_sizes, color='#f4a58a', opacity=0.85,
                        line=dict(color='#e07050', width=1)),
            customdata=leaf_customdata,
            text=leaf_hover, hoverinfo='text', showlegend=False))

        fig.update_layout(**_base_layout(
            xaxis=dict(visible=False), yaxis=dict(visible=False, scaleanchor='x'),
            height=520, margin=dict(l=10, r=10, t=10, b=10)))

        # Build lineage edge coordinate dict for JS crosstalk (polyline points)
        import json as _json
        lineage_edges = {}
        for i, sid in enumerate(id_to_seq.keys()):
            edge_idxs = lineages.get(sid, [])
            if edge_idxs:
                coords = []
                for eidx in edge_idxs:
                    coords.append([[round(px, 6), round(py, 6)] for px, py in edges[eidx]])
                lineage_edges[i] = coords
        lineage_json = _json.dumps(lineage_edges)

        return fig, lineage_json

    def _log2fc_color(val, vmin, vmax):
        """Map a log2fc value to a color on a 3-stop gradient:
        #2d0a1e (dark) → #8b2252 (mid) → #e83e8c (bright pink)."""
        if vmax <= vmin:
            return 'rgb(139,34,82)'
        t = max(0.0, min(1.0, (val - vmin) / (vmax - vmin)))
        stops = [(0x2d, 0x0a, 0x1e), (0x8b, 0x22, 0x52), (0xe8, 0x3e, 0x8c)]
        if t <= 0.5:
            s = t * 2
            r = int(stops[0][0] + s * (stops[1][0] - stops[0][0]))
            g = int(stops[0][1] + s * (stops[1][1] - stops[0][1]))
            b = int(stops[0][2] + s * (stops[1][2] - stops[0][2]))
        else:
            s = (t - 0.5) * 2
            r = int(stops[1][0] + s * (stops[2][0] - stops[1][0]))
            g = int(stops[1][1] + s * (stops[2][1] - stops[1][1]))
            b = int(stops[1][2] + s * (stops[2][2] - stops[1][2]))
        return f'rgb({r},{g},{b})'

    def _build_cluster_phylo(cluster_df, annotations_df, rounds_data_dict, round_names_list, max_clusters=100):
        """Build a radial NJ tree of top enriched cluster leads using Gonnet matrix.
        Returns (fig, matrix_name, lineage_json, family_to_idx, l2fc_range) or
        (None, None, None, None, None)."""
        try:
            from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
            from Bio.Align import MultipleSeqAlignment
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
        except ImportError:
            return None, None, None, None, None

        if cluster_df is None or cluster_df.empty:
            return None, None, None, None, None

        # Get top clusters sorted by last-round CPM
        last_rnd = round_names_list[-1]
        last_cpm_key = f'{last_rnd}_cpm'
        top_cl = cluster_df.head(max_clusters).copy()
        if len(top_cl) < 3:
            return None, None, None, None, None

        # Map family → representative full protein sequence (prefer full seq over CDR3)
        fam_to_seq = {}
        if not annotations_df.empty:
            for fam in top_cl['clone_family']:
                fam_rows = annotations_df[annotations_df['clone_family'] == fam]
                if not fam_rows.empty:
                    fam_to_seq[fam] = fam_rows.iloc[0]['Sequence']

        # Build sequence list; fall back to CDR3 if no full sequence
        seq_list = []
        meta = []  # (family, cdr3, last_cpm, cluster_size, log2fc)
        family_sizes = annotations_df.groupby('clone_family').size().to_dict() if not annotations_df.empty else {}
        for _, row in top_cl.iterrows():
            fam = row['clone_family']
            cdr3 = row.get('representative_CDR3', '')
            seq = fam_to_seq.get(fam, cdr3)
            if not seq or len(seq) < 3:
                continue
            cpm = row.get(last_cpm_key, 0)
            size = family_sizes.get(fam, 1)
            l2fc = float(row.get('log2fc', 0))
            seq_list.append(seq)
            meta.append((fam, cdr3, float(cpm), int(size), l2fc))

        if len(seq_list) < 3:
            return None, None, None, None, None

        # Pad to equal length and build MSA
        max_len = max(len(s) for s in seq_list)
        records = []
        for i, seq in enumerate(seq_list):
            padded = seq + '-' * (max_len - len(seq))
            records.append(SeqRecord(Seq(padded), id=f"c{i}", description=''))

        aln = MultipleSeqAlignment(records)

        # Try Gonnet, fall back to blosum62, then identity
        matrix_used = 'Gonnet'
        try:
            dm = DistanceCalculator('gonnet1992').get_distance(aln)
        except Exception:
            matrix_used = 'BLOSUM62'
            try:
                dm = DistanceCalculator('blosum62').get_distance(aln)
            except Exception:
                matrix_used = 'identity'
                dm = DistanceCalculator('identity').get_distance(aln)

        tree = DistanceTreeConstructor().nj(dm)
        positions, edges, lineages = _radial_layout(tree)

        fig = go.Figure()

        # Edges (elbow polylines)
        edge_x, edge_y = [], []
        for pts in edges:
            for (px, py) in pts:
                edge_x.append(px)
                edge_y.append(py)
            edge_x.append(None)
            edge_y.append(None)
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y, mode='lines',
            line=dict(color='#9ca3af', width=1),
            showlegend=False, hoverinfo='skip'))

        # Compute log2fc color range across all leaves
        all_l2fc = [m[4] for m in meta]
        l2fc_min, l2fc_max = min(all_l2fc), max(all_l2fc)

        # Leaf nodes — track customdata (leaf index) for JS crosstalk
        leaf_x, leaf_y, leaf_sizes, leaf_hover, leaf_customdata = [], [], [], [], []
        leaf_colors = []
        leaf_idx_map = {}  # maps sequential leaf position → meta index i
        for i, (fam, cdr3, cpm, size, l2fc) in enumerate(meta):
            sid = f"c{i}"
            if sid not in positions:
                continue
            x, y = positions[sid]
            leaf_idx_map[len(leaf_x)] = i
            leaf_x.append(x)
            leaf_y.append(y)
            log_cpm = np.log10(max(cpm, 1))
            sz = max(4, min(35, (log_cpm - 2.5) * 12))
            leaf_sizes.append(sz)
            leaf_customdata.append(i)
            leaf_colors.append(_log2fc_color(l2fc, l2fc_min, l2fc_max))
            leaf_hover.append(
                f'CDR3: {cdr3}<br>Family: {fam}<br>'
                f'CPM: {cpm:,.0f}<br>Cluster size: {size}<br>'
                f'Log2FC: {l2fc:.2f}')

        fig.add_trace(go.Scatter(
            x=leaf_x, y=leaf_y, mode='markers',
            marker=dict(size=leaf_sizes, color=leaf_colors, opacity=0.85,
                        line=dict(color='rgba(0,0,0,0.15)', width=1)),
            customdata=leaf_customdata,
            text=leaf_hover, hoverinfo='text', showlegend=False))

        fig.update_layout(**_base_layout(
            xaxis=dict(visible=False), yaxis=dict(visible=False, scaleanchor='x'),
            height=560, margin=dict(l=10, r=10, t=10, b=10)))

        # Build lineage edge coordinate dict for JS crosstalk (polyline points)
        import json as _json
        lineage_edges = {}
        for i, (fam, cdr3, cpm, size, l2fc) in enumerate(meta):
            sid = f"c{i}"
            edge_idxs = lineages.get(sid, [])
            if edge_idxs:
                coords = []
                for eidx in edge_idxs:
                    coords.append([[round(px, 6), round(py, 6)] for px, py in edges[eidx]])
                lineage_edges[i] = coords
        lineage_json = _json.dumps(lineage_edges)

        # Map family name → meta index for table row tagging
        family_to_idx = {fam: i for i, (fam, _, _, _, _) in enumerate(meta)}

        return fig, matrix_used, lineage_json, family_to_idx, (l2fc_min, l2fc_max)

    # ── Tab 1: Overview ──────────────────────────────────────────────────

    # Metric cards
    total_rounds = len(rds)
    total_valid = sum(all_stats[r]['valid_proteins'] for r in rds)
    all_unique = set()
    for r in rds:
        all_unique.update(rounds_data[r].keys())
    n_unique = len(all_unique)
    top_clone_freq = top_df.iloc[0][f'{last_r}_freq'] * 100 if len(top_df) > 0 else 0
    last_shannon = diversity_df[diversity_df['Round'] == last_r]['Shannon_diversity'].values
    last_shannon = float(last_shannon[0]) if len(last_shannon) > 0 else 0
    overall_fold = top_df[fold_col].median() if fold_col in top_df.columns and len(top_df) > 0 else 0
    n_families = annotations['clone_family'].nunique() if not annotations.empty else 0
    mean_prod = productivity_df['pct_productive'].mean() if productivity_df is not None and not productivity_df.empty else 0

    metrics_html = '<div class="metrics">'
    for label, val, cls in [
        ('Total Rounds', str(total_rounds), ''),
        ('Total Valid Reads', f'{total_valid:,}', 'teal'),
        ('Unique Sequences', f'{n_unique:,}', 'orange'),
        ('Top Clone Freq', f'{top_clone_freq:.2f}%', ''),
        ('Shannon Diversity', f'{last_shannon:.2f}', 'teal'),
        ('Median Enrichment', f'{overall_fold:.1f}x' if not np.isinf(overall_fold) else 'N/A', 'orange'),
        ('CDR3 Families', str(n_families), ''),
        ('Mean Productivity', f'{mean_prod:.1f}%', 'teal'),
    ]:
        metrics_html += (f'<div class="metric-card"><div class="metric-label">{label}</div>'
                         f'<div class="metric-value {cls}">{val}</div></div>')
    metrics_html += '</div>'

    # Read counts bar chart
    fig_rc = go.Figure()
    valid_reads = [all_stats[r]['valid_proteins'] for r in rds]
    fig_rc.add_trace(go.Bar(x=labels, y=valid_reads, marker_color='#4ECDC4',
                            text=[f'{v:,.0f}' for v in valid_reads], textposition='outside',
                            hovertemplate='%{x}: %{y:,.0f} reads<extra></extra>'))
    for thresh, color, name in [(2_000_000,'#e83e8c','2M ideal'),
                                 (500_000,'#fd7e14','500K'),
                                 (100_000,'#ffc107','100K min')]:
        if max(valid_reads) > thresh * 0.3:
            fig_rc.add_hline(y=thresh, line_dash='dash', line_color=color,
                             annotation_text=name, annotation_position='top left',
                             annotation_font_color=color)
    fig_rc.update_layout(**_base_layout(xaxis_title='Round', yaxis_title='Valid Reads'))
    chart_read_counts = _to_div(fig_rc, 'chart_read_counts')

    # Productivity funnel
    fig_pf = go.Figure()
    if all_stats:
        stage_colors = {'Total':'#264653','Merged':'#2A9D8F','Qual Pass':'#E9C46A',
                        'No Stop':'#F4A261','Valid':'#E76F51'}
        for stage_name, key_fn, color in [
            ('Total', lambda s: s['total_reads'], '#264653'),
            ('Merged', lambda s: s['merged'], '#2A9D8F'),
            ('Qual Pass', lambda s: s['merged'] - s.get('qual_fail', 0), '#E9C46A'),
            ('No Stop', lambda s: s['valid_proteins'] + s['length_fail'], '#F4A261'),
            ('Valid', lambda s: s['valid_proteins'], '#E76F51'),
        ]:
            vals = [key_fn(all_stats[r]) for r in rds]
            fig_pf.add_trace(go.Bar(name=stage_name, x=labels, y=vals, marker_color=color,
                                    hovertemplate=f'{stage_name}: ' + '%{y:,.0f}<extra></extra>'))
    fig_pf.update_layout(**_base_layout(barmode='group',
                          xaxis_title='Round', yaxis_title='Read Count'))
    chart_productivity = _to_div(fig_pf, 'chart_productivity')

    # Per-round radial phylogenetic trees + data tables
    treemap_html = ''
    seq_to_ann = {}
    if not annotations.empty:
        for _, a in annotations.iterrows():
            seq_to_ann[a['Sequence']] = a

    # Bubble-size legend HTML (matches _build_round_phylo sizing: (log10-2.5)*12)
    _legend_items = ''
    for lv in (3.5, 4.0, 4.5, 5.0):
        d = max(4, min(35, (lv - 2.5) * 12))
        r = d / 2
        _legend_items += (f'<span class="phylo-legend-item">'
                          f'<svg width="{d+2}" height="{d+2}">'
                          f'<circle cx="{r+1}" cy="{r+1}" r="{r}"/></svg>'
                          f'{lv}</span>')
    _legend_html = (f'<div class="phylo-legend">'
                    f'<span class="phylo-legend-label">Abundance: log10(CPM)</span>'
                    f'{_legend_items}</div>')

    for rnd in rds[1:]:
        rl = rnd.replace('Round', 'R')
        rnd_counter = rounds_data[rnd]
        total_rnd = sum(rnd_counter.values())
        if total_rnd == 0:
            continue
        top100 = rnd_counter.most_common(100)
        seqs = [s for s, _ in top100]
        cnts = {s: c for s, c in top100}

        # Build the data table rows
        table_id = f'phylo-table-{rl}'
        trows = ''
        for lidx, (seq, cnt) in enumerate(top100):
            cpm = cnt / total_rnd * 1e6
            prop = cnt / total_rnd
            ann = seq_to_ann.get(seq, {})
            cdr3 = (ann.get('CDR3', '') if isinstance(ann, dict)
                    else (ann['CDR3'] if 'CDR3' in ann.index else ''))
            family = (ann.get('clone_family', '') if isinstance(ann, dict)
                      else (ann['clone_family'] if 'clone_family' in ann.index else ''))
            trows += (f'<tr data-lidx="{lidx}">'
                      f'<td>{cdr3 if cdr3 else seq[:20]}</td>'
                      f'<td data-val="{cnt}">{cnt:,}</td>'
                      f'<td data-val="{prop:.6f}">{prop:.4f}</td>'
                      f'<td data-val="{cpm:.2f}">{cpm:,.2f}</td>'
                      f'<td>{family}</td></tr>')

        table_html = f'''<div class="phylo-paginated">
          <div class="table-controls">
            <input type="text" class="phylo-search" placeholder="Search...">
            <button class="btn-export" data-table="{table_id}">CSV</button>
            <button class="btn-export-xl" data-table="{table_id}">Excel</button>
            <span class="phylo-entries-bar">Show
              <select class="phylo-page-size"><option value="5" selected>5</option>
                <option value="10">10</option><option value="25">25</option>
                <option value="100">100</option></select> entries</span>
          </div>
          <div style="overflow-x:auto">
          <table class="enrich-table" id="{table_id}">
            <thead><tr>
              <th>CDR3<span class="sort-arrow">&#x25B2;</span></th>
              <th>Count<span class="sort-arrow">&#x25B2;</span></th>
              <th>Proportion<span class="sort-arrow">&#x25B2;</span></th>
              <th>CPM<span class="sort-arrow">&#x25B2;</span></th>
              <th>Family<span class="sort-arrow">&#x25B2;</span></th>
            </tr></thead>
            <tbody>{trows}</tbody>
          </table></div>
          <div class="phylo-page-controls">
            <div class="phylo-page-info"></div>
            <div class="phylo-page-nav"></div>
          </div>
        </div>'''

        try:
            fig_tree, rnd_lineage_json = _build_round_phylo(seqs, cnts, total_rnd, seq_to_ann, rl)
            _phylo_fig_legend = _legend("Radial phylogenetic tree of the top 100 most abundant sequences in this round, built from pairwise BLOSUM62 alignment distances using neighbor-joining. Bubble size is proportional to log\u2081\u2080(CPM). Tight clusters indicate sequence families with shared CDR3 regions. Dominant clones appear as large bubbles. Comparing trees across rounds reveals which families expand during panning.")
            if fig_tree is not None:
                rnd_lineage_script = f'<script>window.__phyloLineage_{rl}={rnd_lineage_json};</script>'
                treemap_html += (f'<div class="chart-card">'
                                 f'<div class="chart-title">Most abundant sequences in {rnd}</div>'
                                 f'{_phylo_fig_legend}'
                                 f'{_legend_html}'
                                 f'{_to_div(fig_tree, f"chart_phylo_{rl}")}'
                                 f'{rnd_lineage_script}'
                                 f'{table_html}</div>')
            else:
                treemap_html += (f'<div class="chart-card">'
                                 f'<div class="chart-title">Most abundant sequences in {rnd}</div>'
                                 f'{_phylo_fig_legend}'
                                 f'{_placeholder("🌳", "Too few sequences for tree", "Need ≥ 3 sequences")}'
                                 f'{table_html}</div>')
        except Exception as e:
            log.debug(f"Phylo tree for {rl} failed: {e}")
            _phylo_fig_legend_fallback = _legend("Radial phylogenetic tree of the top 100 most abundant sequences in this round, built from pairwise BLOSUM62 alignment distances using neighbor-joining. Bubble size is proportional to log\u2081\u2080(CPM). Tight clusters indicate sequence families with shared CDR3 regions. Dominant clones appear as large bubbles. Comparing trees across rounds reveals which families expand during panning.")
            treemap_html += (f'<div class="chart-card">'
                             f'<div class="chart-title">Most abundant sequences in {rnd}</div>'
                             f'{_phylo_fig_legend_fallback}'
                             f'{_placeholder("🌳", "Biopython required for phylogenetic trees", "pip install biopython")}'
                             f'{table_html}</div>')

    tab1_html = f'''
    <div id="tab-overview" class="tab-content active">
      {metrics_html}
      <div class="chart-row">
        <div class="chart-card"><div class="chart-title">Read Counts Per Round</div>{_legend("Bars show the number of valid protein-coding reads per panning round after quality filtering, merging, and translation. Reference lines mark recommended thresholds: green (2M ideal), yellow (500K acceptable), red (100K minimum). Low counts reduce statistical power for enrichment analysis and may indicate library preparation or sequencing issues. Counts should be roughly comparable across rounds.")}{chart_read_counts}</div>
        <div class="chart-card"><div class="chart-title">Productivity Funnel</div>{_legend("Stacked bars show how reads are filtered at each processing step: total raw reads, successfully merged pairs, quality-filtered, no stop codons, and valid protein length. A large drop at the merge step may indicate poor read overlap or adapter contamination. High stop-codon rates can suggest frameshifts. The ratio of valid proteins to total reads is the overall productivity rate.")}{chart_productivity}</div>
      </div>
      {treemap_html}
    </div>'''

    # ── Tab 2: Enrichment ─────────────────────────────────────────────────

    # Enrichment trajectories
    fig_et = go.Figure()
    for idx, row in top_df.head(10).iterrows():
        freqs = [row.get(f'{r}_freq', 0) * 100 for r in rds]
        fig_et.add_trace(go.Scatter(
            x=labels, y=freqs, mode='lines+markers', name=f'Rank {idx+1}',
            line=dict(shape='spline', width=2, color=DASHBOARD_PALETTE[idx % 10]),
            hovertemplate=f'Rank {idx+1}<br>%{{x}}: %{{y:.4f}}%<extra></extra>'))
    fig_et.update_layout(**_base_layout(xaxis_title='Round', yaxis_title='Frequency (%)'))
    chart_trajectories = _to_div(fig_et, 'chart_trajectories')

    # Log2FC vs Abundance scatter
    chart_log2fc = ''
    if len(top_df) > 0 and last_cpm_col in top_df.columns and log2fc_col in top_df.columns:
        plot_df = top_df[top_df[last_cpm_col] > 0].copy()
        if not plot_df.empty:
            fig_lf = go.Figure()
            fig_lf.add_trace(go.Scatter(
                x=plot_df[last_cpm_col], y=plot_df[log2fc_col], mode='markers',
                marker=dict(size=8, color='#45B7D1', opacity=0.7, line=dict(width=0.5, color='black')),
                text=[f'Rank {i+1}' for i in range(len(plot_df))],
                hovertemplate='%{text}<br>CPM: %{x:,.1f}<br>Log2FC: %{y:.2f}<extra></extra>'))
            fig_lf.add_hline(y=1.5, line_dash='dash', line_color='#e83e8c',
                             annotation_text='log2FC=1.5', annotation_font_color='#e83e8c')
            # Label top 5
            for i in range(min(5, len(plot_df))):
                r = plot_df.iloc[i]
                fig_lf.add_annotation(x=np.log10(r[last_cpm_col]), y=r[log2fc_col],
                                      text=f'Rank {i+1}', showarrow=False,
                                      xref='x', yref='y', font=dict(size=9, color='#e83e8c'),
                                      xshift=12, yshift=8)
            fig_lf.update_layout(**_base_layout(xaxis_title=f'{last_r} CPM',
                                  yaxis_title='Log2 Fold-Change', xaxis_type='log'))
            chart_log2fc = _to_div(fig_lf, 'chart_log2fc')
    if not chart_log2fc:
        chart_log2fc = _placeholder('&#x1F4CA;', 'Log2FC chart not available')

    # Top sequence detail (dual axis)
    chart_top_seq = ''
    if len(top_df) > 0:
        tr = top_df.iloc[0]
        fig_ts = make_subplots(specs=[[{"secondary_y": True}]])
        counts = [tr.get(f'{r}_count', 0) for r in rds]
        freqs = [tr.get(f'{r}_freq', 0) * 100 for r in rds]
        fig_ts.add_trace(go.Bar(x=labels, y=counts, name='Count', marker_color='#45B7D1',
                                hovertemplate='Count: %{y:,.0f}<extra></extra>'), secondary_y=False)
        fig_ts.add_trace(go.Scatter(x=labels, y=freqs, name='Freq %', mode='lines+markers',
                                    line=dict(color='#e83e8c', width=2),
                                    hovertemplate='Freq: %{y:.4f}%<extra></extra>'), secondary_y=True)
        fig_ts.update_layout(**_base_layout())
        fig_ts.update_yaxes(title_text='Read Count', secondary_y=False,
                            gridcolor='rgba(255,255,255,0.08)', color='#e0e0e0')
        fig_ts.update_yaxes(title_text='Frequency (%)', secondary_y=True,
                            gridcolor='rgba(255,255,255,0.08)', color='#e0e0e0')
        chart_top_seq = _to_div(fig_ts, 'chart_top_seq')

    # Enrichment table
    table_html = _build_enrichment_table(top_df, annotations, rds)

    # FASTA download buttons
    import base64
    fasta_dl_html = ''
    if fasta_enrichment or fasta_diversity:
        fasta_dl_html = '<div class="download-bar">'
        if fasta_enrichment:
            b64_e = base64.b64encode(fasta_enrichment.encode()).decode()
            fasta_dl_html += (f'<a class="download-btn" href="data:text/plain;base64,{b64_e}" '
                              f'download="top_enriched.fasta">'
                              f'<svg viewBox="0 0 24 24"><path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/></svg>'
                              f'Download Top Enriched (FASTA)</a>')
        if fasta_diversity:
            b64_d = base64.b64encode(fasta_diversity.encode()).decode()
            fasta_dl_html += (f'<a class="download-btn" href="data:text/plain;base64,{b64_d}" '
                              f'download="top_diverse.fasta">'
                              f'<svg viewBox="0 0 24 24"><path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/></svg>'
                              f'Download Top Diverse (FASTA)</a>')
        fasta_dl_html += '</div>'

    tab2_html = f'''
    <div id="tab-enrichment" class="tab-content">
      {fasta_dl_html}
      <div class="chart-row">
        <div class="chart-card"><div class="chart-title">Enrichment Trajectories</div>{_legend("Line plot tracking the frequency (CPM) of the top 10 most enriched sequences across panning rounds. Sequences that rise steeply are strong binder candidates. Sequences that plateau or decline in later rounds may indicate non-specific enrichment or competition. The trajectory shape helps distinguish true enrichment from stochastic fluctuation.")}{chart_trajectories}</div>
        <div class="chart-card"><div class="chart-title">Log2FC vs Abundance</div>{_legend("Each point is a sequence plotted by its final-round abundance (CPM, x-axis) against its log\u2082 fold-change between the last and first rounds (y-axis). Top-right hits are both abundant and highly enriched \u2014 the most promising candidates. Points with high fold-change but low CPM may be enriched by chance. The top 5 sequences are labeled.")}{chart_log2fc}</div>
      </div>
      <div class="chart-card"><div class="chart-title">Top Sequence Detail</div>{_legend("Dual-axis chart for the single most enriched sequence. Bars show the raw read count per round (left axis) and the line shows frequency percentage (right axis). This reveals whether enrichment is driven by true amplification of this clone or merely by depletion of library diversity. Ideally both count and frequency increase together.")}{chart_top_seq}</div>
      <div class="chart-card"><div class="chart-title">Enrichment Table</div>{_legend("Sortable table of the top N enriched sequences ranked by fold-change. Columns include per-round counts, frequencies, CPM values, overall fold enrichment, log\u2082FC, CDR3 sequence, and clone family assignment. Use the search box to filter by CDR3 or family. Click column headers to sort.")}{table_html}</div>
    </div>'''

    # ── Tab 3: Diversity ──────────────────────────────────────────────────

    # Shannon + Dominance
    chart_div_dom = ''
    if not diversity_df.empty:
        fig_dd = make_subplots(rows=1, cols=2, subplot_titles=('Shannon Diversity', 'Dominance'))
        dl = diversity_df['Round'].str.replace('Round', 'R')
        fig_dd.add_trace(go.Bar(x=dl.tolist(), y=diversity_df['Shannon_diversity'].tolist(),
                                marker_color='#4ECDC4', name='Shannon',
                                hovertemplate='%{x}: %{y:.2f} bits<extra></extra>'), row=1, col=1)
        fig_dd.add_trace(go.Bar(x=dl.tolist(), y=diversity_df['Top1_dominance_%'].tolist(),
                                marker_color='#e83e8c', name='Top 1%',
                                hovertemplate='Top1: %{y:.1f}%<extra></extra>'), row=1, col=2)
        fig_dd.add_trace(go.Bar(x=dl.tolist(), y=diversity_df['Top10_dominance_%'].tolist(),
                                marker_color='#FF784F', name='Top 10%', opacity=0.6,
                                hovertemplate='Top10: %{y:.1f}%<extra></extra>'), row=1, col=2)
        fig_dd.update_layout(**_base_layout(barmode='overlay'))
        chart_div_dom = _to_div(fig_dd, 'chart_div_dom')

    # Convergence heatmap
    chart_conv = ''
    if convergence_df is not None and not convergence_df.empty:
        n_rnd = len(round_names)
        heatmap_mat = np.zeros((n_rnd, n_rnd))
        rnd_to_idx = {r: i for i, r in enumerate(round_names)}
        for _, row in convergence_df.iterrows():
            i, j = rnd_to_idx[row['Round_A']], rnd_to_idx[row['Round_B']]
            heatmap_mat[i, j] = row['Jaccard']
            heatmap_mat[j, i] = row['Jaccard']
        np.fill_diagonal(heatmap_mat, 1.0)
        rl = [r.replace('Round', 'R') for r in round_names]
        text_mat = [[f'{heatmap_mat[i,j]:.3f}' for j in range(n_rnd)] for i in range(n_rnd)]
        fig_cv = go.Figure(go.Heatmap(z=heatmap_mat.tolist(), x=rl, y=rl,
                                       text=text_mat, texttemplate='%{text}',
                                       colorscale='YlOrRd', zmin=0, zmax=1,
                                       hovertemplate='%{x} vs %{y}: %{z:.3f}<extra></extra>'))
        fig_cv.update_layout(**_base_layout())
        chart_conv = _to_div(fig_cv, 'chart_convergence')
    if not chart_conv:
        chart_conv = _placeholder('&#x1F50D;', 'Convergence data not available')

    # Saturation curves
    chart_sat = ''
    if saturation_df is not None and not saturation_df.empty:
        fig_sat = go.Figure()
        for idx_r, rnd in enumerate(sorted(saturation_df['Round'].unique(),
                                           key=lambda x: int(x.replace('Round', '')))):
            rnd_data = saturation_df[saturation_df['Round'] == rnd]
            fig_sat.add_trace(go.Scatter(
                x=rnd_data['sample_size'], y=rnd_data['number_unique'],
                mode='lines', name=rnd.replace('Round', 'R'),
                line=dict(color=DASHBOARD_PALETTE[idx_r % 10], width=2),
                hovertemplate=f'{rnd.replace("Round","R")}<br>Samples: %{{x:,.0f}}<br>Unique: %{{y:,.0f}}<extra></extra>'))
        if diversity_est is not None and not diversity_est.empty:
            for _, de_row in diversity_est.iterrows():
                fig_sat.add_hline(y=de_row['approx_diversity'], line_dash='dot',
                                  line_color='#a0a0b0', opacity=0.5,
                                  annotation_text=f'{de_row["Round"].replace("Round","R")} est',
                                  annotation_font_size=9, annotation_font_color='#a0a0b0')
        fig_sat.update_layout(**_base_layout(xaxis_title='Sample Size',
                              yaxis_title='Unique Sequences'))
        chart_sat = _to_div(fig_sat, 'chart_saturation')
    if not chart_sat:
        chart_sat = _placeholder('&#x1F4C8;', 'Saturation data not available')

    # Clone family dynamics
    chart_clone_dyn = ''
    if clone_tracking is not None and not clone_tracking.empty:
        freq_cols = [f'{r}_freq' for r in round_names if f'{r}_freq' in clone_tracking.columns]
        if freq_cols:
            top_fams = clone_tracking.head(8)
            fig_cd = go.Figure()
            x_labels = [r.replace('Round', 'R') for r in round_names if f'{r}_freq' in clone_tracking.columns]
            for idx_f, (_, frow) in enumerate(top_fams.iterrows()):
                vals = [float(frow.get(c, 0)) * 100 for c in freq_cols]
                fig_cd.add_trace(go.Scatter(
                    x=x_labels, y=vals, name=frow['clone_family'],
                    mode='lines', stackgroup='one',
                    line=dict(color=DASHBOARD_PALETTE[idx_f % 10]),
                    hovertemplate=f'{frow["clone_family"]}<br>%{{x}}: %{{y:.3f}}%<extra></extra>'))
            fig_cd.update_layout(**_base_layout(xaxis_title='Round',
                                  yaxis_title='Frequency (%)'))
            chart_clone_dyn = _to_div(fig_cd, 'chart_clone_dyn')
    if not chart_clone_dyn:
        chart_clone_dyn = _placeholder('&#x1F9EC;', 'Clone tracking data not available')

    tab3_html = f'''
    <div id="tab-diversity" class="tab-content">
      <div class="chart-row">
        <div class="chart-card"><div class="chart-title">Shannon Diversity &amp; Dominance</div>{_legend("<b>Left:</b> Shannon diversity index per round \u2014 higher values indicate more even sequence distributions. A steep drop between rounds signals strong selection. <b>Right:</b> Top 1% and Top 10% dominance \u2014 the fraction of total reads captured by the most abundant 1% or 10% of unique sequences. Rising dominance confirms convergence toward a few dominant clones during panning.")}{chart_div_dom}</div>
        <div class="chart-card"><div class="chart-title">Convergence Heatmap</div>{_legend("Jaccard similarity between each pair of rounds, measuring the overlap in unique sequence sets. Values range from 0 (no shared sequences) to 1 (identical repertoires). High similarity between consecutive late rounds indicates the library has converged. Low similarity between early and late rounds confirms successful selection pressure.")}{chart_conv}</div>
      </div>
      <div class="chart-row">
        <div class="chart-card"><div class="chart-title">Saturation Curves</div>{_legend("Rarefaction curves showing the number of unique sequences discovered as a function of sampling depth for each round. A curve that plateaus has been sequenced to sufficient depth \u2014 most unique sequences have been observed. Curves still rising linearly indicate under-sampling. Horizontal dashed lines show estimated total diversity using the Chao1 estimator.")}{chart_sat}</div>
        <div class="chart-card"><div class="chart-title">Clone Family Dynamics</div>{_legend("Stacked area chart tracking the top 8 clone families by frequency across rounds. Families that expand to dominate later rounds are strong enrichment candidates. The &quot;Other&quot; category represents all remaining families. A single family dominating &gt;50% of late rounds suggests very strong selection. Competing families may indicate multiple binding epitopes or modes.")}{chart_clone_dyn}</div>
      </div>
    </div>'''

    # ── Tab 4: Clusters ───────────────────────────────────────────────────

    # CDR3 length distribution
    chart_cdr3_len = ''
    if not annotations.empty:
        cdr3_lens = annotations[annotations['CDR3_length'] > 0]['CDR3_length']
        if len(cdr3_lens) > 0:
            med = float(cdr3_lens.median())
            fig_cl = go.Figure()
            fig_cl.add_trace(go.Histogram(
                x=cdr3_lens.tolist(), marker_color='#45B7D1', opacity=0.8,
                hovertemplate='Length: %{x}<br>Count: %{y}<extra></extra>'))
            fig_cl.add_vline(x=med, line_dash='dash', line_color='#e83e8c',
                             annotation_text=f'Median: {med:.0f}', annotation_font_color='#e83e8c')
            fig_cl.update_layout(**_base_layout(xaxis_title='CDR3 Length (aa)',
                                  yaxis_title='Count'))
            chart_cdr3_len = _to_div(fig_cl, 'chart_cdr3_len')
    if not chart_cdr3_len:
        chart_cdr3_len = _placeholder('&#x1F4CF;', 'CDR3 length data not available')

    # CDR3 length vs enrichment bubble
    chart_cdr3_enrich = ''
    if not annotations.empty and len(top_df) > 0 and fold_col in top_df.columns:
        freq_col = f'{last_r}_freq'
        if freq_col in top_df.columns:
            merged_cdr3 = top_df.merge(annotations[['Sequence', 'CDR3_length']], on='Sequence', how='inner')
            merged_cdr3 = merged_cdr3[merged_cdr3['CDR3_length'] > 0]
            if not merged_cdr3.empty:
                max_freq = merged_cdr3[freq_col].max()
                sizes = (merged_cdr3[freq_col] / max_freq * 40 + 5).tolist() if max_freq > 0 else [10] * len(merged_cdr3)
                fig_ce = go.Figure()
                fig_ce.add_trace(go.Scatter(
                    x=merged_cdr3['CDR3_length'].tolist(), y=merged_cdr3[fold_col].tolist(),
                    mode='markers', marker=dict(size=sizes, color='#4ECDC4', opacity=0.6,
                                                line=dict(width=0.5, color='black')),
                    text=[f'Rank {i+1}' for i in range(len(merged_cdr3))],
                    hovertemplate='%{text}<br>CDR3 len: %{x}<br>Fold: %{y:.1f}x<extra></extra>'))
                fig_ce.update_layout(**_base_layout(xaxis_title='CDR3 Length (aa)',
                                      yaxis_title='Enrichment Fold', yaxis_type='log'))
                chart_cdr3_enrich = _to_div(fig_ce, 'chart_cdr3_enrich')
    if not chart_cdr3_enrich:
        chart_cdr3_enrich = _placeholder('&#x1F4CA;', 'CDR3 enrichment data not available')

    # Cluster PCA — one trace per family for filter toggles, alpseq styling
    chart_pca = ''
    pca_table_html = ''
    pca_filter_html = ''
    if cluster_pca is not None and not cluster_pca.empty:
        max_cpm = cluster_pca['cpm'].max() if cluster_pca['cpm'].max() > 0 else 1
        fig_pca = go.Figure()
        families = cluster_pca['clone_family'].unique()

        # Classify sequences: enrichment top 100 vs abundance top 100
        enrich_top100 = set()
        if fold_col in top_df.columns and len(top_df) > 0:
            enrich_sorted = top_df.sort_values(fold_col, ascending=False).head(100)
            enrich_top100 = set(enrich_sorted['Sequence'].tolist())
        abund_top100 = set()
        last_rnd_counter = rounds_data.get(last_r, Counter())
        for sq, _ in last_rnd_counter.most_common(100):
            abund_top100.add(sq)

        for idx_f, fam in enumerate(families):
            sub = cluster_pca[cluster_pca['clone_family'] == fam]
            sizes = (np.log10(sub['cpm'].clip(lower=1)) / np.log10(max(max_cpm, 10)) * 25 + 6).tolist()
            # Classify each sequence
            cats = []
            for sq in sub['Sequence'].tolist() if 'Sequence' in sub.columns else [''] * len(sub):
                in_e = sq in enrich_top100
                in_a = sq in abund_top100
                if in_e and in_a:
                    cats.append('both')
                elif in_a:
                    cats.append('diversity')
                elif in_e:
                    cats.append('enrichment')
                else:
                    cats.append('neither')
            fig_pca.add_trace(go.Scatter(
                x=sub['PC1'].tolist(), y=sub['PC2'].tolist(),
                mode='markers', name=str(fam),
                customdata=list(zip(
                    sub['rank'].tolist() if 'rank' in sub.columns else list(range(len(sub))),
                    cats)),
                marker=dict(size=sizes, color='#DC7F9B', opacity=0.7,
                            line=dict(width=1, color='#760a2a')),
                text=[f'{fam}<br>CPM: {c:,.0f}' for c in sub['cpm']],
                hovertemplate='%{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>'))
        fig_pca.update_layout(**_base_layout(xaxis_title='PC1', yaxis_title='PC2',
                              dragmode='lasso', showlegend=len(families) <= 15))
        chart_pca = _to_div(fig_pca, 'chart_pca')

        # PCA filter checkboxes — alpseq style
        pca_filter_html = '''<div class="filter-bar" id="pca-filters">
          <span style="font-weight:600;color:var(--text);font-size:13px;margin-right:8px">Show top 100 sequences</span>
          <label><input type="checkbox" value="both" checked>In both top 100s</label>
          <label><input type="checkbox" value="diversity" checked>In diversity top 100 only</label>
          <label><input type="checkbox" value="enrichment" checked>In enrichment top 100 only</label>
          <label><input type="checkbox" value="neither" checked>In neither top 100</label>
        </div>'''

        # PCA-linked table
        pca_rows = []
        for i, (_, pr) in enumerate(cluster_pca.iterrows()):
            rank = int(pr.get('rank', i + 1)) if 'rank' in pr.index else i + 1
            cdr3 = str(pr.get('CDR3', ''))[:12] if 'CDR3' in pr.index else ''
            fam = str(pr.get('clone_family', ''))
            cpm = float(pr.get('cpm', 0))
            enr = float(pr.get('enrichment', 0))
            is_top = 'Yes' if rank <= 100 else ''
            pca_rows.append(
                f'<tr data-rank="{rank}" data-family="{fam}">'
                f'<td data-val="{rank}">{rank}</td>'
                f'<td>{cdr3}</td><td>{fam}</td>'
                f'<td data-val="{cpm:.1f}">{cpm:,.0f}</td>'
                f'<td data-val="{enr:.2f}">{enr:.2f}</td>'
                f'<td>{is_top}</td></tr>')
        pca_th = ''.join(f'<th>{c}<span class="sort-arrow">&#x25B2;</span></th>'
                         for c in ['Rank','CDR3','Family','CPM','Log2FC','Top 100'])
        pca_table_html = f'''<div class="phylo-paginated">
        <div class="table-controls">
          <input type="text" class="phylo-search" placeholder="Search PCA table...">
          <button class="btn-export" data-table="pca-table">CSV</button>
          <button class="btn-export-xl" data-table="pca-table">Excel</button>
          <span class="phylo-entries-bar">Show
            <select class="phylo-page-size"><option value="10" selected>10</option>
              <option value="25">25</option><option value="50">50</option>
              <option value="100">100</option></select> entries</span>
        </div>
        <div style="overflow-x:auto">
        <table class="enrich-table" id="pca-table">
          <thead><tr>{pca_th}</tr></thead>
          <tbody>{''.join(pca_rows)}</tbody>
        </table></div>
        <div class="phylo-page-controls">
          <div class="phylo-page-info"></div>
          <div class="phylo-page-nav"></div>
        </div>
        </div>'''
    if not chart_pca:
        chart_pca = _placeholder('&#x1F52C;', 'PCA not available', 'Requires: biopython, scikit-learn')

    # Variant heatmap
    chart_variant = ''
    if not annotations.empty and len(top_df) > 0:
        top_seq = top_df.iloc[0]['Sequence']
        top_ann = annotations[annotations['Sequence'] == top_seq]
        if not top_ann.empty:
            dom_family = top_ann.iloc[0]['clone_family']
            family_seqs = annotations[annotations['clone_family'] == dom_family]['Sequence'].tolist()
            variants = top_df[top_df['Sequence'].isin(family_seqs)].head(20)
            if len(variants) >= 2:
                seqs = variants['Sequence'].tolist()
                ref = seqs[0]
                max_len = max(len(s) for s in seqs)
                mat = np.zeros((len(seqs), max_len), dtype=int)
                for i, seq in enumerate(seqs):
                    for j in range(max_len):
                        if j < len(seq) and j < len(ref):
                            mat[i, j] = 0 if seq[j] == ref[j] else 1
                        else:
                            mat[i, j] = 2
                var_cols = [j for j in range(max_len) if np.any(mat[:, j] != 0)]
                if var_cols:
                    sub_mat = mat[:, var_cols]
                    cell_text = []
                    for i, seq in enumerate(seqs):
                        row_t = [seq[j] if j < len(seq) else '-' for j in var_cols]
                        cell_text.append(row_t)
                    ylabels = [f'Rank {int(variants.iloc[i].name) + 1}' for i in range(len(seqs))]
                    xlabels = [str(c + 1) for c in var_cols]
                    fig_vh = go.Figure(go.Heatmap(
                        z=sub_mat.tolist(), x=xlabels, y=ylabels,
                        text=cell_text, texttemplate='%{text}',
                        colorscale=[[0,'#E8F5E9'],[0.5,'#FF8A80'],[1,'#BDBDBD']],
                        showscale=False,
                        hovertemplate='Pos %{x}: %{text}<extra></extra>'))
                    fig_vh.update_layout(**_base_layout(
                        xaxis_title='Position', yaxis=dict(autorange='reversed')))
                    chart_variant = _to_div(fig_vh, 'chart_variant')
    if not chart_variant:
        chart_variant = _placeholder('&#x1F9EC;', 'Variant heatmap not available',
                                     'Requires >=2 variants in dominant family')

    # Enriched clusters: radial phylo tree + data table from clone_tracking
    cluster_tree_html = ''
    cluster_table_html = ''
    if clone_tracking is not None and not clone_tracking.empty:
        family_sizes = annotations.groupby('clone_family').size().to_dict() if not annotations.empty else {}

        # ── Radial tree of top 100 enriched clusters ──
        cl_fam_to_idx = {}
        try:
            cl_fig, cl_matrix, cl_lineage_json, cl_fam_to_idx_ret, cl_l2fc_range = _build_cluster_phylo(
                clone_tracking, annotations, rounds_data, round_names, max_clusters=100)
            if cl_fam_to_idx_ret:
                cl_fam_to_idx = cl_fam_to_idx_ret
            if cl_fig is not None:
                cl_subtitle = f'using {cl_matrix}'
                lineage_script = f'<script>window.__clusterLineage={cl_lineage_json};</script>'
                # Log2FC color bar legend
                l2fc_lo, l2fc_hi = cl_l2fc_range if cl_l2fc_range else (0, 0)
                l2fc_bar_html = (
                    f'<div class="phylo-legend" style="margin-top:6px">'
                    f'<span class="phylo-legend-label">Log2FC:</span>'
                    f'<span style="font-size:11px;color:var(--text-dim);font-family:var(--mono)">{l2fc_lo:.1f}</span>'
                    f'<span style="display:inline-block;width:120px;height:12px;border-radius:4px;'
                    f'background:linear-gradient(to right,#2d0a1e,#8b2252,#e83e8c);'
                    f'vertical-align:middle;margin:0 4px"></span>'
                    f'<span style="font-size:11px;color:var(--text-dim);font-family:var(--mono)">{l2fc_hi:.1f}</span>'
                    f'</div>')
                cluster_tree_html = (
                    f'<div style="color:var(--text-dim);font-size:12px;font-family:var(--mono);'
                    f'margin-bottom:4px">{cl_subtitle}</div>'
                    f'{_legend_html}'
                    f'{_to_div(cl_fig, "chart_cluster_phylo")}'
                    f'{l2fc_bar_html}'
                    f'{lineage_script}')
            else:
                cluster_tree_html = _placeholder(
                    '&#x1F333;', 'Too few clusters for tree', 'Need >= 3 enriched clusters')
        except Exception as e:
            log.debug(f"Cluster phylo tree failed: {e}")
            cluster_tree_html = _placeholder(
                '&#x1F333;', 'Biopython required for phylogenetic trees',
                'pip install biopython')

        # ── Data table ──
        ct_rows = ''
        for _, crow in clone_tracking.iterrows():
            fam = crow['clone_family']
            cdr3_raw = crow.get('representative_CDR3', '')
            cdr3 = str(cdr3_raw) if pd.notna(cdr3_raw) and cdr3_raw != '' else ''
            # Fall back to truncated sequence when CDR3 annotation is empty
            if not cdr3:
                if not annotations.empty:
                    fam_rows = annotations[annotations['clone_family'] == fam]
                    if not fam_rows.empty:
                        seq = str(fam_rows.iloc[0]['Sequence'])
                        cdr3 = (seq[:20] + '\u2026') if len(seq) > 20 else seq
                if not cdr3:
                    cdr3 = fam
            size = family_sizes.get(fam, 1)
            cidx = cl_fam_to_idx.get(fam)
            cidx_attr = f' data-cidx="{cidx}"' if cidx is not None else ''
            cells = (f'<td title="{cdr3}">{cdr3}</td>'
                     f'<td data-val="{size}">{size}</td>')
            for rnd in round_names:
                cpm = crow.get(f'{rnd}_cpm', 0)
                cells += f'<td data-val="{cpm:.2f}">{cpm:,.2f}</td>'
            enr = crow.get('enrichment', 0)
            fold_str = f'{enr:.2f}' if not np.isinf(enr) else 'de novo'
            cells += f'<td data-val="{enr if not np.isinf(enr) else 1e9}">{fold_str}</td>'
            l2fc = crow.get('log2fc', 0)
            cells += f'<td data-val="{l2fc:.2f}">{l2fc:.2f}</td>'
            ct_rows += f'<tr{cidx_attr}>{cells}</tr>'

        ct_headers = (['cluster_lead', 'cluster_size']
                      + [f'{r.replace("Round", "R")}_CPM' for r in round_names]
                      + ['Fold', 'Log2FC'])
        ct_th = ''.join(f'<th>{h}<span class="sort-arrow">&#x25B2;</span></th>' for h in ct_headers)
        cluster_table_html = f'''<div class="phylo-paginated">
          <div class="table-controls">
            <input type="text" class="phylo-search" placeholder="Search clusters...">
            <button class="btn-export" data-table="cluster-enrich-table">CSV</button>
            <button class="btn-export-xl" data-table="cluster-enrich-table">Excel</button>
            <span class="phylo-entries-bar">Show
              <select class="phylo-page-size"><option value="10" selected>10</option>
                <option value="25">25</option><option value="50">50</option>
                <option value="100">100</option></select> entries</span>
          </div>
          <div style="overflow-x:auto">
          <table class="enrich-table" id="cluster-enrich-table">
            <thead><tr>{ct_th}</tr></thead>
            <tbody>{ct_rows}</tbody>
          </table></div>
          <div class="phylo-page-controls">
            <div class="phylo-page-info"></div>
            <div class="phylo-page-nav"></div>
          </div>
        </div>'''
    if not cluster_tree_html and not cluster_table_html:
        cluster_tree_html = _placeholder('&#x1F4CA;', 'No cluster tracking data available')

    tab4_html = f'''
    <div id="tab-clusters" class="tab-content">
      <div class="chart-row">
        <div class="chart-card"><div class="chart-title">CDR3 Length Distribution</div>{_legend("Histogram of CDR3 loop lengths (amino acids) across all annotated sequences. The vertical dashed line marks the median. Nanobody CDR3 loops are typically 10\u201325 residues. A narrow distribution may indicate limited library diversity or selection for a preferred CDR3 length. Bimodal distributions can indicate distinct structural classes.")}{chart_cdr3_len}</div>
        <div class="chart-card"><div class="chart-title">CDR3 Length vs Enrichment</div>{_legend("Bubble scatter where each bubble is a sequence positioned by CDR3 length (x-axis) and enrichment fold (y-axis). Bubble size reflects frequency. This reveals whether enrichment favors particular CDR3 lengths \u2014 often a specific length range dominates among top binders, reflecting the structural requirements of the target epitope.")}{chart_cdr3_enrich}</div>
      </div>
      <div class="chart-card">
        <div class="chart-title">PCA plot &amp; top 100s</div>
        {_legend("Principal component analysis of the top enriched sequences based on pairwise BLOSUM62 alignment distances. Points are colored by clone family. Clusters of same-colored points indicate related sequences. Use the checkboxes to filter by whether sequences appear in the diversity top 100, enrichment top 100, both, or neither. Lasso-select points to highlight them in the linked table below.")}
        {pca_filter_html}
        {chart_pca}
        {pca_table_html}
      </div>
      <div class="chart-card">
        <div class="chart-title">Enriched Clusters</div>
        {_legend("Phylogenetic tree and table of clone family representatives (cluster leads). The tree shows structural relationships between the most enriched families. The table includes cluster size, per-round CPM, fold enrichment, and log\u2082FC. Large clusters with high enrichment represent the most promising candidate families for downstream characterization.")}
        {cluster_tree_html}
        {cluster_table_html}
      </div>
      <div class="chart-card"><div class="chart-title">Variant Heatmap</div>{_legend("Heatmap of sequence variants within the most dominant clone family. Only variable positions are shown. Green = match to consensus, red = mismatch, gray = gap. This reveals the mutational landscape within a family \u2014 positions with frequent variation are likely tolerant to substitution, while conserved positions may be critical for binding.")}{chart_variant}</div>
    </div>'''

    # ── Tab 5: Scoring (conditional) ──────────────────────────────────────
    has_ablang = ablang_df is not None and not ablang_df.empty
    has_igblast = gene_usage_df is not None and not gene_usage_df.empty
    show_scoring = has_ablang or has_igblast

    chart_ppl = ''
    ppl_table_html = ''
    if has_ablang and len(top_df) > 0 and fold_col in top_df.columns:
        merged_ab = top_df.merge(ablang_df, on='Sequence', how='inner')
        if not merged_ab.empty:
            fam_map = dict(zip(annotations['Sequence'], annotations['clone_family']))
            cdr3_map = dict(zip(annotations['Sequence'], annotations['CDR3'])) if not annotations.empty and 'CDR3' in annotations.columns else {}
            merged_ab['clone_family'] = merged_ab['Sequence'].map(fam_map).fillna('unknown')
            merged_ab['CDR3'] = merged_ab['Sequence'].map(cdr3_map).fillna('')
            merged_ab = merged_ab.reset_index(drop=True)
            fig_ppl = go.Figure()
            families = merged_ab['clone_family'].unique()
            for idx_f, fam in enumerate(families):
                sub = merged_ab[merged_ab['clone_family'] == fam]
                fig_ppl.add_trace(go.Scatter(
                    x=sub[fold_col].tolist(), y=sub['ablang_pseudo_ppl'].tolist(),
                    mode='markers', name=fam if len(families) <= 10 else None,
                    customdata=sub.index.tolist(),
                    marker=dict(size=8, color=DASHBOARD_PALETTE[idx_f % 10], opacity=0.7,
                                line=dict(width=0.5, color='black')),
                    hovertemplate=f'{fam}<br>Fold: %{{x:.1f}}<br>PPL: %{{y:.3f}}<extra></extra>'))
            fig_ppl.update_layout(**_base_layout(xaxis_title='Enrichment Fold (log)',
                                  yaxis_title='AbLang Pseudo-PPL',
                                  xaxis_type='log', showlegend=len(families) <= 10))
            chart_ppl = _to_div(fig_ppl, 'chart_ppl')

            # Build PPL table
            ppl_rows = ''
            for pidx, prow in merged_ab.iterrows():
                cdr3 = str(prow.get('CDR3', ''))[:15]
                fam = str(prow.get('clone_family', ''))
                fold = float(prow.get(fold_col, 0))
                ppl = float(prow.get('ablang_pseudo_ppl', 0))
                mean_p = float(prow.get('ablang_mean_prob', 0))
                cpm = float(prow.get(last_cpm_col, 0)) if last_cpm_col in prow.index else 0
                fold_str = f'{fold:.2f}' if not np.isinf(fold) else 'de novo'
                fold_val = fold if not np.isinf(fold) else 1e9
                ppl_rows += (f'<tr data-pidx="{pidx}">'
                             f'<td title="{cdr3}">{cdr3}</td>'
                             f'<td>{fam}</td>'
                             f'<td data-val="{fold_val:.2f}">{fold_str}</td>'
                             f'<td data-val="{ppl:.4f}">{ppl:.3f}</td>'
                             f'<td data-val="{mean_p:.4f}">{mean_p:.3f}</td>'
                             f'<td data-val="{cpm:.2f}">{cpm:,.0f}</td></tr>')
            ppl_th = ''.join(f'<th>{h}<span class="sort-arrow">&#x25B2;</span></th>'
                             for h in ['CDR3', 'Family', 'Fold', 'PPL', 'Mean Prob', 'CPM'])
            ppl_table_html = f'''<div class="phylo-paginated">
              <div class="table-controls">
                <input type="text" class="phylo-search" placeholder="Search...">
                <button class="btn-export" data-table="ppl-table">CSV</button>
                <button class="btn-export-xl" data-table="ppl-table">Excel</button>
                <span class="phylo-entries-bar">Show
                  <select class="phylo-page-size"><option value="10" selected>10</option>
                    <option value="25">25</option><option value="50">50</option>
                    <option value="100">100</option></select> entries</span>
              </div>
              <div style="overflow-x:auto">
              <table class="enrich-table" id="ppl-table">
                <thead><tr>{ppl_th}</tr></thead>
                <tbody>{ppl_rows}</tbody>
              </table></div>
              <div class="phylo-page-controls">
                <div class="phylo-page-info"></div>
                <div class="phylo-page-nav"></div>
              </div>
            </div>'''
    if not chart_ppl:
        chart_ppl = _placeholder('&#x1F9EA;', 'AbLang data not available', 'Run with --ablang flag')

    chart_gene = ''
    if has_igblast:
        gene_types = [gt for gt in ['V', 'D', 'J'] if gt in gene_usage_df['gene_type'].values]
        if gene_types:
            fig_gu = make_subplots(rows=1, cols=len(gene_types))
            for idx_g, gt in enumerate(gene_types):
                gt_data = gene_usage_df[gene_usage_df['gene_type'] == gt].sort_values('percentage', ascending=True)
                fig_gu.add_trace(go.Bar(
                    y=gt_data['gene'].tolist(), x=gt_data['percentage'].tolist(),
                    orientation='h', marker_color=DASHBOARD_PALETTE[idx_g % 10],
                    hovertemplate='%{y}: %{x:.1f}%<extra></extra>',
                    showlegend=False), row=1, col=idx_g + 1)
            fig_gu.update_layout(**_base_layout())
            chart_gene = _to_div(fig_gu, 'chart_gene_usage')
    if not chart_gene:
        chart_gene = _placeholder('&#x1F9EC;', 'IgBLAST data not available', 'Run with --igblast flag')

    tab5_html = ''
    if show_scoring:
        tab5_html = f'''
        <div id="tab-scoring" class="tab-content">
          <div class="chart-card">
            <div class="chart-title">PPL vs Enrichment</div>
            {_legend("Scatter plot of AbLang pseudo-perplexity (y-axis) versus enrichment fold (x-axis), colored by clone family. AbLang PPL measures how &quot;antibody-like&quot; a sequence is \u2014 lower PPL means the sequence better matches natural antibody language patterns. Ideal candidates appear in the bottom-right: highly enriched and predicted to be well-folded. High PPL outliers may have unusual frameworks or mutations that warrant inspection.")}
            {chart_ppl}
            {ppl_table_html}
          </div>
          <div class="chart-card"><div class="chart-title">V/D/J Gene Usage</div>{_legend("Horizontal bar charts showing the percentage usage of germline V, D, and J gene segments among annotated sequences. Dominant V genes indicate the germline framework most compatible with binding. Comparing gene usage to the naive library can reveal germline-level selection biases introduced by panning.")}{chart_gene}</div>
        </div>'''

    # ── Tab 6: Methods (always shown) ───────────────────────────────────
    mermaid_src = '\n'.join([
        'graph TD',
        'A[Raw FASTQ] --> B[Merge and QC]',
        'B --> C[Translate]',
        'C --> D[Length Filter]',
        'D --> E[Enrichment Analysis]',
        'E --> F[CDR Annotation]',
        'F --> G[CDR3 Clustering]',
        'G --> H[Chimera Detection]',
        'H --> I[Clone Tracking]',
        'I --> J[Saturation Analysis]',
        'J --> K{Optional Modules}',
        'K -->|ablang| L[AbLang Scoring]',
        'K -->|igblast| M[IgBLAST V/D/J]',
        'K --> N[Cluster PCA]',
        'L --> O[Dashboard and Reports]',
        'M --> O',
        'N --> O',
        'J --> O',
    ])
    tab6_html = f'''
    <div id="tab-methods" class="tab-content">
      <div class="methods-text">
        <h3>Analysis Pipeline</h3>
        <p>Nabbit processes paired-end Illumina FASTQ data through a multi-step pipeline.
        Raw reads are merged by overlap assembly, quality-filtered, and translated to protein sequences.
        Enrichment is computed across panning rounds, followed by CDR annotation, CDR3-based clustering,
        chimera detection, and cross-round clone tracking. Optional modules include AbLang language-model
        scoring, IgBLAST V/D/J gene assignment, and cluster PCA visualization.</p>
      </div>
      <div class="mermaid-wrap">
        <pre class="mermaid">{mermaid_src}</pre>
      </div>
    </div>'''

    # ── Assemble HTML ─────────────────────────────────────────────────────
    scoring_tab = '<div class="tab" data-tab="tab-scoring">Scoring</div>' if show_scoring else ''

    html = f'''<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Nabbit Dashboard</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=DM+Sans:wght@400;500;600;700&display=swap" rel="stylesheet">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.min.js"></script>
<script>mermaid.initialize({{startOnLoad:false,theme:'dark'}});</script>
<style>{DASHBOARD_CSS}</style>
</head>
<body>
<div class="container">
<div class="header">
  <div class="header-left">
    <div class="logo"><svg xmlns="http://www.w3.org/2000/svg" viewBox="100 90 300 280" width="40" height="40" style="vertical-align:middle;margin-right:10px"><defs><filter id="glow" x="-20%" y="-20%" width="140%" height="140%"><feGaussianBlur stdDeviation="5" result="blur"/><feComposite in="SourceGraphic" in2="blur" operator="over"/></filter></defs><ellipse cx="250" cy="330" rx="110" ry="35" fill="#2F3542"/><path d="M180 330C180 230,320 230,320 330Z" fill="#FF9F43"/><g transform="rotate(-20 200 160)"><ellipse cx="200" cy="160" rx="20" ry="65" fill="#FF9F43" stroke="#E17055" stroke-width="2"/><ellipse cx="200" cy="165" rx="10" ry="45" fill="#FFD1D1"/></g><g transform="rotate(25 300 160)"><ellipse cx="300" cy="160" rx="20" ry="65" fill="#FF9F43" stroke="#E17055" stroke-width="2"/><ellipse cx="300" cy="165" rx="10" ry="45" fill="#FFD1D1"/></g><ellipse cx="250" cy="245" rx="70" ry="60" fill="#FF9F43"/><ellipse cx="230" cy="265" rx="25" ry="20" fill="#FFF"/><ellipse cx="270" cy="265" rx="25" ry="20" fill="#FFF"/><path d="M240 255L260 255L250 265Z" fill="#FF4757"/><ellipse cx="220" cy="235" rx="8" ry="14" fill="#2F3542"/><circle cx="218" cy="230" r="3" fill="#FFF"/><ellipse cx="280" cy="235" rx="8" ry="14" fill="#2F3542"/><circle cx="278" cy="230" r="3" fill="#FFF"/><ellipse cx="185" cy="335" rx="22" ry="16" fill="#FF9F43" stroke="#E17055" stroke-width="1.5"/><path d="M175 345L175 330M195 345L195 330" stroke="#E17055" stroke-width="2" stroke-linecap="round"/><g filter="url(#glow)"><path d="M315 285L340 235L325 185M340 235L380 200" stroke="#1E90FF" stroke-width="18" stroke-linecap="round" stroke-linejoin="round" fill="none"/><path d="M315 285L340 235L325 185M340 235L380 200" stroke="#70A1FF" stroke-width="8" stroke-linecap="round" stroke-linejoin="round" fill="none"/><path d="M315 285L340 235L325 185M340 235L380 200" stroke="#FFF" stroke-width="3" stroke-linecap="round" stroke-linejoin="round" fill="none"/><path d="M380 160L383 168L391 171L383 174L380 182L377 174L369 171L377 168Z" fill="#70A1FF"/><circle cx="310" cy="180" r="4" fill="#70A1FF"/><circle cx="395" cy="220" r="3" fill="#70A1FF"/></g><ellipse cx="310" cy="285" rx="22" ry="18" fill="#FF9F43" stroke="#E17055" stroke-width="1.5" transform="rotate(-15 310 285)"/><path d="M305 295L315 275M320 295L330 275" stroke="#E17055" stroke-width="2" stroke-linecap="round"/></svg><span>Nabbit</span></div>
    <div class="subtitle">Nanobody NGS Enrichment Dashboard &middot; {time.strftime('%Y-%m-%d %H:%M:%S')}</div>
  </div>
  <div class="header-right">
    <button class="theme-btn" id="theme-toggle">Light Mode</button>
  </div>
</div>
<div class="tabs">
  <div class="tab active" data-tab="tab-overview">Overview</div>
  <div class="tab" data-tab="tab-enrichment">Enrichment</div>
  <div class="tab" data-tab="tab-diversity">Diversity</div>
  <div class="tab" data-tab="tab-clusters">Clusters</div>
  {scoring_tab}
  <div class="tab" data-tab="tab-methods">Methods</div>
</div>
{tab1_html}
{tab2_html}
{tab3_html}
{tab4_html}
{tab5_html}
{tab6_html}
</div>
<script>{DASHBOARD_JS}</script>
</body>
</html>'''

    out_path = os.path.join(output_dir, 'dashboard.html')
    with open(out_path, 'w') as f:
        f.write(html)
    log.info(f"Interactive dashboard: {out_path}")
    return out_path


def _build_enrichment_table(top_df, annotations, rds):
    """Build a sortable/filterable HTML table for the enrichment data."""
    seq_to_ann = {}
    if not annotations.empty:
        for _, a in annotations.iterrows():
            seq_to_ann[a['Sequence']] = a

    header_cols = ['Rank', 'Len']
    for r in rds:
        rl = r.replace('Round', 'R')
        header_cols.extend([f'{rl} Cnt', f'{rl} %', f'{rl} CPM'])
    header_cols.extend(['Fold', 'Log2FC', 'CDR3', 'CDR3 Len', 'Family'])

    first_r, last_r = rds[0], rds[-1]
    fold_col = f'{first_r}_to_{last_r}_fold'
    log2fc_col = f'{first_r}_to_{last_r}_log2fc'

    rows_html = []
    for idx, row in top_df.head(50).iterrows():
        seq = row['Sequence']
        ann = seq_to_ann.get(seq, {})
        cdr3 = ann.get('CDR3', '') if isinstance(ann, dict) else (ann['CDR3'] if 'CDR3' in ann.index else '')
        cdr3_len = ann.get('CDR3_length', 0) if isinstance(ann, dict) else (ann['CDR3_length'] if 'CDR3_length' in ann.index else 0)
        family = ann.get('clone_family', '') if isinstance(ann, dict) else (ann['clone_family'] if 'clone_family' in ann.index else '')

        cdr3_display = cdr3[:12] + '...' if len(str(cdr3)) > 12 else str(cdr3)
        cells = [f'<td data-val="{idx+1}">{idx+1}</td>',
                 f'<td data-val="{int(row["Length"])}">{int(row["Length"])}</td>']
        for r in rds:
            cnt = int(row.get(f'{r}_count', 0))
            freq = row.get(f'{r}_freq', 0) * 100
            cpm = row.get(f'{r}_cpm', 0)
            cells.append(f'<td data-val="{cnt}">{cnt:,}</td>')
            cells.append(f'<td data-val="{freq:.4f}">{freq:.3f}</td>')
            cells.append(f'<td data-val="{cpm:.1f}">{cpm:,.0f}</td>')

        fold = row.get(fold_col, 0) if fold_col in row.index else 0
        l2fc = row.get(log2fc_col, 0) if log2fc_col in row.index else 0
        fold_str = 'de novo' if np.isinf(fold) else f'{fold:.1f}'
        cells.append(f'<td data-val="{fold if not np.isinf(fold) else 1e9}">{fold_str}</td>')
        cells.append(f'<td data-val="{l2fc:.2f}">{l2fc:.2f}</td>')
        cells.append(f'<td title="{cdr3}">{cdr3_display}</td>')
        cells.append(f'<td data-val="{cdr3_len}">{cdr3_len}</td>')
        cells.append(f'<td>{family}</td>')
        rows_html.append(f'<tr data-rank="{idx+1}" data-family="{family}">' + ''.join(cells) + '</tr>')

    th_html = ''.join(f'<th>{c}<span class="sort-arrow">&#x25B2;</span></th>' for c in header_cols)

    return f'''<div class="phylo-paginated">
    <div class="table-controls">
      <input type="text" class="phylo-search" placeholder="Search sequences...">
      <button class="btn-export" data-table="enrichment-table">Export CSV</button>
      <label style="font-size:12px;color:var(--text-dim);cursor:pointer;display:flex;align-items:center;gap:4px;margin-left:8px">
        <input type="checkbox" id="group-by-family">Group by family
      </label>
      <span class="phylo-entries-bar">Show
        <select class="phylo-page-size"><option value="10">10</option>
          <option value="25" selected>25</option><option value="50">50</option>
          <option value="100">100</option></select> entries</span>
    </div>
    <div style="overflow-x:auto">
    <table class="enrich-table" id="enrichment-table">
      <thead><tr>{th_html}</tr></thead>
      <tbody>{''.join(rows_html)}</tbody>
    </table>
    </div>
    <div class="phylo-page-controls">
      <div class="phylo-page-info"></div>
      <div class="phylo-page-nav"></div>
    </div>
    </div>'''


# ===========================================================================
# WEB LAUNCHER (--serve mode)
# ===========================================================================

LAUNCHER_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Nabbit — Launcher</title>
<link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500;600&display=swap" rel="stylesheet">
<style>
:root{
  --bg:#0c0e13;--surface:#14171e;--surface2:#1a1e28;--border:#252a36;
  --text:#e2e4e9;--text-dim:#7a8094;--accent:#a78bfa;
  --teal:#06b6d4;--orange:#f97316;--green:#34d399;--red:#f87171;
  --mono:'IBM Plex Mono',monospace;
}
[data-theme=light]{
  --bg:#f5f5f7;--surface:#ffffff;--surface2:#eeeef0;--border:#d4d4d8;
  --text:#1a1a2e;--text-dim:#6b7280;--accent:#7c3aed;
  --teal:#0891b2;--orange:#ea580c;--green:#059669;--red:#dc2626;
}
*{margin:0;padding:0;box-sizing:border-box}
*::-webkit-scrollbar{width:8px;height:8px}
*::-webkit-scrollbar-track{background:var(--surface);border-radius:4px}
*::-webkit-scrollbar-thumb{background:var(--border);border-radius:4px}
*::-webkit-scrollbar-thumb:hover{background:var(--text-dim)}
*::-webkit-scrollbar-corner{background:var(--surface)}
*{scrollbar-width:thin;scrollbar-color:var(--border) var(--surface)}
input[type=checkbox]{-webkit-appearance:none;appearance:none;width:16px;height:16px;border:1.5px solid var(--border);
  border-radius:4px;background:var(--surface2);cursor:pointer;position:relative;vertical-align:middle;
  flex-shrink:0;transition:all .15s}
input[type=checkbox]:hover{border-color:var(--accent)}
input[type=checkbox]:checked{background:var(--accent);border-color:var(--accent)}
input[type=checkbox]:checked::after{content:'';position:absolute;left:4.5px;top:1.5px;width:4px;height:8px;
  border:solid #fff;border-width:0 2px 2px 0;transform:rotate(45deg)}
select{-webkit-appearance:none;appearance:none;background:var(--surface2);border:1px solid var(--border);
  color:var(--text);padding:6px 28px 6px 10px;border-radius:6px;font-size:12px;font-family:var(--mono);
  cursor:pointer;transition:all .15s;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%237a8094'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:right 8px center}
select:hover{border-color:var(--accent)}
select:focus{outline:none;border-color:var(--accent)}
select option{background:var(--surface2);color:var(--text)}
body{background:var(--bg);color:var(--text);font-family:'DM Sans',sans-serif;line-height:1.5;min-height:100vh}
.container{max-width:900px;margin:0 auto;padding:32px 24px}
.header{margin-bottom:32px;display:flex;align-items:center;justify-content:space-between}
.header-left{display:flex;align-items:flex-start;flex-direction:column;gap:4px}
.logo{font-size:28px;font-weight:700;letter-spacing:-0.5px}
.logo span{color:var(--accent)}
.subtitle{font-size:14px;color:var(--text-dim);font-family:var(--mono)}
.theme-btn{background:var(--surface);border:1px solid var(--border);color:var(--text-dim);
  padding:8px 16px;border-radius:8px;cursor:pointer;font-size:12px;font-family:var(--mono);
  font-weight:500;transition:all .2s}
.theme-btn:hover{border-color:var(--accent);color:var(--text)}
.section{background:var(--surface);border:1px solid var(--border);border-radius:10px;padding:24px;margin-bottom:20px}
.section-title{font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:1.2px;
  color:var(--text-dim);margin-bottom:16px;font-family:var(--mono)}
.btn{background:var(--surface2);border:1px solid var(--border);color:var(--text-dim);
  padding:10px 20px;border-radius:8px;cursor:pointer;font-size:12px;font-family:var(--mono);
  font-weight:500;transition:all .2s;white-space:nowrap}
.btn:hover{border-color:var(--accent);color:var(--text)}
.btn-accent{background:var(--accent);border-color:var(--accent);color:#fff;font-weight:600;
  font-size:14px;padding:12px 32px}
.btn-accent:hover{opacity:.9;color:#fff}
.btn-accent:disabled{opacity:.4;cursor:default}
.btn-green{background:var(--green);border-color:var(--green);color:#fff;font-weight:600;
  font-size:14px;padding:12px 32px;text-decoration:none;display:inline-block;border-radius:8px;
  font-family:var(--mono);cursor:pointer;transition:all .2s}
.btn-green:hover{opacity:.9;color:#fff}
.drop-zone{border:2px dashed var(--border);border-radius:10px;padding:40px 24px;text-align:center;
  cursor:pointer;transition:all .2s;margin-bottom:16px}
.drop-zone:hover,.drop-zone.drag-over{border-color:var(--accent);background:rgba(167,139,250,0.05)}
.drop-zone-icon{margin-bottom:12px;opacity:.5}
.drop-zone-text{font-size:13px;color:var(--text-dim);font-family:var(--mono);margin-bottom:4px}
.drop-zone-hint{font-size:11px;color:var(--text-dim);font-family:var(--mono);opacity:.6}
.drop-zone-path{font-size:12px;color:var(--accent);font-family:var(--mono);margin-top:10px;
  word-break:break-all}
.manual-path{display:flex;gap:10px;align-items:center}
.manual-path input{background:var(--surface2);border:1px solid var(--border);color:var(--text);
  padding:8px 12px;border-radius:8px;font-size:12px;font-family:var(--mono);flex:1}
.manual-path input:focus{outline:none;border-color:var(--accent)}
table.pairs-table{width:100%;border-collapse:collapse;font-size:12px;font-family:var(--mono)}
table.pairs-table th{background:var(--surface2);color:var(--text-dim);padding:10px 8px;
  text-align:left;border-bottom:1px solid var(--border);font-size:11px;font-weight:600;
  text-transform:uppercase;letter-spacing:0.5px}
table.pairs-table td{padding:8px;border-bottom:1px solid var(--border);color:var(--text-dim)}
table.pairs-table select{font-size:12px}
.btn-remove{background:none;border:1px solid var(--border);color:var(--text-dim);width:26px;height:26px;
  border-radius:6px;cursor:pointer;font-size:16px;line-height:1;display:inline-flex;align-items:center;
  justify-content:center;transition:all 0.15s}
.btn-remove:hover{background:var(--red);color:#fff;border-color:var(--red)}
.options-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}
.option-item{display:flex;align-items:center;gap:8px}
.option-item label{font-size:12px;color:var(--text-dim);font-family:var(--mono);cursor:pointer}
.option-item input[type=checkbox]{margin-right:2px}
.option-item input[type=number]{width:70px;background:var(--surface2);border:1px solid var(--border);
  color:var(--text);padding:6px 8px;border-radius:6px;font-size:12px;font-family:var(--mono);
  -moz-appearance:textfield;appearance:textfield}
.option-item input[type=number]::-webkit-inner-spin-button,
.option-item input[type=number]::-webkit-outer-spin-button{-webkit-appearance:none;margin:0}
.num-wrap{position:relative;display:inline-flex;align-items:center}
.num-wrap input[type=number]{padding-right:24px}
.num-spin{position:absolute;right:1px;top:1px;bottom:1px;width:22px;display:flex;flex-direction:column;
  border-left:1px solid var(--border);border-radius:0 5px 5px 0;overflow:hidden}
.num-spin button{flex:1;background:var(--surface2);border:none;color:var(--text-dim);cursor:pointer;
  display:flex;align-items:center;justify-content:center;padding:0;font-size:10px;line-height:1}
.num-spin button:hover{background:var(--accent);color:#fff}
.num-spin button+button{border-top:1px solid var(--border)}
.status-area{margin-top:20px;display:none}
.status-bar{display:flex;align-items:center;gap:12px;margin-bottom:12px}
.spinner{width:20px;height:20px;border:2px solid var(--border);border-top-color:var(--accent);
  border-radius:50%;animation:spin 1s linear infinite}
@keyframes spin{to{transform:rotate(360deg)}}
.status-text{font-size:13px;font-family:var(--mono);color:var(--text-dim)}
.status-text.done{color:var(--green)}
.status-text.error{color:var(--red)}
.log-output{background:var(--surface2);border:1px solid var(--border);border-radius:8px;
  padding:12px;font-size:11px;font-family:var(--mono);color:var(--text-dim);
  max-height:300px;overflow-y:auto;white-space:pre-wrap;line-height:1.6}
.empty-state{text-align:center;padding:32px;color:var(--text-dim);font-size:12px;font-family:var(--mono)}
.run-bar{display:flex;justify-content:flex-end;align-items:center;gap:12px;margin-top:20px}
.custom-select{position:relative;display:inline-block}
.custom-select-trigger{background:var(--surface2);border:1px solid var(--border);color:var(--text);
  padding:6px 28px 6px 10px;border-radius:6px;font-size:12px;font-family:var(--mono);cursor:pointer;
  transition:all .15s;white-space:nowrap;user-select:none;position:relative}
.custom-select-trigger::after{content:'';position:absolute;right:9px;top:50%;transform:translateY(-50%);
  border-left:4px solid transparent;border-right:4px solid transparent;border-top:5px solid var(--text-dim);
  transition:transform .15s}
.custom-select.open .custom-select-trigger{border-color:var(--accent)}
.custom-select.open .custom-select-trigger::after{transform:translateY(-50%) rotate(180deg)}
.custom-select-trigger:hover{border-color:var(--accent)}
.custom-select-options{display:none;position:absolute;bottom:100%;left:0;margin-bottom:4px;
  background:var(--surface2);border:1px solid var(--border);border-radius:6px;min-width:100%;
  z-index:999;overflow:hidden;box-shadow:0 -4px 16px rgba(0,0,0,0.4)}
.custom-select.open .custom-select-options{display:block}
.custom-select.drop-down .custom-select-options{bottom:auto;top:100%;margin-bottom:0;margin-top:4px;
  box-shadow:0 4px 16px rgba(0,0,0,0.4)}
.custom-select-option{padding:6px 12px;font-size:12px;font-family:var(--mono);color:var(--text-dim);
  cursor:pointer;transition:all .1s;white-space:nowrap}
.custom-select-option:hover{background:var(--accent);color:#fff}
.custom-select-option.selected{color:var(--accent)}
.custom-select-option.selected:hover{color:#fff}
</style>
</head>
<body>
<div class="container">
  <div class="header">
    <div class="header-left">
      <div class="logo"><svg xmlns="http://www.w3.org/2000/svg" viewBox="100 90 300 280" width="40" height="40" style="vertical-align:middle;margin-right:10px"><defs><filter id="glow" x="-20%" y="-20%" width="140%" height="140%"><feGaussianBlur stdDeviation="5" result="blur"/><feComposite in="SourceGraphic" in2="blur" operator="over"/></filter></defs><ellipse cx="250" cy="330" rx="110" ry="35" fill="#2F3542"/><path d="M180 330C180 230,320 230,320 330Z" fill="#FF9F43"/><g transform="rotate(-20 200 160)"><ellipse cx="200" cy="160" rx="20" ry="65" fill="#FF9F43" stroke="#E17055" stroke-width="2"/><ellipse cx="200" cy="165" rx="10" ry="45" fill="#FFD1D1"/></g><g transform="rotate(25 300 160)"><ellipse cx="300" cy="160" rx="20" ry="65" fill="#FF9F43" stroke="#E17055" stroke-width="2"/><ellipse cx="300" cy="165" rx="10" ry="45" fill="#FFD1D1"/></g><ellipse cx="250" cy="245" rx="70" ry="60" fill="#FF9F43"/><ellipse cx="230" cy="265" rx="25" ry="20" fill="#FFF"/><ellipse cx="270" cy="265" rx="25" ry="20" fill="#FFF"/><path d="M240 255L260 255L250 265Z" fill="#FF4757"/><ellipse cx="220" cy="235" rx="8" ry="14" fill="#2F3542"/><circle cx="218" cy="230" r="3" fill="#FFF"/><ellipse cx="280" cy="235" rx="8" ry="14" fill="#2F3542"/><circle cx="278" cy="230" r="3" fill="#FFF"/><ellipse cx="185" cy="335" rx="22" ry="16" fill="#FF9F43" stroke="#E17055" stroke-width="1.5"/><path d="M175 345L175 330M195 345L195 330" stroke="#E17055" stroke-width="2" stroke-linecap="round"/><g filter="url(#glow)"><path d="M315 285L340 235L325 185M340 235L380 200" stroke="#1E90FF" stroke-width="18" stroke-linecap="round" stroke-linejoin="round" fill="none"/><path d="M315 285L340 235L325 185M340 235L380 200" stroke="#70A1FF" stroke-width="8" stroke-linecap="round" stroke-linejoin="round" fill="none"/><path d="M315 285L340 235L325 185M340 235L380 200" stroke="#FFF" stroke-width="3" stroke-linecap="round" stroke-linejoin="round" fill="none"/><path d="M380 160L383 168L391 171L383 174L380 182L377 174L369 171L377 168Z" fill="#70A1FF"/><circle cx="310" cy="180" r="4" fill="#70A1FF"/><circle cx="395" cy="220" r="3" fill="#70A1FF"/></g><ellipse cx="310" cy="285" rx="22" ry="18" fill="#FF9F43" stroke="#E17055" stroke-width="1.5" transform="rotate(-15 310 285)"/><path d="M305 295L315 275M320 295L330 275" stroke="#E17055" stroke-width="2" stroke-linecap="round"/></svg><span>Nabbit</span></div>
      <div class="subtitle">Nanobody NGS Enrichment Pipeline</div>
    </div>
    <button class="theme-btn" onclick="toggleTheme()">Toggle theme</button>
  </div>

  <div class="section">
    <div class="section-title">FASTQ Files</div>
    <div class="drop-zone" id="drop-zone" onclick="browseDir()">
      <div class="drop-zone-icon"><svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"><path d="M22 19a2 2 0 0 1-2 2H4a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h5l2 3h9a2 2 0 0 1 2 2z"/><line x1="12" y1="11" x2="12" y2="17"/><polyline points="9 14 12 11 15 14"/></svg></div>
      <div class="drop-zone-text">Click to browse for a FASTQ directory</div>
      <div class="drop-zone-hint">or drop .fastq.gz files here</div>
      <div class="drop-zone-path" id="drop-zone-path"></div>
    </div>
    <div class="manual-path">
      <input type="text" id="fastq-dir" placeholder="Or type a path manually...">
      <button class="btn" onclick="scanManual()">Scan</button>
    </div>
    <div id="scan-error" style="color:var(--red);font-size:12px;font-family:var(--mono);display:none;margin-top:12px"></div>
    <div id="pairs-container"></div>
  </div>

  <div class="section">
    <div class="section-title">Options</div>
    <div class="options-grid">
      <div class="option-item">
        <label for="opt-threads">Threads</label>
        <div class="num-wrap">
          <input type="number" id="opt-threads" value="4" min="1" max="64">
          <div class="num-spin">
            <button type="button" onclick="stepNum('opt-threads',1)">&#9650;</button>
            <button type="button" onclick="stepNum('opt-threads',-1)">&#9660;</button>
          </div>
        </div>
      </div>
      <div class="option-item">
        <label for="opt-topn">Top N</label>
        <div class="num-wrap">
          <input type="number" id="opt-topn" value="50" min="10" max="500">
          <div class="num-spin">
            <button type="button" onclick="stepNum('opt-topn',1)">&#9650;</button>
            <button type="button" onclick="stepNum('opt-topn',-1)">&#9660;</button>
          </div>
        </div>
      </div>
      <div class="option-item">
        <input type="checkbox" id="opt-ablang">
        <label for="opt-ablang">AbLang scoring</label>
      </div>
      <div class="option-item">
        <input type="checkbox" id="opt-igblast">
        <label for="opt-igblast">IgBLAST annotation</label>
      </div>
    </div>
  </div>

  <div class="run-bar">
    <button class="btn btn-accent" id="run-btn" onclick="runPipeline()" disabled>Run Nabbit</button>
  </div>

  <div class="status-area" id="status-area">
    <div class="section">
      <div class="section-title">Pipeline Status</div>
      <div class="status-bar">
        <div class="spinner" id="status-spinner"></div>
        <div class="status-text" id="status-text">Starting pipeline...</div>
      </div>
      <div id="dashboard-link-area" style="display:none;margin-bottom:12px;text-align:center">
        <a class="btn-green" id="dashboard-link" href="#" target="_blank">Open Dashboard</a>
      </div>
      <div class="log-output" id="log-output"></div>
    </div>
  </div>
</div>

<script>
let pairs = [];
let manualFiles = [];
let manualMode = false;
let selectedDir = '';
let currentRunId = null;
let pollTimer = null;

function initCustomSelects(){
  document.querySelectorAll('select').forEach(sel=>{
    if(sel.dataset.customized)return;
    sel.dataset.customized='1';
    sel.style.display='none';
    const wrap=document.createElement('div');
    wrap.className='custom-select';
    sel.parentNode.insertBefore(wrap,sel);
    wrap.appendChild(sel);
    const trigger=document.createElement('div');
    trigger.className='custom-select-trigger';
    const opts=document.createElement('div');
    opts.className='custom-select-options';
    function buildOpts(){
      opts.innerHTML='';
      Array.from(sel.options).forEach(o=>{
        const d=document.createElement('div');
        d.className='custom-select-option'+(o.selected?' selected':'');
        d.textContent=o.textContent;d.dataset.value=o.value;
        d.addEventListener('click',e=>{
          e.stopPropagation();
          sel.value=o.value;sel.dispatchEvent(new Event('change',{bubbles:true}));
          opts.querySelectorAll('.custom-select-option').forEach(x=>x.classList.remove('selected'));
          d.classList.add('selected');
          trigger.textContent=o.textContent;
          wrap.classList.remove('open');
        });
        opts.appendChild(d);
      });
      const cur=sel.options[sel.selectedIndex];
      trigger.textContent=cur?cur.textContent:'';
    }
    buildOpts();
    wrap.appendChild(trigger);wrap.appendChild(opts);
    trigger.addEventListener('click',e=>{
      e.stopPropagation();
      document.querySelectorAll('.custom-select.open').forEach(w=>{if(w!==wrap)w.classList.remove('open');});
      wrap.classList.toggle('open');
      if(wrap.classList.contains('open')){
        const r=wrap.getBoundingClientRect();
        const spaceBelow=window.innerHeight-r.bottom;
        wrap.classList.toggle('drop-down',spaceBelow>200);
      }
    });
  });
  document.addEventListener('click',()=>{document.querySelectorAll('.custom-select.open').forEach(w=>w.classList.remove('open'));});
}

function stepNum(id, dir) {
  const inp = document.getElementById(id);
  const min = Number(inp.min), max = Number(inp.max), step = Number(inp.step) || 1;
  let val = Number(inp.value) + dir * step;
  if (!isNaN(min) && val < min) val = min;
  if (!isNaN(max) && val > max) val = max;
  inp.value = val;
}

function toggleTheme() {
  const h = document.documentElement;
  h.setAttribute('data-theme', h.getAttribute('data-theme') === 'light' ? '' : 'light');
}

/* --- Drop zone --- */
const dropZone = document.getElementById('drop-zone');
dropZone.addEventListener('dragover', e => { e.preventDefault(); dropZone.classList.add('drag-over'); });
dropZone.addEventListener('dragleave', () => dropZone.classList.remove('drag-over'));
dropZone.addEventListener('drop', async e => {
  e.preventDefault();
  dropZone.classList.remove('drag-over');
  const files = Array.from(e.dataTransfer.files);
  if (files.length === 0) return;
  const fqFiles = files.filter(f =>
    /\\.(fastq|fq)(\\.gz)?$/i.test(f.name));
  if (fqFiles.length === 0) {
    showError('No FASTQ files found in the dropped items');
    return;
  }
  const droppedNames = fqFiles.map(f => f.name);
  document.getElementById('drop-zone-path').textContent = 'Locating files...';
  try {
    const resp = await fetch('/api/locate', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({filenames: droppedNames})
    });
    const data = await resp.json();
    if (data.dir) {
      selectedDir = data.dir;
      document.getElementById('fastq-dir').value = data.dir;
      document.getElementById('drop-zone-path').textContent = data.dir;
      scanDir(data.dir, droppedNames);
    } else {
      document.getElementById('drop-zone-path').textContent = '';
      showError(data.error || 'Could not locate files on disk.');
    }
  } catch (err) {
    document.getElementById('drop-zone-path').textContent = '';
    showError('Could not locate files automatically. Please enter the FASTQ directory path manually.');
  }
});

async function browseDir() {
  try {
    const resp = await fetch('/api/browse', {method: 'POST'});
    const data = await resp.json();
    if (data.dir) {
      selectedDir = data.dir;
      document.getElementById('fastq-dir').value = data.dir;
      document.getElementById('drop-zone-path').textContent = data.dir;
      scanDir(data.dir);
    }
  } catch (e) {
    showError('Browse failed: ' + e.message);
  }
}

function scanManual() {
  const dir = document.getElementById('fastq-dir').value.trim();
  if (!dir) return;
  selectedDir = dir;
  document.getElementById('drop-zone-path').textContent = dir;
  scanDir(dir);
}

function showError(msg) {
  const el = document.getElementById('scan-error');
  el.textContent = msg;
  el.style.display = 'block';
}

async function scanDir(dir, filterNames) {
  const errEl = document.getElementById('scan-error');
  errEl.style.display = 'none';
  try {
    const payload = {dir: dir};
    if (filterNames) payload.filenames = filterNames;
    const resp = await fetch('/api/scan', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(payload)
    });
    const data = await resp.json();
    if (data.error) { showError(data.error); return; }
    pairs = data.pairs || [];
    if (pairs.length > 0) {
      manualMode = false;
      manualFiles = [];
      renderPairs();
    } else if (data.files && data.files.length > 0) {
      /* Try client-side auto-pairing by common prefix */
      const pairRe = /^(.+)[_.]([12])\\.(?:fastq|fq)(?:\\.gz)?$/i;
      const r1r2Re = /^(.+?)[._]R([12])[._]/i;
      const groups = {};
      const ungrouped = [];
      data.files.forEach(f => {
        let m = pairRe.exec(f.name);
        if (m) { groups[m[1]] = groups[m[1]] || {}; groups[m[1]][m[2]] = f; return; }
        m = r1r2Re.exec(f.name);
        if (m) { groups[m[1]] = groups[m[1]] || {}; groups[m[1]][m[2]] = f; return; }
        ungrouped.push(f);
      });
      const autoPairs = [];
      Object.keys(groups).sort().forEach((prefix, idx) => {
        const g = groups[prefix];
        if (g['1'] && g['2']) {
          autoPairs.push({round: idx+1, r1: g['1'].path, r2: g['2'].path,
            r1_name: g['1'].name, r2_name: g['2'].name});
        } else {
          Object.values(g).forEach(f => ungrouped.push(f));
        }
      });
      if (autoPairs.length > 0 && ungrouped.length === 0) {
        pairs = autoPairs;
        manualMode = false;
        manualFiles = [];
        renderPairs();
      } else {
        /* Fall back to manual file assignment */
        manualMode = true;
        pairs = [];
        manualFiles = data.files.map((f, i) => ({
          path: f.path, name: f.name, round: Math.floor(i/2)+1,
          read: (i % 2 === 0) ? 'R1' : 'R2'
        }));
        renderManualFiles();
      }
    } else {
      manualMode = false;
      manualFiles = [];
      renderPairs();
    }
  } catch (e) { showError('Failed to scan directory: ' + e.message); }
}

function renderPairs() {
  const c = document.getElementById('pairs-container');
  if (pairs.length === 0) {
    c.innerHTML = '<div class="empty-state" style="margin-top:16px">No FASTQ pairs found in this directory</div>';
    document.getElementById('run-btn').disabled = true;
    return;
  }
  let html = '<table class="pairs-table" style="margin-top:16px"><thead><tr><th>Round</th><th>R1</th><th>R2</th><th></th></tr></thead><tbody>';
  pairs.forEach((p, i) => {
    let opts = '<option value="0"' + (p.round === 0 ? ' selected' : '') + '>Round 0 (naive)</option>';
    for (let r = 1; r <= 10; r++) {
      opts += '<option value="' + r + '"' + (r === p.round ? ' selected' : '') + '>Round ' + r + '</option>';
    }
    html += '<tr><td><select onchange="pairs[' + i + '].round=parseInt(this.value)">' + opts + '</select></td>';
    html += '<td title="' + p.r1 + '">' + p.r1_name + '</td>';
    html += '<td title="' + p.r2 + '">' + p.r2_name + '</td>';
    html += '<td><button class="btn-remove" onclick="removePair(' + i + ')" title="Remove this pair">&times;</button></td></tr>';
  });
  html += '</tbody></table>';
  c.innerHTML = html;
  document.getElementById('run-btn').disabled = false;
  initCustomSelects();
}

function removePair(index) {
  pairs.splice(index, 1);
  renderPairs();
}

function renderManualFiles() {
  const c = document.getElementById('pairs-container');
  if (manualFiles.length === 0) {
    c.innerHTML = '<div class="empty-state" style="margin-top:16px">No FASTQ files to assign</div>';
    document.getElementById('run-btn').disabled = true;
    return;
  }
  let html = '<div style="margin-top:12px;margin-bottom:8px;font-size:12px;color:var(--text-dim);font-family:var(--mono)">' +
    'Auto-pairing failed. Assign each file to a round and read manually:</div>';
  html += '<table class="pairs-table"><thead><tr><th>File</th><th>Round</th><th>Read</th><th></th></tr></thead><tbody>';
  manualFiles.forEach((f, i) => {
    let roundOpts = '<option value="0"' + (f.round === 0 ? ' selected' : '') + '>Round 0 (naive)</option>';
    for (let r = 1; r <= 10; r++) {
      roundOpts += '<option value="' + r + '"' + (r === f.round ? ' selected' : '') + '>Round ' + r + '</option>';
    }
    const r1Sel = f.read === 'R1' ? ' selected' : '';
    const r2Sel = f.read === 'R2' ? ' selected' : '';
    html += '<tr><td title="' + f.path + '">' + f.name + '</td>';
    html += '<td><select onchange="manualFiles[' + i + '].round=parseInt(this.value)">' + roundOpts + '</select></td>';
    html += '<td><select onchange="manualFiles[' + i + '].read=this.value">' +
      '<option value="R1"' + r1Sel + '>R1</option><option value="R2"' + r2Sel + '>R2</option></select></td>';
    html += '<td><button class="btn-remove" onclick="removeManualFile(' + i + ')" title="Remove this file">&times;</button></td></tr>';
  });
  html += '</tbody></table>';
  c.innerHTML = html;
  document.getElementById('run-btn').disabled = false;
  initCustomSelects();
}

function removeManualFile(index) {
  manualFiles.splice(index, 1);
  renderManualFiles();
}

async function runPipeline() {
  /* If in manual mode, validate and build pairs from manualFiles */
  if (manualMode) {
    const grouped = {};
    manualFiles.forEach(f => {
      const key = f.round;
      if (!grouped[key]) grouped[key] = {};
      if (grouped[key][f.read]) {
        showError('Round ' + key + ' has duplicate ' + f.read + ' assignments');
        return;
      }
      grouped[key][f.read] = f;
    });
    const built = [];
    const errors = [];
    for (const rnd of Object.keys(grouped).map(Number).sort((a,b)=>a-b)) {
      const g = grouped[rnd];
      if (!g.R1) errors.push('Round ' + rnd + ' is missing R1');
      if (!g.R2) errors.push('Round ' + rnd + ' is missing R2');
      if (g.R1 && g.R2) {
        built.push({round: rnd, r1: g.R1.path, r2: g.R2.path,
          r1_name: g.R1.name, r2_name: g.R2.name});
      }
    }
    if (errors.length > 0) { showError(errors.join('; ')); return; }
    if (built.length === 0) { showError('No valid pairs could be built from assignments'); return; }
    pairs = built;
  }
  if (pairs.length === 0) return;
  const btn = document.getElementById('run-btn');
  btn.disabled = true;
  const statusArea = document.getElementById('status-area');
  statusArea.style.display = 'block';
  document.getElementById('status-spinner').style.display = 'block';
  document.getElementById('status-text').className = 'status-text';
  document.getElementById('status-text').textContent = 'Starting pipeline...';
  document.getElementById('dashboard-link-area').style.display = 'none';
  document.getElementById('log-output').textContent = '';

  const config = {
    dir: selectedDir || document.getElementById('fastq-dir').value.trim(),
    pairs: pairs,
    threads: parseInt(document.getElementById('opt-threads').value) || 4,
    top_n: parseInt(document.getElementById('opt-topn').value) || 50,
    ablang: document.getElementById('opt-ablang').checked,
    igblast: document.getElementById('opt-igblast').checked,
  };

  try {
    const resp = await fetch('/api/run', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(config)
    });
    const data = await resp.json();
    if (data.error) {
      document.getElementById('status-text').textContent = 'Error: ' + data.error;
      document.getElementById('status-text').className = 'status-text error';
      document.getElementById('status-spinner').style.display = 'none';
      btn.disabled = false;
      return;
    }
    currentRunId = data.run_id;
    pollTimer = setInterval(pollStatus, 2000);
  } catch (e) {
    document.getElementById('status-text').textContent = 'Failed to start: ' + e.message;
    document.getElementById('status-text').className = 'status-text error';
    document.getElementById('status-spinner').style.display = 'none';
    btn.disabled = false;
  }
}

async function pollStatus() {
  if (!currentRunId) return;
  try {
    const resp = await fetch('/api/status/' + currentRunId);
    const data = await resp.json();
    const logEl = document.getElementById('log-output');
    if (data.log) {
      logEl.textContent = data.log;
      logEl.scrollTop = logEl.scrollHeight;
    }
    if (data.status === 'done') {
      clearInterval(pollTimer);
      document.getElementById('status-spinner').style.display = 'none';
      document.getElementById('status-text').textContent = 'Pipeline complete!';
      document.getElementById('status-text').className = 'status-text done';
      if (data.dashboard) {
        const link = document.getElementById('dashboard-link');
        link.href = data.dashboard;
        document.getElementById('dashboard-link-area').style.display = 'block';
      }
      document.getElementById('run-btn').disabled = false;
    } else if (data.status === 'error') {
      clearInterval(pollTimer);
      document.getElementById('status-spinner').style.display = 'none';
      document.getElementById('status-text').textContent = 'Pipeline failed (exit code ' + (data.exit_code || '?') + ')';
      document.getElementById('status-text').className = 'status-text error';
      document.getElementById('run-btn').disabled = false;
    }
  } catch (e) { /* ignore transient fetch errors */ }
}
</script>
</body>
</html>"""

_RUNS = {}  # run_id -> {process, output_dir, log_file, log_fh, closed}


class NabbitHandler:
    """HTTP request handler for the Nabbit web launcher."""

    # Import here so the class definition doesn't fail if http.server is
    # somehow absent, though it's stdlib.
    from http.server import BaseHTTPRequestHandler as _Base

    class _Handler(_Base):

        def do_GET(self):
            if self.path == '/' or self.path == '':
                self._send_html(LAUNCHER_HTML)
            elif self.path.startswith('/api/status/'):
                self._handle_status()
            elif self.path.startswith('/results/'):
                self._serve_result_file()
            else:
                self._send_json({'error': 'Not found'}, 404)

        def do_POST(self):
            if self.path == '/api/scan':
                self._handle_scan()
            elif self.path == '/api/run':
                self._handle_run()
            elif self.path == '/api/browse':
                self._handle_browse()
            elif self.path == '/api/locate':
                self._handle_locate()
            else:
                self._send_json({'error': 'Not found'}, 404)

        def _read_body(self):
            length = int(self.headers.get('Content-Length', 0))
            return json.loads(self.rfile.read(length)) if length else {}

        def _send_json(self, data, code=200):
            body = json.dumps(data).encode()
            self.send_response(code)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _send_html(self, html):
            body = html.encode()
            self.send_response(200)
            self.send_header('Content-Type', 'text/html')
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        # --- /api/browse ---
        def _handle_browse(self):
            """Open a native OS directory picker and return the selected path."""
            chosen = None
            try:
                if sys.platform == 'darwin':
                    r = subprocess.run(
                        ['osascript', '-e',
                         'POSIX path of (choose folder with prompt "Select FASTQ Directory")'],
                        capture_output=True, text=True, timeout=120)
                    if r.returncode == 0 and r.stdout.strip():
                        chosen = r.stdout.strip().rstrip('/')
                else:
                    # Linux / other: try tkinter
                    r = subprocess.run(
                        [sys.executable, '-c',
                         'import tkinter as tk; from tkinter import filedialog; '
                         'root=tk.Tk(); root.withdraw(); '
                         'print(filedialog.askdirectory(title="Select FASTQ Directory")); '
                         'root.destroy()'],
                        capture_output=True, text=True, timeout=120)
                    if r.returncode == 0 and r.stdout.strip():
                        chosen = r.stdout.strip()
            except Exception:
                pass
            self._send_json({'dir': chosen})

        # --- /api/locate (find dropped files on disk) ---
        def _handle_locate(self):
            try:
                body = self._read_body()
            except Exception as e:
                self._send_json({'error': f'Bad request: {e}'}, 400)
                return
            filenames = body.get('filenames', [])
            if not filenames:
                self._send_json({'error': 'No filenames provided'})
                return
            target = filenames[0]
            found_dir = None
            try:
                if sys.platform == 'darwin':
                    r = subprocess.run(
                        ['mdfind', '-name', target],
                        capture_output=True, text=True, timeout=15)
                    if r.returncode == 0:
                        best_dir, best_hits = None, 0
                        for line in r.stdout.strip().split('\n'):
                            if not line:
                                continue
                            candidate = os.path.dirname(line)
                            hits = sum(1 for f in filenames
                                       if os.path.isfile(os.path.join(candidate, f)))
                            if hits == len(filenames):
                                found_dir = candidate
                                break
                            if hits > best_hits:
                                best_dir, best_hits = candidate, hits
                        else:
                            if best_dir and best_hits >= len(filenames) // 2:
                                found_dir = best_dir
                else:
                    # Linux fallback: search common locations
                    home = os.path.expanduser('~')
                    for root_dir in [home, '/data', '/tmp']:
                        if not os.path.isdir(root_dir):
                            continue
                        r = subprocess.run(
                            ['find', root_dir, '-maxdepth', '5',
                             '-name', target, '-type', 'f'],
                            capture_output=True, text=True, timeout=15)
                        if r.returncode == 0 and r.stdout.strip():
                            first_hit = r.stdout.strip().split('\n')[0]
                            found_dir = os.path.dirname(first_hit)
                            break
            except subprocess.TimeoutExpired:
                self._send_json({'dir': None, 'error':
                    'File search timed out. Please enter the path manually.'})
                return
            except Exception as e:
                self._send_json({'dir': None, 'error':
                    f'File search failed: {e}. Please enter the path manually.'})
                return
            if found_dir:
                self._send_json({'dir': found_dir})
            else:
                self._send_json({'dir': None, 'error':
                    'Could not locate files on disk. Please use Browse or enter the path manually.'})

        # --- /api/scan ---
        def _handle_scan(self):
            body = self._read_body()
            dir_path = body.get('dir', '')
            filter_names = set(body.get('filenames', []))
            if not dir_path or not os.path.isdir(dir_path):
                self._send_json({'error': f'Directory not found: {dir_path}'})
                return
            try:
                fastq_map = discover_fastqs(dir_path)
                pairs = []
                for rnd in sorted(fastq_map):
                    r1 = fastq_map[rnd]['R1']
                    r2 = fastq_map[rnd]['R2']
                    # If filenames were provided (drag-drop), only include
                    # pairs where both files are in the dropped set
                    if filter_names:
                        if os.path.basename(r1) not in filter_names or \
                           os.path.basename(r2) not in filter_names:
                            continue
                    pairs.append({
                        'round': rnd,
                        'r1': r1,
                        'r2': r2,
                        'r1_name': os.path.basename(r1),
                        'r2_name': os.path.basename(r2),
                    })
                # If no pairs found, return individual FASTQ files for manual assignment
                raw_files = []
                if not pairs:
                    dp = Path(dir_path)
                    for ext in ('*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq'):
                        for fp in dp.rglob(ext):
                            name = fp.name
                            if filter_names and name not in filter_names:
                                continue
                            raw_files.append({'path': str(fp.resolve()), 'name': name})
                    raw_files.sort(key=lambda x: x['name'])
                self._send_json({'pairs': pairs, 'files': raw_files})
            except Exception as e:
                self._send_json({'error': str(e)})

        # --- /api/run ---
        def _handle_run(self):
            body = self._read_body()
            pairs = body.get('pairs', [])
            if not pairs:
                self._send_json({'error': 'No FASTQ pairs configured'})
                return

            run_id = f"{int(time.time())}_{uuid.uuid4().hex[:6]}"
            output_dir = os.path.join(tempfile.gettempdir(), f'nabbit_{run_id}')
            os.makedirs(output_dir, exist_ok=True)

            # Write sample map so the pipeline knows the round assignments
            sample_map = os.path.join(output_dir, 'sample_map.tsv')
            with open(sample_map, 'w') as f:
                f.write('# round\tR1\tR2\n')
                for p in pairs:
                    f.write(f"{p['round']}\t{p['r1']}\t{p['r2']}\n")

            cmd = [
                sys.executable, os.path.abspath(__file__),
                '--fastq-dir', body.get('dir', '.'),
                '--output-dir', output_dir,
                '--sample-map', sample_map,
                '--threads', str(body.get('threads', 4)),
                '--top-n', str(body.get('top_n', 50)),
            ]
            if body.get('ablang'):
                cmd.append('--ablang')
            if body.get('igblast'):
                cmd.append('--igblast')

            log_file = os.path.join(output_dir, 'nabbit.log')
            log_fh = open(log_file, 'w')
            proc = subprocess.Popen(cmd, stdout=log_fh, stderr=subprocess.STDOUT)

            _RUNS[run_id] = {
                'process': proc,
                'output_dir': output_dir,
                'log_file': log_file,
                'log_fh': log_fh,
                'closed': False,
            }
            self._send_json({'run_id': run_id, 'output_dir': output_dir})

        # --- /api/status/<id> ---
        def _handle_status(self):
            run_id = self.path.split('/api/status/')[-1]
            run = _RUNS.get(run_id)
            if not run:
                self._send_json({'error': 'Run not found'}, 404)
                return

            proc = run['process']
            log_text = ''
            try:
                with open(run['log_file'], 'r') as f:
                    f.seek(0, 2)
                    size = f.tell()
                    f.seek(max(0, size - 51200))  # last ~50 KB
                    log_text = f.read()
            except Exception:
                pass

            poll = proc.poll()
            if poll is None:
                self._send_json({'status': 'running', 'log': log_text})
            else:
                if not run['closed']:
                    run['log_fh'].close()
                    run['closed'] = True
                if poll == 0:
                    dashboard_path = os.path.join(run['output_dir'], 'dashboard.html')
                    dashboard_url = f'/results/{run_id}/dashboard.html' if os.path.isfile(dashboard_path) else None
                    self._send_json({
                        'status': 'done', 'log': log_text,
                        'dashboard': dashboard_url,
                        'output_dir': run['output_dir'],
                    })
                else:
                    self._send_json({'status': 'error', 'exit_code': poll, 'log': log_text})

        # --- /results/<id>/... ---
        def _serve_result_file(self):
            parts = self.path.split('/', 3)  # ['', 'results', run_id, filename]
            if len(parts) < 4:
                self._send_json({'error': 'Not found'}, 404)
                return
            run_id = parts[2]
            filename = parts[3]
            run = _RUNS.get(run_id)
            if not run:
                self._send_json({'error': 'Run not found'}, 404)
                return

            filepath = os.path.join(run['output_dir'], filename)
            # Prevent path traversal
            real_path = os.path.realpath(filepath)
            real_outdir = os.path.realpath(run['output_dir'])
            if not real_path.startswith(real_outdir + os.sep) and real_path != real_outdir:
                self._send_json({'error': 'Access denied'}, 403)
                return
            if not os.path.isfile(filepath):
                self._send_json({'error': 'File not found'}, 404)
                return

            ext = os.path.splitext(filepath)[1].lower()
            content_types = {
                '.html': 'text/html', '.css': 'text/css', '.js': 'application/javascript',
                '.json': 'application/json', '.tsv': 'text/tab-separated-values',
                '.csv': 'text/csv', '.png': 'image/png', '.svg': 'image/svg+xml',
            }
            ctype = content_types.get(ext, 'application/octet-stream')
            with open(filepath, 'rb') as f:
                data = f.read()
            self.send_response(200)
            self.send_header('Content-Type', ctype)
            self.send_header('Content-Length', str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def log_message(self, format, *args):
            pass  # suppress default request logging


def start_server(port=8080):
    """Start the Nabbit web launcher on localhost."""
    from http.server import HTTPServer
    import webbrowser
    handler = NabbitHandler._Handler
    server = HTTPServer(('', port), handler)
    log.info(f"Nabbit launcher running at http://localhost:{port}")
    webbrowser.open(f'http://localhost:{port}')
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        log.info("Shutting down server")
        server.shutdown()


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Nabbit — Nanobody NGS Enrichment Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--fastq-dir', default=None)
    parser.add_argument('--output-dir', default=None)
    parser.add_argument('--sample-map', default=None)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--min-overlap', type=int, default=20)
    parser.add_argument('--min-protein-len', type=int, default=90)
    parser.add_argument('--max-protein-len', type=int, default=180)
    parser.add_argument('--min-avg-qual', type=float, default=20.0)
    parser.add_argument('--top-n', type=int, default=50)
    parser.add_argument('--cdr3-cluster-threshold', type=float, default=0.8)
    parser.add_argument('--rounds', default=None, help='Comma-separated round numbers')
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--no-dashboard', action='store_true', help='Skip interactive HTML dashboard')
    parser.add_argument('--ablang', action='store_true', help='Score with AbLang (pip install ablang)')
    parser.add_argument('--ablang-device', default='cpu', help='cpu or cuda')
    parser.add_argument('--igblast', action='store_true', help='Enable IgBLAST V/D/J annotation')
    parser.add_argument('--igblast-db', default='./igblast_refs',
                        help='IgBLAST database directory (default: ./igblast_refs)')
    parser.add_argument('--igblast-chunk-size', type=int, default=0,
                        help='Chunk size for IgBLAST (0=no chunking, recommended: 25000)')
    parser.add_argument('--serve', action='store_true', help='Start web launcher UI')
    parser.add_argument('--port', type=int, default=8080, help='Port for web server (default: 8080)')

    args = parser.parse_args()

    if args.serve:
        start_server(args.port)
        return

    if not args.fastq_dir or not args.output_dir:
        parser.error('--fastq-dir and --output-dir are required (unless using --serve)')

    rounds_list = [int(x) for x in args.rounds.split(',')] if args.rounds else None
    os.makedirs(args.output_dir, exist_ok=True)

    # --- 1. Discover FASTQs ---
    log.info("=" * 60)
    log.info("STEP 1: Discovering FASTQ files")
    fastq_map = discover_fastqs(args.fastq_dir, args.sample_map, rounds_list)
    if not fastq_map:
        log.error("No valid FASTQ pairs found!"); sys.exit(1)
    for rnd in sorted(fastq_map):
        log.info(f"  Round {rnd}: {os.path.basename(fastq_map[rnd]['R1'])}")

    # --- 2. Process reads ---
    log.info("=" * 60)
    log.info("STEP 2: Processing reads")
    rounds_data = {}
    all_stats = {}
    tasks = [(rnd, fastq_map[rnd]['R1'], fastq_map[rnd]['R2'],
              args.min_overlap, args.min_protein_len, args.max_protein_len, args.min_avg_qual)
             for rnd in sorted(fastq_map)]

    if args.threads > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            for future in as_completed({ex.submit(_process_round, t): t[0] for t in tasks}):
                rnd, counter, stats = future.result()
                rounds_data[f"Round{rnd}"] = counter
                all_stats[f"Round{rnd}"] = stats
                log.info(f"  Round {rnd}: {stats['valid_proteins']:,} proteins ({len(counter):,} unique)")
    else:
        for t in tasks:
            rnd, counter, stats = _process_round(t)
            rounds_data[f"Round{rnd}"] = counter
            all_stats[f"Round{rnd}"] = stats
            log.info(f"  Round {rnd}: {stats['valid_proteins']:,} proteins ({len(counter):,} unique)")

    pd.DataFrame(all_stats).T.to_csv(os.path.join(args.output_dir, 'processing_stats.tsv'), sep='\t')

    # --- 3. Enrichment ---
    log.info("=" * 60)
    log.info("STEP 3: Enrichment analysis")
    full_df, top_df = compute_enrichment(rounds_data, args.top_n)
    top_df.to_csv(os.path.join(args.output_dir, 'enrichment_table.tsv'), sep='\t', index=False)
    full_df.to_csv(os.path.join(args.output_dir, 'enrichment_full.tsv.gz'), sep='\t', index=False, compression='gzip')

    # --- 4. CDR annotation ---
    log.info("=" * 60)
    log.info("STEP 4: CDR annotation")
    annotator = VHHAnnotator()
    ann_rows = []
    for i, row in top_df.iterrows():
        a = annotator.annotate(row['Sequence'])
        a['Sequence'] = row['Sequence']
        a['Rank'] = i + 1
        a['Length'] = int(row['Length'])
        ann_rows.append(a)
    annotations = pd.DataFrame(ann_rows)
    good = (annotations['annotation_quality'] == 'good').sum()
    log.info(f"  {good}/{len(annotations)} fully annotated")

    # --- 5. Clone clustering ---
    log.info("=" * 60)
    log.info("STEP 5: Clone clustering")
    annotations = cluster_by_cdr3(annotations, args.cdr3_cluster_threshold)
    annotations = detect_chimeras(annotations)
    n_chimeras = annotations['is_potential_chimera'].sum()
    log.info(f"  {annotations['clone_id'].nunique()} clones, {annotations['clone_family'].nunique()} families")
    if n_chimeras:
        log.info(f"  ⚠ {n_chimeras} potential chimeras")
    annotations.to_csv(os.path.join(args.output_dir, 'cdr_annotations.tsv'), sep='\t', index=False)

    # --- 6. Cross-round tracking ---
    log.info("=" * 60)
    log.info("STEP 6: Clone tracking")
    clone_tracking = track_clones(rounds_data, annotations, args.top_n)
    clone_tracking.to_csv(os.path.join(args.output_dir, 'clone_tracking.tsv'), sep='\t', index=False)
    log.info(f"  Tracking {len(clone_tracking)} families")

    # --- 7. AbLang scoring ---
    ablang_df = None
    if args.ablang:
        log.info("=" * 60)
        log.info("STEP 7: AbLang scoring")
        ablang_df = score_with_ablang(top_df['Sequence'].tolist(), device=args.ablang_device)
        if not ablang_df.empty:
            ablang_df.to_csv(os.path.join(args.output_dir, 'ablang_scores.tsv'), sep='\t', index=False)
            log.info(f"  Scored {len(ablang_df)} sequences")

    # --- 8. Diversity ---
    log.info("=" * 60)
    log.info("STEP 8: Diversity & convergence")
    diversity_df = compute_diversity(rounds_data)
    diversity_df.to_csv(os.path.join(args.output_dir, 'diversity_stats.tsv'), sep='\t', index=False)
    convergence_df = compute_convergence(rounds_data)
    convergence_df.to_csv(os.path.join(args.output_dir, 'convergence_analysis.tsv'), sep='\t', index=False)

    # --- 8b. Saturation analysis ---
    log.info("=" * 60)
    log.info("STEP 8b: Saturation analysis")
    saturation_df = compute_saturation(rounds_data, annotations)
    saturation_df.to_csv(os.path.join(args.output_dir, 'saturation_data.tsv'), sep='\t', index=False)
    diversity_est = estimate_diversity(saturation_df)
    diversity_est.to_csv(os.path.join(args.output_dir, 'diversity_estimates.tsv'), sep='\t', index=False)
    log.info(f"  Saturation: {len(saturation_df)} data points")
    if not diversity_est.empty:
        for _, row in diversity_est.iterrows():
            log.info(f"  {row['Round']}: ~{int(row['approx_diversity']):,} unique CDR3s")

    # --- 8c. Productivity stats ---
    log.info("=" * 60)
    log.info("STEP 8c: Productivity stats")
    productivity_df = compute_productivity_stats(all_stats)
    productivity_df.to_csv(os.path.join(args.output_dir, 'productivity_stats.tsv'), sep='\t', index=False)
    log.info(f"  Productivity stats for {len(productivity_df)} rounds")

    # --- 8d. IgBLAST annotation (optional) ---
    igblast_df = None
    gene_usage_df = None
    if args.igblast:
        log.info("=" * 60)
        log.info("STEP 8d: IgBLAST V/D/J annotation")
        # Build sequence list from top_df with annotation headers
        ig_seqs = []
        for i, row in top_df.iterrows():
            ann_row = annotations[annotations['Sequence'] == row['Sequence']]
            cid = ann_row.iloc[0]['clone_id'] if not ann_row.empty else f'seq_{i}'
            ig_seqs.append((cid, row['Sequence']))
        igblast_df, gene_usage_df = run_igblast(ig_seqs, args.igblast_db, args.threads, args.igblast_chunk_size)
        if igblast_df is not None and not igblast_df.empty:
            igblast_df.to_csv(os.path.join(args.output_dir, 'igblast_annotations.tsv'), sep='\t', index=False)
            log.info(f"  Annotated {len(igblast_df)} sequences")
        if gene_usage_df is not None and not gene_usage_df.empty:
            gene_usage_df.to_csv(os.path.join(args.output_dir, 'vdj_gene_usage.tsv'), sep='\t', index=False)
            log.info(f"  Gene usage: {len(gene_usage_df)} entries")

    # --- 8e. Cluster PCA (optional) ---
    cluster_pca = None
    try:
        log.info("=" * 60)
        log.info("STEP 8e: Cluster PCA")
        cluster_pca = compute_cluster_pca(top_df, annotations)
        if cluster_pca is not None:
            log.info(f"  PCA: {len(cluster_pca)} sequences projected")
    except Exception as e:
        log.warning(f"Cluster PCA failed: {e}")

    # --- 8f. CDR3 phylogenetic tree (optional) ---
    try:
        log.info("=" * 60)
        log.info("STEP 8f: CDR3 phylogenetic tree")
        build_cdr3_tree(annotations, args.output_dir)
    except Exception as e:
        log.warning(f"CDR3 tree failed: {e}")

    # --- 9. Outputs ---
    log.info("=" * 60)
    log.info("STEP 9: Writing outputs")
    fasta_entries = []
    for i, row in top_df.iterrows():
        ann_row = annotations[annotations['Sequence'] == row['Sequence']]
        cid = ann_row.iloc[0]['clone_id'] if not ann_row.empty else 'unk'
        cdr3_len = int(ann_row.iloc[0]['CDR3_length']) if not ann_row.empty else 0
        last = sorted(rounds_data.keys(), key=lambda x: int(x.replace('Round', '')))[-1]
        fasta_entries.append((
            f"rank{i+1}_{cid}_len{int(row['Length'])}_CDR3len{cdr3_len}",
            row['Sequence']))
    write_fasta(fasta_entries, os.path.join(args.output_dir, 'top_enriched_sequences.fasta'))

    # Write diverse sequences FASTA
    diverse_entries = select_diverse_sequences(top_df, annotations, n=args.top_n)
    write_fasta(diverse_entries, os.path.join(args.output_dir, 'top_diverse_sequences.fasta'))
    log.info(f"  Diverse sequences: {len(diverse_entries)}")

    write_report(top_df, annotations, clone_tracking, diversity_df, rounds_data, args.output_dir,
                 ablang_df=ablang_df, productivity_df=productivity_df,
                 saturation_df=saturation_df, diversity_est=diversity_est,
                 igblast_df=igblast_df)
    if not args.no_plots:
        generate_plots(top_df, diversity_df, rounds_data, annotations, args.output_dir,
                       ablang_df=ablang_df, clone_tracking=clone_tracking,
                       all_stats=all_stats, saturation_df=saturation_df,
                       diversity_est=diversity_est, igblast_df=igblast_df,
                       gene_usage_df=gene_usage_df, cluster_pca=cluster_pca)

    # --- 9b. Interactive Dashboard ---
    if not args.no_dashboard:
        try:
            # Build FASTA strings for dashboard download buttons
            fasta_enrich_str = '\n'.join(f'>{h}\n{s}' for h, s in fasta_entries)
            fasta_diverse_str = '\n'.join(f'>{h}\n{s}' for h, s in diverse_entries)
            dashboard_path = generate_dashboard(
                top_df=top_df, diversity_df=diversity_df, rounds_data=rounds_data,
                annotations=annotations, output_dir=args.output_dir,
                all_stats=all_stats, clone_tracking=clone_tracking,
                saturation_df=saturation_df, diversity_est=diversity_est,
                productivity_df=productivity_df, convergence_df=convergence_df,
                ablang_df=ablang_df, igblast_df=igblast_df,
                gene_usage_df=gene_usage_df, cluster_pca=cluster_pca,
                fasta_enrichment=fasta_enrich_str,
                fasta_diversity=fasta_diverse_str,
            )
            log.info(f"  Dashboard: {dashboard_path}")
        except Exception as e:
            log.warning(f"Dashboard generation failed: {e}")

    log.info("=" * 60)
    log.info(f"DONE — results in {args.output_dir}")


if __name__ == '__main__':
    main()
