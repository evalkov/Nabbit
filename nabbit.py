#!/usr/bin/env python3
"""
Nabbit — Nanobody NGS Enrichment Analysis Pipeline
=================================================
HPC pipeline for analyzing nanobody phage/yeast display enrichment
from paired-end Illumina sequencing data (fastq.gz).

Features:
  - Paired-end read merging, translation, quality filtering
  - Enrichment analysis with fold-change and diversity metrics
  - CDR annotation (IMGT-based VHH motif anchoring)
  - Clone clustering by CDR3 identity/similarity
  - Cross-round clone tracking with chimera detection
  - AbLang pseudo-perplexity scoring (optional, --ablang)

Dependencies:
  module load python/3.11
  pip install --user numpy pandas scipy matplotlib
  pip install --user ablang   # optional, for --ablang flag

Usage:
  python nanobody_ngs_enrichment.py \\
      --fastq-dir /data/$USER/nanobody_sequencing \\
      --output-dir /data/$USER/results \\
      --threads 8 --ablang
"""

import argparse
import gzip
import logging
import os
import re
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
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
            row[f'{rnd}_count'] = c
            row[f'{rnd}_freq'] = c / total if total > 0 else 0
        f1 = row.get(f'{round_names[0]}_freq', 0)
        fl = row.get(f'{round_names[-1]}_freq', 0)
        row['enrichment'] = fl / f1 if f1 > 0 else (float('inf') if fl > 0 else 0)
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

    # Fold-changes
    first = round_names[0]
    ps = 1e-8
    df[f'{first}_to_{last_round}_fold'] = (df[f'{last_round}_freq'] + ps) / (df[f'{first}_freq'] + ps)
    for i in range(len(round_names) - 1):
        rp, rn = round_names[i], round_names[i+1]
        df[f'{rp}_to_{rn}_fold'] = (df[f'{rn}_freq'] + ps) / (df[f'{rp}_freq'] + ps)

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
# Plots
# ===========================================================================
def generate_plots(top_df, diversity_df, rounds_data, annotations, output_dir):
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

    log.info(f"Plots saved to {output_dir}")


# ===========================================================================
# Report
# ===========================================================================
def write_report(top_df, annotations, clone_tracking, diversity_df, rounds_data, output_dir,
                 ablang_df=None):
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
# MAIN
# ===========================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Nabbit — Nanobody NGS Enrichment Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--fastq-dir', required=True)
    parser.add_argument('--output-dir', required=True)
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
    parser.add_argument('--ablang', action='store_true', help='Score with AbLang (pip install ablang)')
    parser.add_argument('--ablang-device', default='cpu', help='cpu or cuda')

    args = parser.parse_args()
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
    compute_convergence(rounds_data).to_csv(os.path.join(args.output_dir, 'convergence_analysis.tsv'), sep='\t', index=False)

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
    write_report(top_df, annotations, clone_tracking, diversity_df, rounds_data, args.output_dir, ablang_df)
    if not args.no_plots:
        generate_plots(top_df, diversity_df, rounds_data, annotations, args.output_dir)

    log.info("=" * 60)
    log.info(f"DONE — results in {args.output_dir}")


if __name__ == '__main__':
    main()
