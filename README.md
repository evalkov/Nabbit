# Nanobody NGS Enrichment Pipeline v2.0

Analyze nanobody phage/yeast display selection from paired-end Illumina fastq.gz.

## What it does

1. Merges paired reads (R1+R2 overlap), translates, quality-filters
2. Counts unique protein sequences per round, computes enrichment
3. Annotates CDR1/2/3 and FR1-4 regions (IMGT-based VHH motif anchoring)
4. Clusters sequences by CDR3 identity into clonotypes and families
5. Tracks clone families across rounds, flags chimeras
6. Optionally scores with AbLang (antibody language model pseudo-perplexity)
7. Outputs TSVs, FASTA, plots, human-readable report

## Dependencies

```bash
module load python/3.11
# numpy, pandas, scipy, matplotlib should be available
# For AbLang scoring (optional):
pip install --user ablang
```

## Usage

```bash
# Basic
sbatch run_pipeline.sh /data/$USER/fastqs /data/$USER/results

# Specific rounds
sbatch run_pipeline.sh /data/$USER/fastqs /data/$USER/results "1,2,3,5,6"

# With AbLang scoring
sbatch run_pipeline.sh /data/$USER/fastqs /data/$USER/results "" --ablang

# Interactive
sinteractive --cpus-per-task=8 --mem=32g
module load python/3.11
python nanobody_ngs_enrichment.py \
    --fastq-dir /data/$USER/fastqs \
    --output-dir /data/$USER/results \
    --threads 8 --ablang
```

## FASTQ naming

Expects `{round_number}_R{1,2}_001.fastq.gz` or Illumina-style naming.
Or use `--sample-map map.tsv` (tab-separated: `round  R1_path  R2_path`).

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--fastq-dir` | required | Directory with fastq.gz files |
| `--output-dir` | required | Output directory |
| `--threads` | 1 | Parallel workers |
| `--min-overlap` | 20 | Min bp overlap for merging |
| `--min-protein-len` | 90 | Min protein length (aa) |
| `--max-protein-len` | 180 | Max protein length (aa) |
| `--min-avg-qual` | 20 | Min average Phred quality |
| `--top-n` | 50 | Top sequences to report |
| `--cdr3-cluster-threshold` | 0.8 | CDR3 identity for family grouping |
| `--rounds` | all | Comma-separated rounds |
| `--no-plots` | off | Skip matplotlib |
| `--ablang` | off | Score with AbLang heavy-chain model |
| `--ablang-device` | cpu | cpu or cuda |

## Outputs

```
results/
├── processing_stats.tsv           # Read counts, merge rates per round
├── enrichment_table.tsv           # Top N with counts/freq/fold-change
├── enrichment_full.tsv.gz         # All sequences (gzipped)
├── cdr_annotations.tsv            # FR1-4, CDR1-3, clone_id, family, chimera flags
├── clone_tracking.tsv             # Clone family frequencies across rounds
├── ablang_scores.tsv              # Pseudo-perplexity scores (if --ablang)
├── diversity_stats.tsv            # Shannon, Simpson, dominance
├── convergence_analysis.tsv       # Round-to-round overlap
├── top_sequences_report.txt       # Human-readable report
├── top_enriched_sequences.fasta   # For AF2/RF2/ANARCI
├── enrichment_trajectories.png
├── top_sequence_detail.png
├── diversity_across_rounds.png
└── cdr3_length_distribution.png
```
