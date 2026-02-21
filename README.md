# Nabbit

Nanobody NGS enrichment analysis pipeline for phage/yeast display selections from paired-end Illumina FASTQ data.

## What it does

1. Merges paired reads (R1+R2 overlap), translates to protein, quality-filters
2. Counts unique protein sequences per panning round, computes enrichment (fold-change, CPM, log2FC)
3. Annotates CDR1/2/3 and FR1-4 regions using IMGT-based VHH motif anchoring
4. Clusters sequences by CDR3 identity into clonotypes and families
5. Tracks clone families across rounds, flags chimeric sequences
6. Computes diversity metrics (Shannon, Simpson, dominance), convergence (Jaccard), and saturation (rarefaction + Chao1)
7. Optionally scores with AbLang (antibody language model pseudo-perplexity)
8. Optionally annotates V/D/J germline gene usage with IgBLAST
9. Generates static plots, an interactive HTML dashboard, and tabular/FASTA outputs

## Requirements

- Python 3.10+
- Operating system: macOS or Linux

### Core dependencies (required)

```
numpy
pandas
scipy
matplotlib
```

### Dashboard (required for interactive HTML report)

```
plotly
```

### Optional dependencies

| Package | Install | Enables |
|---------|---------|---------|
| `ablang` | `pip install ablang` | AbLang pseudo-perplexity scoring (`--ablang`) |
| `biopython` | `pip install biopython` | Phylogenetic trees, cluster PCA |
| `scikit-learn` | `pip install scikit-learn` | Cluster PCA projections |
| `igblast` | `conda install -c bioconda igblast` | V/D/J germline annotation (`--igblast`) |

### Quick install

```bash
# Create a virtual environment (recommended)
python3 -m venv .venv
source .venv/bin/activate

# Install all dependencies
pip install numpy pandas scipy matplotlib plotly

# Optional
pip install ablang biopython scikit-learn
```

### HPC (Biowulf / SLURM)

```bash
module load python/3.11
# numpy, pandas, scipy, matplotlib are typically pre-installed
pip install --user plotly ablang biopython scikit-learn
```

## Usage

### Web launcher (recommended for interactive use)

Start the built-in web UI for configuring and running the pipeline:

```bash
python nabbit.py --serve
```

This opens a browser at `http://localhost:8080` where you can:
- Browse for or drag-and-drop FASTQ files
- Assign round numbers
- Configure options (threads, top N, AbLang, IgBLAST)
- Run the pipeline and view results
- Open the interactive dashboard directly

To use a different port:

```bash
python nabbit.py --serve --port 9090
```

To stop the server, press `Ctrl+C` in the terminal, or kill the process:

```bash
# Find and kill the server
lsof -ti :8080 | xargs kill
```

### Command line

```bash
# Basic
python nabbit.py \
    --fastq-dir /path/to/fastqs \
    --output-dir /path/to/results \
    --threads 8

# With AbLang scoring and IgBLAST annotation
python nabbit.py \
    --fastq-dir /path/to/fastqs \
    --output-dir /path/to/results \
    --threads 8 \
    --ablang \
    --igblast

# Specific rounds only
python nabbit.py \
    --fastq-dir /path/to/fastqs \
    --output-dir /path/to/results \
    --rounds 1,2,3,5
```

### HPC batch job (SLURM)

```bash
# Basic
sbatch run_nabbit.sh /data/$USER/fastqs /data/$USER/results

# Specific rounds
sbatch run_nabbit.sh /data/$USER/fastqs /data/$USER/results "1,2,3,5,6"

# With AbLang
sbatch run_nabbit.sh /data/$USER/fastqs /data/$USER/results "" --ablang
```

## FASTQ naming

Expects paired-end files named `{round_number}_R{1,2}_001.fastq.gz`:

```
1_R1_001.fastq.gz    1_R2_001.fastq.gz    → Round 1
2_R1_001.fastq.gz    2_R2_001.fastq.gz    → Round 2
3_R1_001.fastq.gz    3_R2_001.fastq.gz    → Round 3
```

Alternatively, use a sample map TSV (`--sample-map map.tsv`):

```
# round	R1	R2
1	/path/to/round1_R1.fastq.gz	/path/to/round1_R2.fastq.gz
2	/path/to/round2_R1.fastq.gz	/path/to/round2_R2.fastq.gz
```

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--fastq-dir` | required | Directory containing FASTQ files |
| `--output-dir` | required | Output directory for results |
| `--sample-map` | none | TSV mapping file (round, R1 path, R2 path) |
| `--threads` | 1 | Number of parallel workers |
| `--min-overlap` | 20 | Minimum bp overlap for read merging |
| `--min-protein-len` | 90 | Minimum protein length in amino acids |
| `--max-protein-len` | 180 | Maximum protein length in amino acids |
| `--min-avg-qual` | 20.0 | Minimum average Phred quality score |
| `--top-n` | 50 | Number of top sequences to report |
| `--cdr3-cluster-threshold` | 0.8 | CDR3 identity threshold for family clustering |
| `--rounds` | all | Comma-separated round numbers to process |
| `--no-plots` | off | Skip static matplotlib plots |
| `--no-dashboard` | off | Skip interactive HTML dashboard |
| `--ablang` | off | Enable AbLang pseudo-perplexity scoring |
| `--ablang-device` | cpu | Device for AbLang (`cpu` or `cuda`) |
| `--igblast` | off | Enable IgBLAST V/D/J annotation |
| `--igblast-db` | `./igblast_refs` | Path to IgBLAST reference databases |
| `--igblast-chunk-size` | 0 (auto) | Batch size for IgBLAST queries |
| `--serve` | off | Start the web launcher UI |
| `--port` | 8080 | Port for the web server |

## Outputs

```
results/
├── processing_stats.tsv              # Read counts, merge rates per round
├── enrichment_table.tsv              # Top N sequences with counts/freq/fold-change
├── enrichment_full.tsv.gz            # All sequences (gzipped)
├── cdr_annotations.tsv               # FR1-4, CDR1-3, clone_id, family, chimera flags
├── clone_tracking.tsv                # Clone family frequencies across rounds
├── diversity_stats.tsv               # Shannon, Simpson, dominance metrics
├── convergence_analysis.tsv          # Round-to-round Jaccard overlap
├── top_sequences_report.txt          # Human-readable summary
├── top_enriched_sequences.fasta      # Top enriched sequences for downstream tools
├── top_diverse_sequences.fasta       # Top diverse sequences for downstream tools
├── ablang_scores.tsv                 # Pseudo-perplexity scores (if --ablang)
├── igblast_annotations.tsv           # V/D/J gene segments (if --igblast)
├── dashboard.html                    # Interactive Plotly dashboard
├── enrichment_trajectories.png       # Fold-change trajectories
├── top_sequence_detail.png           # Top sequence count/frequency
├── diversity_across_rounds.png       # Shannon diversity plots
└── cdr3_length_distribution.png      # CDR3 length histograms
```

## Interactive dashboard

The pipeline generates `dashboard.html` with 5 tabs:

- **Overview** — Read counts, productivity funnel, per-round phylogenetic trees
- **Enrichment** — Trajectories, volcano plot, top sequence detail, sortable table
- **Diversity** — Shannon/dominance, convergence heatmap, saturation curves, clone family dynamics
- **Clusters** — CDR3 length analysis, PCA, enriched cluster trees, variant heatmap
- **Scoring** — AbLang PPL vs enrichment, V/D/J gene usage (when `--ablang` or `--igblast` enabled)

Each chart includes an expandable figure legend explaining how to interpret the visualization.

## License

MIT
