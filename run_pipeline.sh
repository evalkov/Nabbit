#!/bin/bash
#SBATCH --job-name=nanobody_ngs
#SBATCH --output=nanobody_ngs_%j.out
#SBATCH --error=nanobody_ngs_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=4:00:00
#SBATCH --partition=norm

# Usage:
#   sbatch run_pipeline.sh /data/$USER/fastqs /data/$USER/results
#   sbatch run_pipeline.sh /data/$USER/fastqs /data/$USER/results "1,2,3,5,6"
#   sbatch run_pipeline.sh /data/$USER/fastqs /data/$USER/results "" --ablang

FASTQ_DIR="${1:?Error: provide FASTQ directory as arg 1}"
OUTPUT_DIR="${2:?Error: provide output directory as arg 2}"
ROUNDS="${3:-}"
EXTRA="${4:-}"

module load python/3.11

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CMD="python ${SCRIPT_DIR}/nanobody_ngs_enrichment.py \
    --fastq-dir $FASTQ_DIR \
    --output-dir $OUTPUT_DIR \
    --threads $SLURM_CPUS_PER_TASK"

[ -n "$ROUNDS" ] && CMD="$CMD --rounds $ROUNDS"
[ -n "$EXTRA" ] && CMD="$CMD $EXTRA"

echo "$(date) | Node: $(hostname) | CPUs: $SLURM_CPUS_PER_TASK"
echo "CMD: $CMD"
echo ""
eval $CMD
echo ""
echo "$(date) | Done"
