#!/bin/bash
#SBATCH --job-name=rmats_tau_sorted
#SBATCH --partition=interruptible_cpu
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=rmats_sorted_log_%j.out
#SBATCH --error=rmats_sorted_log_%j.err

set -euo pipefail

source /scratch/prj/ppn_rnp_networks/users/sidra.bichay/software/Miniconda3/etc/profile.d/conda.sh
conda activate rmats_env

cd /scratch/prj/ppn_rnp_networks/users/sidra.bichay/hnRNPH1_IR_MAF

mkdir -p rmats_output_sorted rmats_tmp_sorted

echo "Starting rMATS analysis with plain sorted BAMs..."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Current directory: $(pwd)"
echo "Using $(which rmats.py)"
echo "Start time: $(date)"

rmats.py \
  --b1 b1_sorted.txt \
  --b2 b2_sorted.txt \
  --gtf gencode.vM25.annotation.gtf \
  --od rmats_output_sorted \
  --tmp rmats_tmp_sorted \
  --task both \
  -t paired \
  --libType fr-unstranded \
  --readLength 100 \
  --nthread "${SLURM_CPUS_PER_TASK}" \
  --tstat "${SLURM_CPUS_PER_TASK}" \
  --allow-clipping \
  --variable-read-length

echo "rMATS analysis complete!"
echo "End time: $(date)"
echo "Output files:"
ls -lh rmats_output_sorted/
