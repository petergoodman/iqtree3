#!/bin/bash
#SBATCH --job-name=sim2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=5gb
#SBATCH --cpus-per-task=94
# One shard of SIM-2 (CODE_PLAN.md 3.8, "SIM-2's plan"). Submit from the run directory, outside the clone:
#   sbatch --account=ACCOUNT --partition=standard --cpus-per-task=94 CLONE/test_scripts/fixedeq/sim/slurm/sim2.sbatch \
#       SHARD SHARDS RUNDIR FROZEN GATE CLONE
# No account or path is written in this file. Slurm runs a copy of this script, so the clone is an argument.

# conda loads through a module, before set -euo pipefail (CLAUDE.md, Puma)
module load anaconda
source ~/.bashrc && conda activate
conda activate iqtree3-fixedeq

set -euo pipefail
if [ "${CONDA_DEFAULT_ENV:-}" != "iqtree3-fixedeq" ]; then
    echo "the iqtree3-fixedeq environment is not active" >&2
    exit 1
fi
if [ "$#" -ne 6 ]; then
    echo "usage: sim2.sbatch SHARD SHARDS RUNDIR FROZEN GATE CLONE" >&2
    exit 2
fi
SHARD=$1 SHARDS=$2 RUNDIR=$3 FROZEN=$4 GATE=$5 CLONE=$6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python "$CLONE/test_scripts/fixedeq/sim/sim2.py" run --out "$RUNDIR" --shard "$SHARD" --shards "$SHARDS" \
    --cores "$SLURM_CPUS_PER_TASK" --frozen "$FROZEN" --gate "$GATE" --git git
