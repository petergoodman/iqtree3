#!/bin/bash
#SBATCH --job-name=sim2-gtr20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=5gb
#SBATCH --cpus-per-task=84
# One shard of SIM-2's seeded GTR20 arm (decision 046; CODE_PLAN.md 3.8, "SIM-2's seeded GTR20 arm"). Submit from
# the run directory, outside the clone, with 84, 84 and 34 cores for shards 1, 2 and 3:
#   sbatch --account=ACCOUNT --partition=standard --cpus-per-task=CORES \
#       CLONE/test_scripts/fixedeq/sim/slurm/sim2_gtr20.sbatch SHARD RUNDIR FROZEN GATE CLONE
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
if [ "$#" -ne 5 ]; then
    echo "usage: sim2_gtr20.sbatch SHARD RUNDIR FROZEN GATE CLONE" >&2
    exit 2
fi
SHARD=$1 RUNDIR=$2 FROZEN=$3 GATE=$4 CLONE=$5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python "$CLONE/test_scripts/fixedeq/sim/sim2_gtr20.py" run --out "$RUNDIR" --shard "$SHARD" \
    --cores "$SLURM_CPUS_PER_TASK" --frozen "$FROZEN" --gate "$GATE" --git git
