#!/bin/bash
#SBATCH --job-name=SashimiPlot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=_data/logs/sashimi_%j.log

# [[_hpc/run_sashimi.sbatch]]
# Usage: sbatch _hpc/run_sashimi.sbatch <cell_line> <comparison> <as_type> <event_id>

if [ "$#" -lt 4 ]; then
    echo "Usage: sbatch _hpc/run_sashimi.sbatch <cell_line> <comparison> <as_type> <event_id>"
    exit 1
fi

CELL_LINE=$1
COMP=$2
AS_TYPE=$3
EVENT_ID=$4

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate cancer_rnaseq

# Define script path
PIPELINE_ROOT="${SLURM_SUBMIT_DIR:-$(dirname "$(realpath "$0")")}"
SASHIMI_SCRIPT="$PIPELINE_ROOT/scripts_downstream/05_02_rmats_sashimi.sh"

echo "========================================================="
echo "   rMATS Sashimi Plot - Job Submission"
echo "   Starting at: $(date)"
echo "   Target: $CELL_LINE | $COMP | $AS_TYPE | ID: $EVENT_ID"
echo "========================================================="

bash "$SASHIMI_SCRIPT" "$CELL_LINE" "$COMP" "$AS_TYPE" "$EVENT_ID"

echo "========================================================="
echo "   Sashimi Plot completed at: $(date)"
echo "========================================================="
