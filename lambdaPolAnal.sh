#!/bin/bash
#SBATCH -D /scratch2/troshin/log
#SBATCH -J polLambdaAnal
#SBATCH -p lustre-test
#SBATCH -a 1-1
#SBATCH --requeue
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH -o /scratch2/troshin/log/%A_%a.log
# Load necessary environment
# Wait for CVMFS and EOS to be available
ls /cvmfs/nica.jinr.ru/sw/os/login.sh
sleep 10

# loading MPDROOT framework
source /cvmfs/nica.jinr.ru/sw/os/login.sh
module add mpddev/v25.09.25-1
export MPDROOT=/scratch2/troshin/mpd_hyperons/mpd
source $MPDROOT/config/env.sh

# Set job identifiers
export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

# Define paths (modify these as needed)
INPUT_DIR="/eos/nica/mpd/users/parfenov/UrQMD/xexe_2.87gev_mf/12948799/files/mcini/"
CONFIG_DIR="/scratch2/troshin/govorun_backup/mpd_uni/latest/dar_flu_lar/global-polarization-/"
OUTPUT_DIR="/scratch2/troshin/mpd_hyperons/mpd_xexe_2.87_enh_lamb_20/afterburner_out/final12948799/"

mkdir -p $OUTPUT_DIR
# Create unique working directory
WORK_DIR="${CONFIG_DIR}lambda_analysis"
mkdir -p "${WORK_DIR}" || { echo "Failed to create working directory"; exit 1; }
cd "${WORK_DIR}" || { echo "Failed to enter working directory"; exit 1; }
echo "Enter working directory"
# Build file names
INPUT_FILE="urqmd_xexe_2.87gev_mf_12948799_${TASK_ID}.mcini.root"
CONFIG_FILE="lamb_config.root"
OUTPUT_FILE="result_urqmd_xexe_2.87gev_mf_12948799_${TASK_ID}.mcini.root"


# Run the analysis with timing measurement
echo "Starting analysis for task ${TASK_ID} at $(date)"
START_TIME=$(date +%s)

cd ${CONFIG_DIR}
root -l -q ${CONFIG_DIR}/generate_dicts.c
# Run the analysis
root -l -b <<EOF
gSystem->Load("${CONFIG_DIR}AutoDict_vector_TVector3__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_UParticle__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_std__vector_ROOT__Math__XYZVector__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_TVector3__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_UParticle__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_std__vector_ROOT__Math__XYZVector__cxx.so")
.L ${CONFIG_DIR}read_unigen_root.cpp
read_unigen_root("${INPUT_DIR}${INPUT_FILE}", "${OUTPUT_DIR}${OUTPUT_FILE}", "${CONFIG_DIR}${CONFIG_FILE}", 0)
.q
EOF

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "Analysis for task ${TASK_ID} completed at $(date)"
echo "Time taken: $((ELAPSED_TIME / 3600)) hours, $(( (ELAPSED_TIME % 3600) / 60 )) minutes, $((ELAPSED_TIME % 60)) seconds"

# Verify and copy results
if [[ -f "${OUTPUT_FILE}" ]]; then
    cp "${OUTPUT_FILE}" "${OUTPUT_DIR}"/
    # cp "${POLARIZATION_OUTPUT}" "${OUTPUT_DIR}"/
else
    echo "Error: Output file was not created!"
    exit 1
fi

echo "Job ${TASK_ID} completed successfully"


