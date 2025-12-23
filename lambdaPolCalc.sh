#!/bin/bash
#SBATCH -D /lustre/home/user/v/vtroshin/mpd_hyperons/mpd_particle_finder/log
#SBATCH -J polLambdaCalc
#SBATCH -p mpd-ice --qos=dirac
#SBATCH -a 1-1
#SBATCH --requeue
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH -o /lustre/home/user/v/vtroshin/mpd_hyperons/mpd_particle_finder/log/%A_%a.log
#SBATCH --exclude=n05p020,n05p014,n05p004,n05p006,n05p009,n05p025
# Load necessary environment
# Wait for CVMFS and EOS to be available
ls /cvmfs/nica.jinr.ru/sw/os/login.sh
sleep 10

# loading MPDROOT framework
source /cvmfs/nica.jinr.ru/sw/os/login.sh
module add mpddev/v25.09.25-1
export MPDROOT=/lustre/home/user/v/vtroshin/mpd
source $MPDROOT/config/env.sh

# Set job identifiers
export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

# Define paths (modify these as needed)

list_dir=/lustre/home/user/v/vtroshin/mpd_hyperons/mpd_uni/afterburner/out/list_dir/
input_list=/lustre/home/user/v/vtroshin/mpd_hyperons/mpd_uni/latest/urqmd_xexe_2.87gev_mf_unigen_50enhanced_lamb_1m_ev_011225.list

split -l 4000 -d -a 4 --additional-suffix=.txt $input_list $list_dir
file_list=$( ls $list_dir | head -n ${TASK_ID} | tail -n 1 )
CONFIG_DIR="/lustre/home/user/v/vtroshin/mpd_hyperons/mpd_uni/latest/dar_flu_lar/global-polarization-/"
OUTPUT_DIR="/lustre/home/user/v/vtroshin/mpd_hyperons/mpd_uni/afterburner/out/analysis/"


OUTPUT_FILE="result_urqmd_xexe_2.87gev_mf_${TASK_ID}_2m_011225.mcini.root"


# Run the analysis with timing measurement
echo "Starting analysis for task ${TASK_ID} at $(date)"
START_TIME=$(date +%s)

#calc_global_polarization("${list_dir}${file_list}","${OUTPUT_DIR}${OUTPUT_FILE}", 50);
# Run the analysis
root -l -b <<EOF
gSystem->Load("${CONFIG_DIR}AutoDict_vector_TVector3__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_UParticle__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_TVector3__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_UParticle__cxx.so")
.L ${CONFIG_DIR}read_unigen_root.cpp
calc_global_polarization("${input_list}","${OUTPUT_DIR}${OUTPUT_FILE}", 50);
.q
EOF

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "Analysis for task ${TASK_ID} completed at $(date)"
echo "Time taken: $((ELAPSED_TIME / 3600)) hours, $(( (ELAPSED_TIME % 3600) / 60 )) minutes, $((ELAPSED_TIME % 60)) seconds"

echo "Job ${TASK_ID} completed successfully"


