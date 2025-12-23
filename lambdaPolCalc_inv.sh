#!/bin/bash
#SBATCH -D /lustre/home/user/v/vtroshin/mpd_hyperons/mpd_particle_finder/log
#SBATCH -J polLambdaCalc
#SBATCH -p mpd-ice --qos=dirac
#SBATCH -a 1-150
#SBATCH --requeue
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH -o /lustre/home/user/v/vtroshin/mpd_hyperons/mpd_particle_finder/log/%A_%a.log
#SBATCH --exclude=n02p002,n02p069,n02p064,n05p010
# Load necessary environment
# Wait for CVMFS and EOS to be available
ls /cvmfs/nica.jinr.ru/sw/os/login.sh
sleep 10

# loading MPDROOT framework
source /cvmfs/nica.jinr.ru/sw/os/login.sh
module add mpddev
export MPDROOT=/lustre/home/user/v/vtroshin/mpd
source $MPDROOT/config/env.sh

# Set job identifiers
export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

# Define paths (modify these as needed)

list_dir=/lustre/home/user/v/vtroshin/mpd_hyperons/mpd_uni/afterburner/out/list_dir/
input_list=/lustre/home/user/v/vtroshin/mpd_hyperons/mpd_uni/latest/urqmd_xexe_2.87gev_mf_unigen_enhanced_lamb_3m_ev_fix6.list

split -l 40 -d -a 4 --additional-suffix=.txt $input_list $list_dir
file_list=$( ls $list_dir | head -n ${TASK_ID} | tail -n 1 )
CONFIG_DIR="/lustre/home/user/v/vtroshin/mpd_hyperons/mpd_uni/latest/dar_flu_lar/global-polarization-/"
OUTPUT_DIR="/lustre/home/user/v/vtroshin/mpd_hyperons/mpd_uni/afterburner/out/analysis/for_merge/"


OUTPUT_FILE="result_urqmd_xexe_2.87gev_mf_${TASK_ID}_inv_0_comb_fix6.mcini.root"


# Run the analysis with timing measurement
echo "Starting analysis for task ${TASK_ID} at $(date)"
START_TIME=$(date +%s)

#calc_global_polarization("${list_dir}${file_list}","${OUTPUT_DIR}${OUTPUT_FILE}", 50);
#calc_global_polarization("${input_list}","${OUTPUT_DIR}${OUTPUT_FILE}", 50);
#calc_pol_vs_Nenh("${list_dir}${file_list}","${OUTPUT_DIR}${OUTPUT_FILE}",vecEnhancedValue);

# Run the analysis
root -l -b <<EOF
gSystem->Load("${CONFIG_DIR}AutoDict_vector_TVector3__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_UParticle__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_std__vector_ROOT__Math__XYZVector__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_TVector3__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_UParticle__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_std__vector_ROOT__Math__XYZVector__cxx.so")
std::vector<Int_t> vecEnhancedValue ={0, 1, 2, 5, 10, 15, 20};
.L ${CONFIG_DIR}read_unigen_root.cpp
calc_global_polarization("${list_dir}${file_list}","${OUTPUT_DIR}${OUTPUT_FILE}", 50);
.q
EOF

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "Analysis for task ${TASK_ID} completed at $(date)"
echo "Time taken: $((ELAPSED_TIME / 3600)) hours, $(( (ELAPSED_TIME % 3600) / 60 )) minutes, $((ELAPSED_TIME % 60)) seconds"

echo "Job ${TASK_ID} completed successfully"


