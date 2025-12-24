# Global Polarization of Λ Hyperons at the MPD Experiment

The afterburner for hyperon global polarization setting and analysis. Utilizes UrQMD modeling data.
## Description

Afterburner of lambda for HIC in UniGen format with realistic parameterization of directed and elliptic flow and global polarization(parapeterized from the exp data).

## Requirements

The afterburner is intended to work with UrQMD generated data in unigen format.

## Usage

First run generate_dicts.c in temporary directory
Second upload dictionaries into root session and run read_unigen_root.cpp
 
```bash
# Example usage of functions

cd ${CONFIG_DIR}
root -l -q ${CONFIG_DIR}/generate_dicts.c

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

# mpd_afterburner_stable
