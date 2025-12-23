# Global Polarization of Λ Hyperons at the MPD Experiment

The afterburner for hyperon global polarization setting and analysis. Utilizes UrQMD modeling data.

## Description

This project consists of two separate parts:

1. `simulate_lambda_decays(TString inputFile, TString outputFile, Int_t flag = 1, Int_t enhanceStat = 1)`  
   - Simulation of Λ decay into proton and pion
   - Global polarization is measured with proton polarization
   - Macro produces secondary particles and then puts polarization in proton

2. `calc_global_polarization(TString InFileName, TString OutFileName)`  
   - Measures global polarization with macro that provides transverse momentum and rapidity binning
   - Fits resulting angular distribution to get the value of global polarization

---

## Requirements

The afterburner is intended to work with UrQMD generated data in unigen format.

## Usage

First run root with commands:
```bash
gInterpreter->GenerateDictionary("vector<UParticle>", "vector;UParticle.h");
gInterpreter->GenerateDictionary("vector<TVector3>", "vector;TVector3.h");

```bash
# Example usage of functions
root -l -q 'simulate_lambda_decays.C("input.root", "output.root")'
root -l -q 'calc_global_polarization.C("output.root", "results.root")'
# Example usage of main macro
root -l -q 'root -x -q -b lambdaPolAnal.cpp'
root -l -q 'root -x -q -b lambdaPolCalc.cpp'
# mpd_afterburner_stable
