

#ifdef __CLING__
//#pragma link C++ class LambdaAnalysisData+;
//#pragma link C++ class UUEvent+;
//#pragma link C++ class std::vector<UParticle>+;
//#pragma link C++ class std::vector<TVector3>+;
#endif

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TProfile.h>
#include <TLegend.h>
//#include <TRotationZ.h>
#include <TROOT.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TChain.h>
#include "UEvent.h"
#include "UParticle.h"
#include "URun.h"

#include "TInterpreter.h"


void simulate_lambda_decays(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat = 0);
void set_lambda_parameterization(TFile* Lambda_yield, Double_t fBVal, UParticle &ULambda); 

Double_t get_mean_polarization(Double_t sNN, Double_t centrality);// Value in % 
Double_t get_random_value(Double_t fMean, Double_t fSigma);
Double_t get_positive_phi(const Double_t& phi);
Double_t get_centrality  (Double_t fBVal);
Double_t get_costh(Double_t alpha, Double_t pol = 0.6);
Double_t get_V2   (Double_t sNN, Double_t centrality, Double_t lambda_pT, Double_t lambda_y);
ROOT::Math::XYZVector get_pol_lambda(UParticle& lambda, Double_t _fpoly, Double_t _fSigmaPol = 0.3);


class TString;
class TClonesArray;
class UParticle;


void simulate_lambda_decays(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat) {


    TChain *inChain = nullptr;
    TFile *outFile = nullptr;
    TTree *outTree = nullptr;
    
    // Important variables that were accidentally removed
    Int_t eventID;
    std::vector<ROOT::Math::XYZVector>* vecPolarization = nullptr;



    // Open input file as TChain
    inChain = new TChain("events");
    inChain->Add(inputFile);

    // Create output file and tree

    outFile = TFile::Open(outputFile, "UPDATE");
    if (!outFile || outFile->IsZombie()) {
        delete outFile;
        outFile = TFile::Open(outputFile, "RECREATE");
        std::cout << "Creating new output file: " << outputFile << std::endl;

    }

    TObjArray* fileElements = inChain->GetListOfFiles();
    if (fileElements && fileElements->GetEntries() > 0) {
        TChainElement* chEl = (TChainElement*)fileElements->At(0);
        TString firstFileName = chEl->GetTitle();

        TFile* fFirst = TFile::Open(firstFileName, "READ");
        if (fFirst && !fFirst->IsZombie()) {
            URun* inRun = nullptr;
            fFirst->GetObject("run", inRun);
            if (inRun) {
                outFile->cd();
                inRun->Write("run", TObject::kOverwrite);
                std::cout << "Copied URun from " << firstFileName << std::endl;
            } else {
                std::cerr << "WARNING: URun object not found in " 
                          << firstFileName << std::endl;
            }
            fFirst->Close();
        }
    }

    // Event variables
    UEvent  *inEvent = nullptr;
    UEvent  *outEvent = new UEvent();
    inChain->SetBranchAddress("event", &inEvent);

    outTree = (TTree*)outFile->Get("decays");
    bool newTree = false;
    if (!outTree) {

        outTree = new TTree("decays", "Lambda decay products");
        newTree = true;
        std::cout << "Creating new output tree" << std::endl;
        
        // Setup branches for new tree
        outTree->Branch("event", &outEvent);
        outTree->Branch("eventID", &eventID, "eventID/I");
        outTree->Branch("LambdaPolarization", &vecPolarization);


    } else {

        
        std::cout << "Appending to existing tree with " << outTree->GetEntries() << " entries" << std::endl;

        // Connect branches for existing tree
        outTree->SetBranchAddress("event", &outEvent);
        outTree->SetBranchAddress("eventID", &eventID);
        outTree->SetBranchAddress("LambdaPolarization", &vecPolarization);

    }


    TFile* Lambda_yield = TFile::Open(confInFile,"READ");


    // Particle instances
    Int_t enhancedFlag = 0;
    Int_t dummy[2];
    UParticle lambda(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.);

    // Physics parameters
    TRandom3* rand = new TRandom3(0);  
    const double mLambda = 1.115683;
    Double_t fSigmaPolVal = 1.5;
    const double anisotropy = 0.732;
    Double_t fBMin = 3.44;
    Double_t fBMax = 7.44;
    Int_t child_null[2] = {0, 0};


    // Process events
    Long64_t nEvents = inChain->GetEntries();

    std::cout << "Events: " << nEvents << std::endl;
    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {

        inChain->GetEntry(iEvent);
        outEvent->Clear(); // Clear previous event
        Int_t lambdaCounter = 0;

        // Copy event header information
        outEvent->SetEventNr(inEvent->GetEventNr());
        outEvent->SetB(inEvent->GetB());
        outEvent->SetPhi(inEvent->GetPhi());
        outEvent->SetNes(inEvent->GetNes());
        outEvent->SetStepNr(inEvent->GetStepNr());
        outEvent->SetStepT(inEvent->GetStepT());

        vecPolarization->clear();
        //first we COUNT lambdas
        for (Int_t i = 0; i < inEvent->GetNpa(); i++) {
            UParticle* part = inEvent->GetParticle(i);
            if (part->GetPdg() == 3122) lambdaCounter++;
        }

        //std::cout << "=====\nLambda counter in event #" << iEvent << " : " << lambdaCounter << " lambdas\n";

        for (Int_t i = 0; i < inEvent->GetNpa(); ++i) {

            UParticle* part = inEvent->GetParticle(i);
            if (part->GetPdg() != 3122) {
                outEvent->AddParticle(*part);
                ROOT::Math::XYZVector nullPol(0, 0, 0);
                vecPolarization->push_back(nullPol);
                continue; 
            }// Select Lambdas (PDG code 3122)
            else {

            
            lambda = *part;
            set_lambda_parameterization(Lambda_yield, inEvent->GetB(), lambda); 
            Double_t fPolY = get_mean_polarization(2.87, get_centrality(inEvent->GetB()));

            ROOT::Math::XYZVector pol = get_pol_lambda(lambda, fPolY/100., fSigmaPolVal);
            vecPolarization->push_back(pol);

            outEvent->AddParticle(lambda);

            //std::cout << "Lambda counter in event #" << iEvent << " b4 enhance : " << ULambda->size() << " lambdas\n"; 



            Double_t fEnhanceStat = enhanceStat;
            while(fEnhanceStat > 1){ //enhancing of lambdas
                enhancedFlag = -9;

                TLorentzVector mom_rand( 1., 1., 1., 1. );
                
                TLorentzVector pos_rand(
                    get_random_value(lambda.X(), 0.03),
                    get_random_value(lambda.Y(), 0.03),
                    get_random_value(lambda.Z(), 0.03),
                    get_random_value(lambda.T(), 0.03)
                );

                UParticle enhancedLambda(lambda);
                enhancedLambda.SetPosition(pos_rand);
                enhancedLambda.SetMate(enhancedFlag);

                set_lambda_parameterization(Lambda_yield, inEvent->GetB(), enhancedLambda); 
                fPolY = get_mean_polarization(2.87, get_centrality(inEvent->GetB()));
                

                ROOT::Math::XYZVector pol = get_pol_lambda(enhancedLambda, fPolY/100., fSigmaPolVal);
                vecPolarization->push_back(pol);

                outEvent->AddParticle(enhancedLambda);

    
                //std::cout << "ooooo Lambda counter in event #" << iEvent << " inside enhance : " << ULambda->size() << " lambdas\n"; 
                fEnhanceStat--;
            }
	}
        }
        //std::cout << "Lambda counter in event #" << iEvent << " after add & enhance : " << ULambda->size() << " lambdas\n"; 

        outTree->Fill();
    }

    outFile->cd();
    outTree->Write("",TObject::kOverwrite);


    delete outEvent;
    delete rand;
    delete inChain;

    outFile->Close();
    Lambda_yield->Close();
}




Double_t get_positive_phi(const Double_t& phi){
    Double_t phi_rot = phi;
    if (phi < 0) phi_rot += 2.*TMath::Pi();
    return phi_rot;
}


Double_t get_costh(Double_t alpha, Double_t pol) {

    TF1* randomCosTheta = new TF1("randomCosTheta", "(1+[0]*x)", -1.0, 1.0);
    randomCosTheta->SetNpx(10000);

    // Set the parameters
    randomCosTheta->SetParameters(alpha);

    // Generate a random value from the distribution
    Double_t value = randomCosTheta->GetRandom();
    // Clean up
    delete randomCosTheta;
    
    return value;

}

Double_t get_random_value(Double_t fMean, Double_t fSigma)
{
    TRandom3* rand = new TRandom3(0);  // Seed with 0 for reproducibility
    Double_t randomVal = rand->Gaus(fMean, fSigma);
    delete rand;
    return randomVal;
}



    void set_lambda_parameterization(TFile* Lambda_yield, Double_t fBVal, UParticle &ULambda){ //Valeriy's function for properly lambda generation

    TRandom3* rand = new TRandom3(0);  // Seed with 0 for reproducibility

	Double_t centrality = get_centrality(fBVal); //centrality for parameterization (b->centrality for Xe+Xe below)
	Double_t sNN = 2.87; // Energy of the collision in center-of-mass system

	TH2F* h_pt_y;
	if(centrality<10) h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT010");
	else if(centrality>10 && centrality <40) h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT");
	else h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT40100");


    Double_t lambda_pT; // pT of Lambda from pT-y TH2F
	Double_t lambda_y;  // rapidity of Lambda from pT-y TH2F

	h_pt_y->GetRandom2(lambda_y,lambda_pT,rand);

	Double_t v1 = (28.8635/(TMath::Power(sNN,2.89092))) 
                * ((-0.0233*centrality+0.5413* TMath::Power(centrality,1./3) ) 
                * (163.536/18.0188 - 163.536/(lambda_pT +18.0188) )* lambda_y + ( -0.0056*centrality+0.377*TMath::Power(centrality,1./3) ) 
                * (0.6653* lambda_pT - 0.6172 * TMath::Power(lambda_pT,2) +0.1154 * TMath::Power(lambda_pT,3) ) * TMath::Power(lambda_y,3) ); // v1 cent-pT-y func
    Double_t v2 = get_V2(sNN, centrality, lambda_pT, lambda_y);
    if( v1  > 1 )  v1 =1;
    else if(v1<-1) v1=-1;
     if( v2  > 1 ) v2 =1;
     if( v2 <-1) v2=-1;

    static TF1 f("f", "[0]*(1+2*[1]*TMath::Cos(x)+2*[2]*TMath::Cos(2*x))+[3]", 0,2*TMath::Pi());
    Double_t a1=1+2*v1+2*v2;
    Double_t a2=1-2*v1+2*v2;
    Double_t a3=1-v1*v1/(4*v2)-2*v2;
    Double_t a=0; //
    if (a1<a) a=a1;  // find analytic minimun to shift
    if (a2<a) a=a2;  // find analytic minimun to shift
    if (a3<a) a=a3;  // find analytic minimun to shift
    f.SetParameter(0,1/(2*TMath::Pi()*(1-a))); // norm
    f.SetParameter(1,v1);  // v1
    f.SetParameter(2,v2);  // v2
    f.SetParameter(3,-a/(2*TMath::Pi()*(1-a))); // shift to have probability
    f.SetNpx(10000);  // to get a better result when using TF1::GetRandom
    Double_t phi=f.GetRandom(rand);

    // Calculate total energy (E) properly
    TLorentzVector vec;
    // Calculate transverse mass (m_T)
    Double_t lambda_mass = 1.115683; // GeV/c² (Lambda mass)
    Double_t mT = sqrt(lambda_pT * lambda_pT + lambda_mass * lambda_mass);
    // Convert rapidity (y) to pseudorapidity (η)
    Double_t pz = mT * sinh(lambda_y); // longitudinal momentum
    Double_t lambda_eta = (pz != 0.0) ? TMath::ATanH(pz / sqrt(pz*pz + lambda_pT*lambda_pT)) : 0.0;
    Double_t fEnergyLambda = sqrt(lambda_pT*lambda_pT * cosh(lambda_eta)*cosh(lambda_eta) + lambda_mass*lambda_mass);
    vec.SetPtEtaPhiE(lambda_pT, lambda_eta, phi, fEnergyLambda);
    ULambda.SetMomentum(vec); //ULambda

    delete rand;


}


Double_t get_V2(Double_t sNN, Double_t centrality, Double_t lambda_pT, Double_t lambda_y) {

    lambda_y = TMath::Abs(lambda_y);
    
    // =============================================
    // 1. Energy-dependent term (depends on sNN)
    // =============================================
    Double_t energy_term = 1.05 * (0.8132 - 11.4 / TMath::Power(sNN, 2.1));

    // =============================================
    // 2. Centrality-dependent term (normalized at centrality=25)
    // =============================================
    Double_t centrality_numerator = 
        -0.01105 * centrality + 0.000162 * TMath::Power(centrality, 2);
    Double_t centrality_denominator = 
        -0.01105 * 25 + 0.000162 * TMath::Power(25, 2);
    Double_t centrality_term = centrality_numerator / centrality_denominator;

    // =============================================
    // 3. Sign flip (-1 factor)
    // =============================================
    Double_t sign_flip = -1;

    // =============================================
    // 4. Rapidity (lambda_y) and pT (lambda_pT) dependent term
    // =============================================
    // 4a. Linear term in lambda_pT (depends on lambda_y)
    Double_t linear_pT_term = 
        (0.32172 * TMath::Power(lambda_y, 1.805) - 0.1578) * lambda_pT;

    // 4b. Nonlinear term (exponential + sinusoidal dependence)
    Double_t exp_argument = 
        lambda_pT * (1.4562 * TMath::Power(lambda_y, 3.28814) + 0.22912);
    Double_t sin_phase = 
        lambda_pT * (-1.004 * TMath::Power(lambda_y, 4.6398) + 3.88097) 
        + (-3.6075 * TMath::Power(lambda_y, 4.6582) + 3.35966);
    
    Double_t nonlinear_term = 
        (0.05838 * TMath::Power(lambda_y, 3.2172) - 0.0179) 
        * TMath::Exp(exp_argument) 
        * TMath::Sin(sin_phase);

    // Combine all terms into final v2 expression
    Double_t v2 = 
        energy_term 
        * centrality_term 
        * sign_flip 
        * (linear_pT_term + nonlinear_term);


    return v2;

}

Double_t get_mean_polarization(Double_t sNN, Double_t centrality){ return (2.8569/ TMath::Power(sNN,0.955513) ) * (2.4702 - 0.0461*centrality + 0.0042 * TMath::Power(centrality, 2)); }// Value in % 

Double_t get_centrality(Double_t fBVal){
    if(fBVal<3.44) {return 5;}
    else if(fBVal<4.88) {return 15;}
    else if(fBVal<5.84) {return 25;}
    else if(fBVal<6.64) {return 35;}
    else if(fBVal<7.44) {return 45;}
    else if(fBVal<8.08) {return 55;}
    else if(fBVal<8.72) {return 65;}
    else if(fBVal<9.36) {return 75;}
    else if(fBVal<9.84) {return 85;}
    else {return 95;}
}


ROOT::Math::XYZVector get_pol_lambda(UParticle& lambda, Double_t _fpoly, Double_t _fSigmaPol){ 
    
    Double_t fPhi = 0;

    Double_t fpolx  = 0; 
    Double_t fpolSx = 0.07;  
    Double_t fpoly  = _fpoly;
    Double_t fpolSy = 0.07; 
    Double_t fpolz  = 0;
    Double_t fpolSz = 0.07;

    TRandom3* rand = new TRandom3(0);  // Seed with 0 for reproducibility
    Double_t polX = rand->Gaus(fpolx,fpolSx);
    Double_t polY = rand->Gaus(fpoly,fpolSy);
    Double_t polZ = rand->Gaus(fpolz,fpolSz);

    ROOT::Math::XYZVector polarizationVec = ROOT::Math::XYZVector(polX, polY, polZ);

    delete rand;
    return polarizationVec;
}

