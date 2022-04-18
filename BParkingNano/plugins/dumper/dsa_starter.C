#include "TChain.h"
#include <iostream>
#include "TProof.h"

void dsa_starter(){
  Bool_t isMC = true;

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_selectedslimmed_looseselectiondsa_loosedeltaptrel.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/big_bparknano_selectedslimmed_looseselectiondsa_loosedeltaptrel.root";
  TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/bparknano_selectedslimmed_looseselectiondsa_loosedeltaptrel.root";

  TString outFileName = "flat_bparknano.root";
  if(isMC) outFileName += "_isMC";

  TChain * c = new TChain("Events");
  c->Add(inFileName);

  c->Process("DSAMuonAnalyser.C+", outFileName);
}
