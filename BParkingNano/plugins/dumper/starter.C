// example of usage: root -l starter.C+\(\"true\"\)

#include "TChain.h"
#include <iostream>
#include "TProof.h"

void starter(TString isMC){
  TString inFileName = "bparknano.root";

  TString outFileName = "flat_bparknano.root";

  TChain* c = new TChain("Events");
  c->Add(inFileName);

  c->Process("NanoDumper.C+", outFileName);

  if(isMC == "true"){
    TChain* c_run = new TChain("Runs");
    c_run->Add(inFileName);

    c_run->Process("NanoRunDumper.C+", outFileName);
  }
}
