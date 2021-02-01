#include "TChain.h"
#include <iostream>
#include "TProof.h"

void starter(){

  TString inFileName = "bparknano.root";

  TString outFileName = "flat_bparknano.root";

  TChain * c = new TChain("Events");
  c->Add(inFileName);

  c->Process("NanoDumper.C", outFileName);

}
