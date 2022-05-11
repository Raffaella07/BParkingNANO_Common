//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  1 11:00:38 2021 by ROOT version 6.12/07
// from TTree Runs/Runs
// found on file: bparknano_withgeninfo_nj79.root
//////////////////////////////////////////////////////////

#ifndef NanoRunDumper_h
#define NanoRunDumper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class NanoRunDumper : public TSelector {
  public :
    TTreeReader     fReader;  //!the tree reader
    TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

    // Readers to access the data (delete the ones you do not need).
    TTreeReaderValue<UInt_t> run = {fReader, "run"};
    TTreeReaderValue<Long64_t> genEventCount = {fReader, "genEventCount"};
    TTreeReaderValue<Double_t> genEventSumw = {fReader, "genEventSumw"};
    TTreeReaderValue<Double_t> genEventSumw2 = {fReader, "genEventSumw2"};
    TTreeReaderValue<UInt_t> nLHEScaleSumw = {fReader, "nLHEScaleSumw"};
    TTreeReaderArray<Double_t> LHEScaleSumw = {fReader, "LHEScaleSumw"};
    TTreeReaderValue<UInt_t> nLHEPdfSumw = {fReader, "nLHEPdfSumw"};
    TTreeReaderArray<Double_t> LHEPdfSumw = {fReader, "LHEPdfSumw"};


    NanoRunDumper(TTree * /*tree*/ =0) { }
    virtual ~NanoRunDumper() { }
    virtual Int_t   Version() const { return 2; }
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
    virtual void    SetOption(const char *option) { fOption = option; }
    virtual void    SetObject(TObject *obj) { fObject = obj; }
    virtual void    SetInputList(TList *input) { fInput = input; }
    virtual TList  *GetOutputList() const { return fOutput; }
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    // output file
    TFile* my_file;  

    // tree to fill
    TTree* run_tree;

    // filling variables
    Float_t the_geneventcount;
    Float_t the_geneventsumw;
    Float_t the_geneventsumw2;

    ClassDef(NanoRunDumper,0);

};

#endif

#ifdef NanoRunDumper_cxx
void NanoRunDumper::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);
}

Bool_t NanoRunDumper::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}


#endif // #ifdef NanoRunDumper_cxx
