#define HNLToMuPiDumper_cxx
// The class definition in HNLToMuPiDumper.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("HNLToMuPiDumper.C")
// root> T->Process("HNLToMuPiDumper.C","some options")
// root> T->Process("HNLToMuPiDumper.C+")
//


#include "HNLToMuPiDumper.h"
#include <TMath.h>
#include <cmath>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
//#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "utils.C"


using namespace std;


void HNLToMuPiDumper::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  cout << " --------------------------" << endl;
  cout << "       HNL->mupi Dumper    " << endl;
  cout << " --------------------------" << endl;
}


void HNLToMuPiDumper::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  TString outFileName = option;

  if(outFileName.Contains("isMC")){
    isMC = true;
    outFileName.Resize(outFileName.Length()-5);
  }
  else isMC = false;

  // check if outputfile exists
  if(gSystem->AccessPathName(outFileName)){
    my_file = new TFile(outFileName, "RECREATE");  
  }
  else{
    my_file = new TFile(outFileName, "UPDATE");  
  }
  my_file->cd();

  // if MC, get the correct content of the GenPart branches
  if(isMC){
    nGenPart = {fReader, "nGenPart"};
    GenPart_eta = {fReader, "GenPart_eta"};
    GenPart_mass = {fReader, "GenPart_mass"};
    GenPart_phi = {fReader, "GenPart_phi"};
    GenPart_pt = {fReader, "GenPart_pt"};
    GenPart_vx = {fReader, "GenPart_vx"};
    GenPart_vy = {fReader, "GenPart_vy"};
    GenPart_vz = {fReader, "GenPart_vz"};
    GenPart_genPartIdxMother = {fReader, "GenPart_genPartIdxMother"};
    GenPart_pdgId = {fReader, "GenPart_pdgId"};
    GenPart_status = {fReader, "GenPart_status"};
    GenPart_statusFlags = {fReader, "GenPart_statusFlags"};
    Muon_genPartIdx = {fReader, "Muon_genPartIdx"};
  }

  // getting the tree ready
  hnl_tree = new TTree("hnl_tree", "hnl_tree");

  hnl_tree->Branch("event", &the_event);
  hnl_tree->Branch("run", &the_run);
  hnl_tree->Branch("lumi", &the_lumi);

  hnl_tree->Branch("hnl_pt", &the_hnl_pt);
  hnl_tree->Branch("hnl_eta", &the_hnl_eta);
  hnl_tree->Branch("hnl_phi", &the_hnl_phi);
  hnl_tree->Branch("hnl_mass", &the_hnl_mass);
  hnl_tree->Branch("hnl_mass_unfitted", &the_hnl_mass_unfitted);
  hnl_tree->Branch("diff_mass", &the_diff_mass);
  hnl_tree->Branch("hnl_charge", &the_hnl_charge);
  hnl_tree->Branch("hnl_cos2d", &the_hnl_cos2d);

  hnl_tree->Branch("mu_pt", &the_mu_pt);
  hnl_tree->Branch("mu_eta", &the_mu_eta);
  hnl_tree->Branch("mu_phi", &the_mu_phi);
  hnl_tree->Branch("mu_charge", &the_mu_charge);
  hnl_tree->Branch("mu_dxy", &the_mu_dxy);
  hnl_tree->Branch("mu_dxysig", &the_mu_dxysig);
  hnl_tree->Branch("mu_dz", &the_mu_dz);
  hnl_tree->Branch("mu_dzsig", &the_mu_dzsig);
  hnl_tree->Branch("mu_looseid", &the_mu_looseid);
  hnl_tree->Branch("mu_mediumid", &the_mu_mediumid);
  hnl_tree->Branch("mu_tightid", &the_mu_tightid);
  hnl_tree->Branch("mu_softid", &the_mu_softid);
  hnl_tree->Branch("mu_istriggering", &the_mu_istriggering);
  hnl_tree->Branch("mu_isslimmed", &the_mu_isslimmed);
  hnl_tree->Branch("mu_isdsa", &the_mu_isdsa);

  hnl_tree->Branch("pi_pt", &the_pi_pt);
  hnl_tree->Branch("pi_eta", &the_pi_eta);
  hnl_tree->Branch("pi_phi", &the_pi_phi);
  hnl_tree->Branch("pi_charge", &the_pi_charge);
  hnl_tree->Branch("pi_dcasig", &the_pi_dcasig);
  hnl_tree->Branch("pi_dxy", &the_pi_dxy);
  hnl_tree->Branch("pi_dz", &the_pi_dz);
  hnl_tree->Branch("pi_dxysig", &the_pi_dxysig);
  hnl_tree->Branch("pi_dzsig", &the_pi_dzsig);

  hnl_tree->Branch("pv_npvs", &the_pv_npvs);

  hnl_tree->Branch("sv_chi2", &the_sv_chi2);
  hnl_tree->Branch("sv_lxy", &the_sv_lxy);
  hnl_tree->Branch("sv_lxysig", &the_sv_lxysig);
  hnl_tree->Branch("sv_prob", &the_sv_prob);

  hnl_tree->Branch("ismatched", &the_ismatched);
  hnl_tree->Branch("mu_ismatched", &the_mu_ismatched);
  hnl_tree->Branch("pi_ismatched", &the_pi_ismatched);
  hnl_tree->Branch("mupi_mass_reco_gen_reldiff", &the_mupi_mass_reco_gen_reldiff);
  hnl_tree->Branch("lxy_reco_gen_reldiff", &the_lxy_reco_gen_reldiff);

  if(isMC){
    hnl_tree->Branch("gen_b_pt", &the_gen_b_pt);
    hnl_tree->Branch("gen_b_eta", &the_gen_b_eta);
    hnl_tree->Branch("gen_b_phi", &the_gen_b_phi);
    hnl_tree->Branch("gen_b_mass", &the_gen_b_mass);
    hnl_tree->Branch("gen_b_pdgid", &the_gen_b_pdgid);
    hnl_tree->Branch("gen_hnl_ct", &the_gen_hnl_ct);
    hnl_tree->Branch("gen_hnl_lxy", &the_gen_hnl_lxy);
    hnl_tree->Branch("gen_hnl_lxyz", &the_gen_hnl_lxyz);
    hnl_tree->Branch("gen_hnl_pt", &the_gen_hnl_pt);
    hnl_tree->Branch("gen_hnl_eta", &the_gen_hnl_eta);
    hnl_tree->Branch("gen_hnl_phi", &the_gen_hnl_phi);
    hnl_tree->Branch("gen_hnl_mass", &the_gen_hnl_mass);
    hnl_tree->Branch("gen_hnl_vx", &the_gen_hnl_vx);
    hnl_tree->Branch("gen_hnl_vy", &the_gen_hnl_vy);
    hnl_tree->Branch("gen_hnl_vz", &the_gen_hnl_vz);
    hnl_tree->Branch("gen_trgmu_pt", &the_gen_trgmu_pt);
    hnl_tree->Branch("gen_trgmu_eta", &the_gen_trgmu_eta);
    hnl_tree->Branch("gen_trgmu_phi", &the_gen_trgmu_phi);
    hnl_tree->Branch("gen_trgmu_vx", &the_gen_trgmu_vx);
    hnl_tree->Branch("gen_trgmu_vy", &the_gen_trgmu_vy);
    hnl_tree->Branch("gen_trgmu_vz", &the_gen_trgmu_vz);
    hnl_tree->Branch("gen_mu_pt", &the_gen_mu_pt);
    hnl_tree->Branch("gen_mu_eta", &the_gen_mu_eta);
    hnl_tree->Branch("gen_mu_phi", &the_gen_mu_phi);
    hnl_tree->Branch("gen_mu_vx", &the_gen_mu_vx);
    hnl_tree->Branch("gen_mu_vy", &the_gen_mu_vy);
    hnl_tree->Branch("gen_mu_vz", &the_gen_mu_vz);
    hnl_tree->Branch("gen_pi_pt", &the_gen_pi_pt);
    hnl_tree->Branch("gen_pi_eta", &the_gen_pi_eta);
    hnl_tree->Branch("gen_pi_phi", &the_gen_pi_phi);
    hnl_tree->Branch("gen_pi_vx", &the_gen_pi_vx);
    hnl_tree->Branch("gen_pi_vy", &the_gen_pi_vy);
    hnl_tree->Branch("gen_pi_vz", &the_gen_pi_vz);
  }
}


Bool_t HNLToMuPiDumper::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetLocalEntry(entry);
  //cout << endl << "--- Entry " << entry << " ---" << endl;

  // for data, we skip the event in case it doesn't pass the lumi mask
  if(!isMC && lumiMask(*run, *luminosityBlock) == false) return false;

  // number of candidates in the event
  UInt_t nCand_sig = *nHNLToMuPi; 

  // branches common to signal and control channels
  the_event = *event; 
  the_run = *run;
  the_lumi = *luminosityBlock; 

  the_pv_npvs = *PV_npvs;

  if(nCand_sig > 0){ // at least one candidate per event

    // selecting the candidate as the one having the largest hnl pt
    // - create candIdx - cos2d pairs
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sig = createPairWithDesc(nCand_sig, HNLToMuPi_hnl_cos2D);
    // - sort it in decreasing cos2d
    stable_sort(pair_candIdx_desc_cos2d_sig.begin(), pair_candIdx_desc_cos2d_sig.end(), sortcansbydesc);

    // - then privilege OS cand over SS ones
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sig, HNLToMuPi_hnl_charge);
    stable_sort(pair_candIdx_desc_cos2d_sign_sig.begin(), pair_candIdx_desc_cos2d_sign_sig.end(), sortcansbydesc_opp);

    // - for signal, priviledge matched candidates
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_matched_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sign_sig, HNLToMuPi_isMatched);
    stable_sort(pair_candIdx_desc_cos2d_sign_matched_sig.begin(), pair_candIdx_desc_cos2d_sign_matched_sig.end(), sortcansbydesc);

    // - finally priviledge candidate with the smallest deltaMassRelDiff(reco, gen)
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_matched_mass_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sign_matched_sig, HNLToMuPi_mupi_mass_reco_gen_reldiff);
    stable_sort(pair_candIdx_desc_cos2d_sign_matched_sig.begin(), pair_candIdx_desc_cos2d_sign_matched_sig.end(), sortcansbydesc_opp);

    // - and select the OS cand with the largest cos2d
    UInt_t selectedCandIdx_sig = pair_candIdx_desc_cos2d_sign_matched_mass_sig[0].first;

    // fill the hnl_tree
    the_hnl_pt = HNLToMuPi_hnl_pt[selectedCandIdx_sig];
    the_hnl_eta = HNLToMuPi_hnl_eta[selectedCandIdx_sig];
    the_hnl_phi = HNLToMuPi_hnl_phi[selectedCandIdx_sig];
    the_hnl_mass = HNLToMuPi_hnl_mass[selectedCandIdx_sig];
    the_hnl_mass_unfitted = HNLToMuPi_hnl_mass_unfitted[selectedCandIdx_sig];
    the_diff_mass = HNLToMuPi_diff_mass[selectedCandIdx_sig];
    the_hnl_charge = HNLToMuPi_hnl_charge[selectedCandIdx_sig];
    the_hnl_cos2d = HNLToMuPi_hnl_cos2D[selectedCandIdx_sig];

    the_mu_pt = HNLToMuPi_fit_mu_pt[selectedCandIdx_sig];
    the_mu_eta = HNLToMuPi_fit_mu_eta[selectedCandIdx_sig];
    the_mu_phi = HNLToMuPi_fit_mu_phi[selectedCandIdx_sig]; 
    the_mu_charge = Muon_charge[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];
    the_mu_dxy = Muon_dxy[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];
    the_mu_dxysig = Muon_dxyS[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];
    the_mu_dzsig = Muon_dzS[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];
    the_mu_looseid = Muon_looseId[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];
    the_mu_mediumid = Muon_mediumId[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];
    the_mu_tightid = Muon_tightId[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];
    the_mu_softid = Muon_softId[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];
    the_mu_istriggering = Muon_isTriggering[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];
    the_mu_isslimmed = Muon_isSlimmedMuon[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];
    the_mu_isdsa = Muon_isDSAMuon[HNLToMuPi_sel_mu_idx[selectedCandIdx_sig]];

    the_pi_pt = HNLToMuPi_fit_pi_pt[selectedCandIdx_sig];
    the_pi_eta = HNLToMuPi_fit_pi_eta[selectedCandIdx_sig];
    the_pi_phi = HNLToMuPi_fit_pi_phi[selectedCandIdx_sig]; 
    the_pi_charge = ProbeTracks_charge[HNLToMuPi_pi_idx[selectedCandIdx_sig]];
    the_pi_dcasig = ProbeTracks_DCASig[HNLToMuPi_pi_idx[selectedCandIdx_sig]];
    the_pi_dxy = ProbeTracks_dxy[HNLToMuPi_pi_idx[selectedCandIdx_sig]];
    the_pi_dz = ProbeTracks_dz[HNLToMuPi_pi_idx[selectedCandIdx_sig]];
    the_pi_dxysig = ProbeTracks_dxyS[HNLToMuPi_pi_idx[selectedCandIdx_sig]];
    the_pi_dzsig = ProbeTracks_dzS[HNLToMuPi_pi_idx[selectedCandIdx_sig]];

    the_sv_chi2 = HNLToMuPi_sv_chi2[selectedCandIdx_sig];
    the_sv_lxy = HNLToMuPi_sv_lxy[selectedCandIdx_sig];
    the_sv_lxysig = HNLToMuPi_sv_lxy_sig[selectedCandIdx_sig];
    the_sv_prob = HNLToMuPi_sv_prob[selectedCandIdx_sig];

    the_ismatched = HNLToMuPi_isMatched[selectedCandIdx_sig];
    the_mu_ismatched = HNLToMuPi_sel_mu_isMatched[selectedCandIdx_sig];
    the_pi_ismatched = HNLToMuPi_pi_isMatched[selectedCandIdx_sig];
    the_mupi_mass_reco_gen_reldiff = HNLToMuPi_mupi_mass_reco_gen_reldiff[selectedCandIdx_sig];
    the_lxy_reco_gen_reldiff = HNLToMuPi_lxy_reco_gen_reldiff[selectedCandIdx_sig];
  } // end at least one candidate in the event
  else{ // we want to save the gen information, even if no reco candidate was built
    the_hnl_pt = -99.;
    the_hnl_eta = -99.;
    the_hnl_phi = -99.;
    the_hnl_mass = -99.;
    the_hnl_charge = -99;
    the_hnl_cos2d = -99.;

    the_mu_pt = -99.;
    the_mu_eta = -99.;
    the_mu_phi = -99.; 
    the_mu_charge = -99;
    the_mu_dxy = -99.;
    the_mu_dxysig = -99.;
    the_mu_dzsig = -99.;
    the_mu_looseid = -99;
    the_mu_mediumid = -99;
    the_mu_tightid = -99;
    the_mu_softid = -99;
    the_mu_istriggering = -99;
    the_mu_isslimmed = -99;
    the_mu_isdsa = -99;

    the_pi_pt = -99.;
    the_pi_eta = -99.;
    the_pi_phi = -99.; 
    the_pi_charge = -99;
    the_pi_dcasig = -99.;
    the_pi_dxy = -99.;
    the_pi_dz = -99.;
    the_pi_dxysig = -99.;
    the_pi_dzsig = -99.;

    the_sv_chi2 = -99.;
    the_sv_lxy = -99.;
    the_sv_lxysig = -99.;
    the_sv_prob = -99.;

    the_ismatched = -99;
    the_mu_ismatched = -99;
    the_pi_ismatched = -99;
    the_mupi_mass_reco_gen_reldiff = -99.;
    the_lxy_reco_gen_reldiff = -99.;
  }

  // getting the displacement at gen level
  if(isMC){
    UInt_t nGen = *nGenPart;

    // find idx of gen particles of interest
    int gen_hnl_idx(-99), gen_b_idx(-99), gen_trgmu_idx(-99), gen_mu_idx(-99), gen_pi_idx(-99);

    for(unsigned int iGen(0); iGen < nGen; ++iGen){
      if(abs(GenPart_pdgId[iGen])==9900015){
        gen_hnl_idx = iGen;
        gen_b_idx = GenPart_genPartIdxMother[iGen]; 
        break;
      }
    }
    for(unsigned int iGen(0); iGen < nGen; ++iGen){
      if(abs(GenPart_pdgId[iGen])==13 && GenPart_genPartIdxMother[iGen]==gen_b_idx)
        gen_trgmu_idx = iGen;
      if((abs(GenPart_pdgId[iGen])==13) && GenPart_genPartIdxMother[iGen]==gen_hnl_idx)
        gen_mu_idx = iGen;
      if(abs(GenPart_pdgId[iGen])==211 && GenPart_genPartIdxMother[iGen]==gen_hnl_idx)
        gen_pi_idx = iGen;
    }

    // get quantities that need more than one object
    if(gen_hnl_idx!=-99 && gen_mu_idx!=-99){

      ROOT::Math::PtEtaPhiMVector hnl_p4(GenPart_pt[gen_hnl_idx], GenPart_eta[gen_hnl_idx], GenPart_phi[gen_hnl_idx], GenPart_mass[gen_hnl_idx]);
      Float_t hnl_betagamma = hnl_p4.Beta() * hnl_p4.Gamma();

      the_gen_hnl_lxyz = get3Ddisp(GenPart_vx[gen_hnl_idx], GenPart_vx[gen_mu_idx],
          GenPart_vy[gen_hnl_idx], GenPart_vy[gen_mu_idx],
          GenPart_vz[gen_hnl_idx], GenPart_vz[gen_mu_idx]);

      the_gen_hnl_lxy = get2Ddisp(GenPart_vx[gen_hnl_idx], GenPart_vx[gen_mu_idx],
          GenPart_vy[gen_hnl_idx], GenPart_vy[gen_mu_idx]);

      the_gen_hnl_ct = the_gen_hnl_lxyz / hnl_betagamma;
    }

    // set quantities for each object
    if(gen_b_idx!=-99){
      the_gen_b_pt = GenPart_pt[gen_b_idx];
      the_gen_b_eta = GenPart_eta[gen_b_idx];
      the_gen_b_phi = GenPart_phi[gen_b_idx];
      the_gen_b_mass = GenPart_mass[gen_b_idx];
      the_gen_b_pdgid = GenPart_pdgId[gen_b_idx];
    }
    if(gen_hnl_idx!=-99){
      the_gen_hnl_pt = GenPart_pt[gen_hnl_idx];
      the_gen_hnl_eta = GenPart_eta[gen_hnl_idx];
      the_gen_hnl_phi = GenPart_phi[gen_hnl_idx];
      the_gen_hnl_mass = GenPart_mass[gen_hnl_idx];
      the_gen_hnl_vx = GenPart_vx[gen_hnl_idx];
      the_gen_hnl_vy = GenPart_vy[gen_hnl_idx];
      the_gen_hnl_vz = GenPart_vz[gen_hnl_idx];
    }
    if(gen_trgmu_idx!=-99){
      the_gen_trgmu_pt = GenPart_pt[gen_trgmu_idx];
      the_gen_trgmu_eta = GenPart_eta[gen_trgmu_idx];
      the_gen_trgmu_phi = GenPart_phi[gen_trgmu_idx];
      the_gen_trgmu_vx = GenPart_vx[gen_trgmu_idx];
      the_gen_trgmu_vy = GenPart_vy[gen_trgmu_idx];
      the_gen_trgmu_vz = GenPart_vz[gen_trgmu_idx];
    }
    if(gen_mu_idx!=-99){
      the_gen_mu_pt = GenPart_pt[gen_mu_idx];
      the_gen_mu_eta = GenPart_eta[gen_mu_idx];
      the_gen_mu_phi = GenPart_phi[gen_mu_idx];
      the_gen_mu_vx = GenPart_vx[gen_mu_idx];
      the_gen_mu_vy = GenPart_vy[gen_mu_idx];
      the_gen_mu_vz = GenPart_vz[gen_mu_idx];
    }
    if(gen_pi_idx!=-99){
      the_gen_pi_pt = GenPart_pt[gen_pi_idx];
      the_gen_pi_eta = GenPart_eta[gen_pi_idx];
      the_gen_pi_phi = GenPart_phi[gen_pi_idx];
      the_gen_pi_vx = GenPart_vx[gen_pi_idx];
      the_gen_pi_vy = GenPart_vy[gen_pi_idx];
      the_gen_pi_vz = GenPart_vz[gen_pi_idx];
    }
  }

  hnl_tree->Fill();
  //} // end bpark line fired
  //} // end sound index

  return kTRUE;
  }


void HNLToMuPiDumper::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}


void HNLToMuPiDumper::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  my_file->cd();

  hnl_tree->Write("", TObject::kOverwrite);

  my_file->Close();

  cout << "- End HNL->mupi Dumper -" << endl;
}

