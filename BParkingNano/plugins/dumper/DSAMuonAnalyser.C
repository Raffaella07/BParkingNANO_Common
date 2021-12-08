#define DSAMuonAnalyser_cxx
// The class definition in DSAMuonAnalyser.h has been generated automatically
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
// root> T->Process("DSAMuonAnalyser.C")
// root> T->Process("DSAMuonAnalyser.C","some options")
// root> T->Process("DSAMuonAnalyser.C+")
//


#include "DSAMuonAnalyser.h"
#include <TMath.h>
#include <cmath>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "utils.C"
#include "DataFormats/Math/interface/deltaR.h"


using namespace std;


void DSAMuonAnalyser::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  cout << " --------------------------" << endl;
  cout << "      DSA Muon Analyser    " << endl;
  cout << " --------------------------" << endl;
}


void DSAMuonAnalyser::SlaveBegin(TTree * /*tree*/)
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
  my_file = new TFile(outFileName, "RECREATE");  
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

  // defining histograms
  my_file->cd();

  hist_ncand_perevent = new TH1F("hist_ncand_perevent", "hist_ncand_perevent", 10, 0, 10);
  hist_ncand_matched_perevent = new TH1F("hist_ncand_matched_perevent", "hist_ncand_matched_perevent", 10, 0, 10);
  hist_ncand_atleastone_matched_perevent = new TH1F("hist_ncand_atleastone_matched_perevent", "hist_ncand_atleastone_matched_perevent", 10, 0, 10);

  hist_unique_slimmedcand = new TH1F("hist_unique_slimmedcand", "Number of events having only one matched slimmed candidate", 2, 0, 2);
  hist_unique_dsacand = new TH1F("hist_unique_dsacand", "Number of events having only one matched dsa candidate", 2, 0, 2);
  hist_unique_dsacand_displacement = new TH1F("hist_unique_dsacand_displacement", "Displacement of dsa cand", 80, 0, 100);
  hist_unique_dsacand_displacement_matched = new TH1F("hist_unique_dsacand_displacement_matched", "Displacement of dsa cand if matched", 80, 0, 100);
  hist_unique_dsacand_displacement_nonmatched = new TH1F("hist_unique_dsacand_displacement_nonmatched", "Displacement of dsa cand if non-matched", 80, 0, 100);
  hist_unique_dsacand_ismatchedtoslimmed = new TH1F("hist_unique_dsacand_ismatchedtoslimmed", "Probability of the dsa candidate to have its muon matched to a slimmed muon", 2, 0, 2);
  hist_unique_dsacand_deltaR = new TH1F("hist_unique_dsacand_deltaR", "DeltaR between DSA and slimmed in dsa events", 30, 0, 0.3);
  hist_unique_dsacand_deltaPtRel = new TH1F("hist_unique_dsacand_deltaPtRel", "DeltaPtRel between DSA and slimmed in dsa events", 30, 0, 0.4);
  hist_unique_dsacand_deltadxyRel = new TH1F("hist_unique_dsacand_deltadxyRel", "DeltadxyRel between DSA and slimmed in dsa events", 100, 0, 50);
  hist_unique_dsacand_deltadzRel = new TH1F("hist_unique_dsacand_deltadzRel", "DeltadzRel between DSA and slimmed in dsa events", 100, 0, 50);
  hist_unique_dsacand_deltadxySRel = new TH1F("hist_unique_dsacand_deltadxySRel", "DeltadxySRel between DSA and slimmed in mixed events", 100, 0, 5);
  hist_unique_dsacand_deltadzSRel = new TH1F("hist_unique_dsacand_deltadzSRel", "DeltadzSRel between DSA and slimmed in mixed events", 100, 0, 5);
  hist_unique_dsacand_deltavx = new TH1F("hist_unique_dsacand_vx", "Difference in vx between DSA and slimmed in mixed events", 100, 0, 50);
  hist_unique_dsacand_deltavz = new TH1F("hist_unique_dsacand_vz", "Difference in vz between DSA and slimmed in mixed events", 100, 0, 50);
  hist_unique_dsacand_samecharge = new TH1F("hist_unique_dsacand_samecharge", "Probability of the DSA muon to have the same charge as the slimmed muon", 2, 0, 2);
  hist_unique_dsacand_matchingefficiency_lxygt1 = new TH1F("hist_unique_dsacand_matchingefficiency_lxygt1", "Matching efficiency in signal events with dsa displacement greater than 1 cm", 2, 0, 2);
  hist_unique_dsacand_matchingefficiency_lxygt5 = new TH1F("hist_unique_dsacand_matchingefficiency_lxygt5", "Matching efficiency in signal events with dsa displacement greater than 5 cm", 2, 0, 2);
  hist_unique_dsacand_matchingefficiency_lxygt10 = new TH1F("hist_unique_dsacand_matchingefficiency_lxygt10", "Matching efficiency in signal events with dsa displacement greater than 10 cm", 2, 0, 2);
  hist_unique_dsacand_matchingefficiency_lxygt15 = new TH1F("hist_unique_dsacand_matchingefficiency_lxygt15", "Matching efficiency in signal events with dsa displacement greater than 15 cm", 2, 0, 2);
  hist_unique_dsacand_matchingefficiency_lxygt20 = new TH1F("hist_unique_dsacand_matchingefficiency_lxygt20", "Matching efficiency in signal events with dsa displacement greater than 20 cm", 2, 0, 2);
  hist_double_slimmedcand = new TH1F("hist_double_slimmedcand", "Number of events having exactly two matched slimmed candidates", 2, 0, 2);
  hist_double_dsacand = new TH1F("hist_double_dsacand", "Number of events having exactly two matched dsa candidates", 2, 0, 2);
  hist_double_mixedcand = new TH1F("hist_double_mixedcand", "Number of events having one matched slimmed and dsa candidate", 2, 0, 2);
  hist_double_mixedcand_displacement = new TH1F("hist_double_mixedcand_displacement", "Displacement of dsa cand", 80, 0, 100);
  hist_double_mixedcand_displacement_matched = new TH1F("hist_double_mixedcand_displacement_matched", "Displacement of dsa cand if matched", 80, 0, 100);
  hist_double_mixedcand_samepion = new TH1F("hist_double_mixedcand_samepion", "Probability that the two candidates have the same pion", 2, 0, 2);
  hist_double_mixedcand_sametrgmu = new TH1F("hist_double_mixedcand_sametrgmu", "Probability that the two candidates have the same trigger muon", 2, 0, 2);
  hist_double_mixedcand_deltaR = new TH1F("hist_double_mixedcand_deltaR", "DeltaR between DSA and slimmed in mixed events", 30, 0, 0.3);
  hist_double_mixedcand_deltaPtRel = new TH1F("hist_double_mixedcand_deltaPtRel", "DeltaPtRel between DSA and slimmed in mixed events", 30, 0, 0.4);
  hist_double_mixedcand_deltadxyRel = new TH1F("hist_double_mixedcand_deltadxyRel", "DeltadxyRel between DSA and slimmed in mixed events", 100, 0, 50);
  hist_double_mixedcand_deltadzRel = new TH1F("hist_double_mixedcand_deltadzRel", "DeltadzRel between DSA and slimmed in mixed events", 100, 0, 50);
  hist_double_mixedcand_deltadxySRel = new TH1F("hist_double_mixedcand_deltadxySRel", "DeltadxySRel between DSA and slimmed in mixed events", 100, 0.9, 1.1);
  hist_double_mixedcand_deltadzSRel = new TH1F("hist_double_mixedcand_deltadzSRel", "DeltadzSRel between DSA and slimmed in mixed events", 100, 0.9, 1.1);
  hist_double_mixedcand_deltavx = new TH1F("hist_double_mixedcand_vx", "Difference in vx between DSA and slimmed in mixed events", 100, 0, 50);
  hist_double_mixedcand_deltavz = new TH1F("hist_double_mixedcand_vz", "Difference in vz between DSA and slimmed in mixed events", 100, 0, 50);
  hist_double_mixedcand_samecharge = new TH1F("hist_double_mixedcand_samecharge", "Probability of the DSA muon to have the same charge as the slimmed muon", 2, 0, 2);
  hist_double_mixedcand_ismatchedtoslimmed = new TH1F("hist_double_mixedcand_ismatchedtoslimmed", "Probability of the dsa candidate to have its muon matched to a slimmed muon", 2, 0, 2);
  hist_double_mixedcand_ismatchedtocorrectslimmed = new TH1F("hist_double_mixedcand_ismatchedtocorrectslimmed", "Probability for the dsa candidate to have its muon matched to the slimmed muon of the slimmed candidate", 2, 0, 2);
  hist_double_mixedcand_matchingefficiency_lxygt1 = new TH1F("hist_double_mixedcand_matchingefficiency_lxygt1", "Matching efficiency in signal events with dsa displacement greater than 1 cm", 2, 0, 2);
  hist_double_mixedcand_matchingefficiency_lxygt5 = new TH1F("hist_double_mixedcand_matchingefficiency_lxygt5", "Matching efficiency in signal events with dsa displacement greater than 5 cm", 2, 0, 2);
  hist_double_mixedcand_matchingefficiency_lxygt10 = new TH1F("hist_double_mixedcand_matchingefficiency_lxygt10", "Matching efficiency in signal events with dsa displacement greater than 10 cm", 2, 0, 2);
  hist_double_mixedcand_matchingefficiency_lxygt15 = new TH1F("hist_double_mixedcand_matchingefficiency_lxygt15", "Matching efficiency in signal events with dsa displacement greater than 15 cm", 2, 0, 2);
  hist_double_mixedcand_matchingefficiency_lxygt20 = new TH1F("hist_double_mixedcand_matchingefficiency_lxygt20", "Matching efficiency in signal events with dsa displacement greater than 20 cm", 2, 0, 2);
  hist_triple_oneslimmed_twodsa = new TH1F("hist_triple_oneslimmed_twodsa", "Number of events having one slimmed and two dsa candidates", 2, 0, 2);
}


Bool_t DSAMuonAnalyser::Process(Long64_t entry)
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

  // number of candidates in the event
  UInt_t nCand_sig = *nBToMuMuPi; 

  //   ----- Histograms -----  //

  // number of matched candidates in the event (only saved if nCand non zero)
  UInt_t nMatchedCand_sig = 0;

  hist_ncand_perevent->Fill(nCand_sig);

  // get index of the matched candidates
  vector<int> matched_idx;
  if(nCand_sig > 0){
    for(unsigned int iCand(0); iCand < nCand_sig; ++iCand){
      if(BToMuMuPi_isMatched[iCand] == 1){
        ++nMatchedCand_sig;
        matched_idx.push_back(iCand);
      }
    }

    hist_ncand_matched_perevent->Fill(nMatchedCand_sig);
    if(nMatchedCand_sig>0) hist_ncand_atleastone_matched_perevent->Fill(nMatchedCand_sig);

    // fetch information on what kind of candidates are built (slimmed, DSA, slimmed-slimmed, DSA-DSA, slimmed-DSA)
    if(nMatchedCand_sig==1 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[0]]]==0) hist_unique_slimmedcand->Fill(1.);
    if(nMatchedCand_sig==1 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[0]]]==1) hist_unique_dsacand->Fill(1.);
    if(nMatchedCand_sig==2 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[0]]]==0 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[1]]]==0) hist_double_slimmedcand->Fill(1.);
    if(nMatchedCand_sig==2 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[0]]]==1 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[1]]]==1) hist_double_dsacand->Fill(1.);

    // check how much we would loose if not considering events with 3 candidates (conclusion: events with more than 2 candidates can be ignored)
    if(nMatchedCand_sig==3 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[0]]]==0 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[1]]]==1 && 
        Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[2]]]==1) hist_triple_oneslimmed_twodsa->Fill(1.);

    bool mixed_candidate = false;
    // by construction the dsa candidate has to come in the second position
    if(nMatchedCand_sig==2 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[0]]]==0 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[1]]]==1){
      hist_double_mixedcand->Fill(1.);
      mixed_candidate = true;
    }

    // matching criterias
    float max_deltaR = 0.1;
    float max_deltaPtRel = 0.2;

    // study how likely it is that the DSA is matched to the slimmed muon in the mixed case
    float deltaR = -99.;
    float deltaPtRel = -99;
    float deltadxyRel = -99;
    float deltadzRel = -99;
    float deltadxySRel = -99;
    float deltadzSRel = -99;
    float deltavx = -99;
    float deltavz = -99;
    int slimmed_idx = -1;
    int dsa_idx = -1;
    if(mixed_candidate){
      slimmed_idx = BToMuMuPi_sel_mu_idx[matched_idx[0]];
      dsa_idx = BToMuMuPi_sel_mu_idx[matched_idx[1]];

      hist_double_mixedcand_displacement->Fill(BToMuMuPi_sv_lxy[matched_idx[1]]);;

      // check that trigger muon and pion are the same in the two candidates
      if(BToMuMuPi_pi_idx[matched_idx[0]] == BToMuMuPi_pi_idx[matched_idx[1]]){
        hist_double_mixedcand_samepion->Fill(1.);
      }
      else{
        hist_double_mixedcand_samepion->Fill(0.);
      }
      if(BToMuMuPi_trg_mu_idx[matched_idx[0]] == BToMuMuPi_trg_mu_idx[matched_idx[1]]){
        hist_double_mixedcand_sametrgmu->Fill(1.);
      }
      else{
        hist_double_mixedcand_sametrgmu->Fill(0.);
      }

      // fetch the relative info between the slimmed and DSA muon
      deltaR = reco::deltaR(Muon_eta[dsa_idx], Muon_phi[dsa_idx], Muon_eta[slimmed_idx], Muon_phi[slimmed_idx]);
      deltaPtRel = fabs(Muon_pt[dsa_idx] - Muon_pt[slimmed_idx]) / Muon_pt[slimmed_idx];
      deltadxyRel = fabs(Muon_dxy[dsa_idx] - Muon_dxy[slimmed_idx]) / Muon_dxy[slimmed_idx];
      deltadzRel = fabs(Muon_dz[dsa_idx] - Muon_dz[slimmed_idx]) / Muon_dz[slimmed_idx];
      deltadxySRel = fabs(Muon_dxyS[dsa_idx] - Muon_dxyS[slimmed_idx]) / Muon_dxyS[slimmed_idx];
      deltadzSRel = fabs(Muon_dzS[dsa_idx] - Muon_dzS[slimmed_idx]) / Muon_dzS[slimmed_idx];
      deltavx = fabs(Muon_vx[dsa_idx] - Muon_vx[slimmed_idx]);
      deltavz = fabs(Muon_vz[dsa_idx] - Muon_vz[slimmed_idx]);
      //if(deltaR < 0.1 && deltaPtRel < 0.2){
      //if(BToMuMuPi_sv_lxy[matched_idx[1]]>20){
        hist_double_mixedcand_deltaR->Fill(deltaR);
        hist_double_mixedcand_deltaPtRel->Fill(deltaPtRel);
        hist_double_mixedcand_deltadxyRel->Fill(deltadxyRel);
        hist_double_mixedcand_deltadzRel->Fill(deltadzRel);
        hist_double_mixedcand_deltadxySRel->Fill(deltadxySRel);
        hist_double_mixedcand_deltadzSRel->Fill(deltadzSRel);
        hist_double_mixedcand_deltavx->Fill(deltavx);
        hist_double_mixedcand_deltavz->Fill(deltavz);
      //}

      // check if the DSA muon is matched to the slimmed muon of the slimmed candidate
      bool matched_toslimmed = false;
      std::vector<pair<int, float>> pairs_slimmedIdx_deltaPtRel;
      for(unsigned int iMuon(0); iMuon<*nMuon; ++iMuon){
        if(Muon_isDSAMuon[iMuon]==1) continue;
        float deltaR = reco::deltaR(Muon_eta[dsa_idx], Muon_phi[dsa_idx], Muon_eta[iMuon], Muon_phi[iMuon]);
        if(deltaR < max_deltaR){
          pair<int, float> pairs_slimmedIdx_deltaPtRel_tmp;
          pairs_slimmedIdx_deltaPtRel_tmp.first = iMuon;
          pairs_slimmedIdx_deltaPtRel_tmp.second = fabs(Muon_pt[dsa_idx] - Muon_pt[iMuon]) / Muon_pt[iMuon];
          pairs_slimmedIdx_deltaPtRel.push_back(pairs_slimmedIdx_deltaPtRel_tmp);
        }
      }
      sort(pairs_slimmedIdx_deltaPtRel.begin(), pairs_slimmedIdx_deltaPtRel.end(), [](const pair<int, float> &pair_i, const pair<int, float> &pair_j){
        return pair_i.second < pair_j.second;
      });

      if(pairs_slimmedIdx_deltaPtRel.size() > 0 && pairs_slimmedIdx_deltaPtRel[0].second < max_deltaPtRel){
        matched_toslimmed = true;
        hist_double_mixedcand_ismatchedtoslimmed->Fill(1.);
        int match_idx = pairs_slimmedIdx_deltaPtRel[0].first;
        if(matched_toslimmed && match_idx == slimmed_idx) hist_double_mixedcand_ismatchedtocorrectslimmed->Fill(1.);
        else hist_double_mixedcand_ismatchedtocorrectslimmed->Fill(0.);

        if(Muon_charge[dsa_idx] == Muon_charge[match_idx]) hist_double_mixedcand_samecharge->Fill(1.);
        else hist_double_mixedcand_samecharge->Fill(0.);
        hist_double_mixedcand_displacement_matched->Fill(BToMuMuPi_sv_lxy[matched_idx[1]]);;
      }
      else{
        hist_double_mixedcand_ismatchedtoslimmed->Fill(0.);
      }

      // study the efficiency in different displacement bins
      float dsa_displacement = BToMuMuPi_sv_lxy[matched_idx[1]];
      if(dsa_displacement > 1){
        if(matched_toslimmed){
          hist_double_mixedcand_matchingefficiency_lxygt1->Fill(1.);
        }
        else{
          hist_double_mixedcand_matchingefficiency_lxygt1->Fill(0.);
        }
      }
      if(dsa_displacement > 5){
        if(matched_toslimmed){
          hist_double_mixedcand_matchingefficiency_lxygt5->Fill(1.);
        }
        else{
          hist_double_mixedcand_matchingefficiency_lxygt5->Fill(0.);
        }
      }
      if(dsa_displacement > 10){
        if(matched_toslimmed){
          hist_double_mixedcand_matchingefficiency_lxygt10->Fill(1.);
        }
        else{
          hist_double_mixedcand_matchingefficiency_lxygt10->Fill(0.);
        }
      }
      if(dsa_displacement > 15){
        if(matched_toslimmed){
          hist_double_mixedcand_matchingefficiency_lxygt15->Fill(1.);
        }
        else{
          hist_double_mixedcand_matchingefficiency_lxygt15->Fill(0.);
        }
      }
      if(dsa_displacement > 20){
        if(matched_toslimmed){
          hist_double_mixedcand_matchingefficiency_lxygt20->Fill(1.);
        }
        else{
          hist_double_mixedcand_matchingefficiency_lxygt20->Fill(0.);
        }
      }
    }

    // now going on the side of the presumably non-matched DSA, i.e events with no slimmed candidate
    // check the matching efficiency of those DSA candidates
    if(nMatchedCand_sig==1 && Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[matched_idx[0]]]==1){
      dsa_idx = BToMuMuPi_sel_mu_idx[matched_idx[0]];
      hist_unique_dsacand_displacement->Fill(BToMuMuPi_sv_lxy[matched_idx[0]]);;

      bool ismatched = false;
      std::vector<pair<int, float>> pairs_slimmedIdx_deltaPtRel;
      for(unsigned int iMuon(0); iMuon<*nMuon; ++iMuon){
        if(Muon_isDSAMuon[iMuon]==1) continue;
        float deltaR = reco::deltaR(Muon_eta[dsa_idx], Muon_phi[dsa_idx], Muon_eta[iMuon], Muon_phi[iMuon]);
        if(deltaR < max_deltaR){
          pair<int, float> pairs_slimmedIdx_deltaPtRel_tmp;
          pairs_slimmedIdx_deltaPtRel_tmp.first = iMuon;
          pairs_slimmedIdx_deltaPtRel_tmp.second = fabs(Muon_pt[dsa_idx] - Muon_pt[iMuon]) / Muon_pt[iMuon];
          pairs_slimmedIdx_deltaPtRel.push_back(pairs_slimmedIdx_deltaPtRel_tmp);
        }
      }
      sort(pairs_slimmedIdx_deltaPtRel.begin(), pairs_slimmedIdx_deltaPtRel.end(), [](const pair<int, float> &pair_i, const pair<int, float> &pair_j){
        return pair_i.second < pair_j.second;
      });
      if(pairs_slimmedIdx_deltaPtRel.size() > 0 && pairs_slimmedIdx_deltaPtRel[0].second < max_deltaPtRel){
        ismatched = true;
        hist_unique_dsacand_ismatchedtoslimmed->Fill(1.);
        int match_idx = pairs_slimmedIdx_deltaPtRel[0].first;
        deltaR = reco::deltaR(Muon_eta[dsa_idx], Muon_phi[dsa_idx], Muon_eta[match_idx], Muon_phi[match_idx]);
        deltaPtRel = fabs(Muon_pt[dsa_idx] - Muon_pt[match_idx]) / Muon_pt[match_idx];
        deltadxyRel = fabs(Muon_dxy[dsa_idx] - Muon_dxy[match_idx]) / Muon_dxy[match_idx];
        deltadzRel = fabs(Muon_dz[dsa_idx] - Muon_dz[match_idx]) / Muon_dz[match_idx];
        deltadxySRel = fabs(Muon_dxyS[dsa_idx] - Muon_dxyS[match_idx]) / Muon_dxyS[match_idx];
        deltadzSRel = fabs(Muon_dzS[dsa_idx] - Muon_dzS[match_idx]) / Muon_dzS[match_idx];
        deltavx = fabs(Muon_vx[dsa_idx] - Muon_vx[match_idx]);
        deltavz = fabs(Muon_vz[dsa_idx] - Muon_vz[match_idx]);
        //if(deltaR < 0.1 && deltaPtRel < 0.2){
        hist_unique_dsacand_deltaR->Fill(deltaR);
        hist_unique_dsacand_deltaPtRel->Fill(deltaPtRel);
        hist_unique_dsacand_deltadxyRel->Fill(deltadxyRel);
        hist_unique_dsacand_deltadzRel->Fill(deltadzRel);
        hist_unique_dsacand_deltadxySRel->Fill(deltadxySRel);
        hist_unique_dsacand_deltadzSRel->Fill(deltadzSRel);
        hist_unique_dsacand_deltavx->Fill(deltavx);
        hist_unique_dsacand_deltavz->Fill(deltavz);
        hist_unique_dsacand_displacement_matched->Fill(BToMuMuPi_sv_lxy[matched_idx[0]]);;
        if(Muon_charge[dsa_idx] == Muon_charge[match_idx]) hist_unique_dsacand_samecharge->Fill(1.);
        else hist_unique_dsacand_samecharge->Fill(0.);
        //}
      }
      else{
        hist_unique_dsacand_ismatchedtoslimmed->Fill(0.);
        hist_unique_dsacand_displacement_nonmatched->Fill(BToMuMuPi_sv_lxy[matched_idx[0]]);;
      }
      // study matching efficiency in bins of displacement
      float displacement = BToMuMuPi_sv_lxy[matched_idx[0]];
      if(displacement > 1){
        if(ismatched){
          hist_unique_dsacand_matchingefficiency_lxygt1->Fill(1.);
        }
        else{
          hist_unique_dsacand_matchingefficiency_lxygt1->Fill(0.);
        }
      }
      if(displacement > 5){
        if(ismatched){
          hist_unique_dsacand_matchingefficiency_lxygt5->Fill(1.);
        }
        else{
          hist_unique_dsacand_matchingefficiency_lxygt5->Fill(0.);
        }
      }
      if(displacement > 10){
        if(ismatched){
          hist_unique_dsacand_matchingefficiency_lxygt10->Fill(1.);
        }
        else{
          hist_unique_dsacand_matchingefficiency_lxygt10->Fill(0.);
        }
      }
      if(displacement > 15){
        if(ismatched){
          hist_unique_dsacand_matchingefficiency_lxygt15->Fill(1.);
        }
        else{
          hist_unique_dsacand_matchingefficiency_lxygt15->Fill(0.);
        }
      }
      if(displacement > 20){
        if(ismatched){
          hist_unique_dsacand_matchingefficiency_lxygt20->Fill(1.);
        }
        else{
          hist_unique_dsacand_matchingefficiency_lxygt20->Fill(0.);
        }
      }
    }
  } // end at least one candidate

  return kTRUE;
}


void DSAMuonAnalyser::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}


void DSAMuonAnalyser::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  my_file->cd();
  my_file->Write();
  my_file->Close();

  cout << "- End DSA Muon Analyser -" << endl;
}

