#define BToKMuMuDumper_cxx
// The class definition in BToKMuMuDumper.h has been generated automatically
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
// root> T->Process("BToKMuMuDumper.C")
// root> T->Process("BToKMuMuDumper.C","some options")
// root> T->Process("BToKMuMuDumper.C+")
//


#include "BToKMuMuDumper.h"
#include <TMath.h>
#include <cmath>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "utils.C"


using namespace std;


void BToKMuMuDumper::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  cout << " --------------------------" << endl;
  cout << "       B->KMuMu Dumper     " << endl;
  cout << " --------------------------" << endl;
}


void BToKMuMuDumper::SlaveBegin(TTree * /*tree*/)
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

  // getting the control tree ready
  control_tree = new TTree("control_tree", "control_tree");

  control_tree->Branch("event", &the_event);
  control_tree->Branch("run", &the_run);
  control_tree->Branch("lumi", &the_lumi);

  control_tree->Branch("b_pt", &the_ctrl_b_pt);
  control_tree->Branch("b_eta", &the_ctrl_b_eta);
  control_tree->Branch("b_phi", &the_ctrl_b_phi);
  control_tree->Branch("b_mass", &the_ctrl_b_mass);
  control_tree->Branch("b_charge", &the_ctrl_b_charge);
  control_tree->Branch("b_pdgid", &the_ctrl_b_pdgid);
  control_tree->Branch("b_cos2d", &the_ctrl_b_cos2d);
  control_tree->Branch("b_iso03", &the_ctrl_b_iso03);
  control_tree->Branch("b_iso03_close", &the_ctrl_b_iso03_close);
  control_tree->Branch("b_iso04", &the_ctrl_b_iso04);
  control_tree->Branch("b_iso04_close", &the_ctrl_b_iso04_close);

  control_tree->Branch("k_pt", &the_ctrl_k_pt);
  control_tree->Branch("k_eta", &the_ctrl_k_eta);
  control_tree->Branch("k_phi", &the_ctrl_k_phi);
  control_tree->Branch("k_charge", &the_ctrl_k_charge);
  control_tree->Branch("k_iso03", &the_ctrl_k_iso03);
  control_tree->Branch("k_iso03_close", &the_ctrl_k_iso03_close);
  control_tree->Branch("k_iso04", &the_ctrl_k_iso04);
  control_tree->Branch("k_iso04_close", &the_ctrl_k_iso04_close);

  control_tree->Branch("l1_pt", &the_ctrl_l1_pt);
  control_tree->Branch("l1_eta", &the_ctrl_l1_eta);
  control_tree->Branch("l1_phi", &the_ctrl_l1_phi);
  control_tree->Branch("l1_charge", &the_ctrl_l1_charge);
  control_tree->Branch("l1_iso03", &the_ctrl_l1_iso03);
  control_tree->Branch("l1_iso03_close", &the_ctrl_l1_iso03_close);
  control_tree->Branch("l1_iso04", &the_ctrl_l1_iso04);
  control_tree->Branch("l1_iso04_close", &the_ctrl_l1_iso04_close);
  control_tree->Branch("l1_istriggering", &the_ctrl_l1_istriggering);
  control_tree->Branch("l1_fired_hlt_mu7_ip4", &the_ctrl_l1_fired_hlt_mu7_ip4);
  control_tree->Branch("l1_fired_hlt_mu8_ip3", &the_ctrl_l1_fired_hlt_mu8_ip3);
  control_tree->Branch("l1_fired_hlt_mu8_ip5", &the_ctrl_l1_fired_hlt_mu8_ip5);
  control_tree->Branch("l1_fired_hlt_mu8_ip6", &the_ctrl_l1_fired_hlt_mu8_ip6);
  control_tree->Branch("l1_fired_hlt_mu8p5_ip3p5", &the_ctrl_l1_fired_hlt_mu8p5_ip3p5);
  control_tree->Branch("l1_fired_hlt_mu9_ip4", &the_ctrl_l1_fired_hlt_mu9_ip4);
  control_tree->Branch("l1_fired_hlt_mu9_ip5", &the_ctrl_l1_fired_hlt_mu9_ip5);
  control_tree->Branch("l1_fired_hlt_mu9_ip6", &the_ctrl_l1_fired_hlt_mu9_ip6);
  control_tree->Branch("l1_fired_hlt_mu10p5_ip3p5", &the_ctrl_l1_fired_hlt_mu10p5_ip3p5);
  control_tree->Branch("l1_fired_hlt_mu12_ip6", &the_ctrl_l1_fired_hlt_mu12_ip6);


  control_tree->Branch("l2_pt", &the_ctrl_l2_pt);
  control_tree->Branch("l2_eta", &the_ctrl_l2_eta);
  control_tree->Branch("l2_phi", &the_ctrl_l2_phi);
  control_tree->Branch("l2_charge", &the_ctrl_l2_charge);
  control_tree->Branch("l2_iso03", &the_ctrl_l2_iso03);
  control_tree->Branch("l2_iso03_close", &the_ctrl_l2_iso03_close);
  control_tree->Branch("l2_iso04", &the_ctrl_l2_iso04);
  control_tree->Branch("l2_iso04_close", &the_ctrl_l2_iso04_close);
  control_tree->Branch("l2_istriggering", &the_ctrl_l2_istriggering);

  control_tree->Branch("dimu_mass", &the_ctrl_dimu_mass);

  control_tree->Branch("sv_x", &the_ctrl_sv_x);
  control_tree->Branch("sv_y", &the_ctrl_sv_y);
  control_tree->Branch("sv_z", &the_ctrl_sv_z);
  control_tree->Branch("sv_lxy", &the_ctrl_sv_lxy);
  control_tree->Branch("sv_lxysig", &the_ctrl_sv_lxysig);
  control_tree->Branch("sv_prob", &the_ctrl_sv_prob);

  control_tree->Branch("ismatched", &the_ctrl_ismatched);

  control_tree->Branch("pv_npvs", &the_pv_npvs);

  control_tree->Branch("weight_hlt", &the_ctrl_weight_hlt);

  control_tree->Branch("hlt_mu7_ip4", &the_hlt_mu7_ip4);
  control_tree->Branch("hlt_mu8_ip6", &the_hlt_mu8_ip6);
  control_tree->Branch("hlt_mu8_ip5", &the_hlt_mu8_ip5);
  control_tree->Branch("hlt_mu8_ip3", &the_hlt_mu8_ip3);
  control_tree->Branch("hlt_mu8p5_ip3p5", &the_hlt_mu8p5_ip3p5);
  control_tree->Branch("hlt_mu9_ip6", &the_hlt_mu9_ip6);
  control_tree->Branch("hlt_mu9_ip5", &the_hlt_mu9_ip5);
  control_tree->Branch("hlt_mu9_ip4", &the_hlt_mu9_ip4);
  control_tree->Branch("hlt_mu10p5_ip3p5", &the_hlt_mu10p5_ip3p5);
  control_tree->Branch("hlt_mu12_ip6", &the_hlt_mu12_ip6);

  // defining histograms
  if(do_fillhistograms){
    my_file->mkdir("control_channel");
    my_file->cd("control_channel");

    ctrlhist_ncand_perevent = new TH1F("ctrlhist_ncand_perevent", "ctrlhist_ncand_perevent", 10, 0, 10);
    ctrlhist_ncand_matched_perevent = new TH1F("ctrlhist_ncand_matched_perevent", "ctrlhist_ncand_matched__perevent", 10, 0, 10);
    ctrlhist_selection_efficiency_bpt_allevents = new TH1F("ctrlhist_selection_efficiency_bpt_allevents", "Efficiency of selection of candidate with largest b pT (all events)", 2, 0, 2);
    ctrlhist_selection_efficiency_bpt_eventswithmultcands = new TH1F("ctrlhist_selection_efficiency_bpt_eventswithmultcands", "Efficiency of selection of candidate with largest b pT (events with multiple candidates)", 2, 0, 2);
    ctrlhist_selection_efficiency_svprob_allevents = new TH1F("ctrlhist_selection_efficiency_svprob_allevents", "Efficiency of selection of candidate with largest vertex prob. (all events)", 2, 0, 2);
    ctrlhist_selection_efficiency_svprob_eventswithmultcands = new TH1F("ctrlhist_selection_efficiency_svprob_eventswithmultcands", "Efficiency of selection of candidate with largest vertex prob. (events with multiple candidates)", 2, 0, 2);
    ctrlhist_selection_efficiency_cos2d_allevents = new TH1F("ctrlhist_selection_efficiency_cos2d_allevents", "Efficiency of selection of candidate with largest vertex cos(angle) (all events)", 2, 0, 2);
    ctrlhist_selection_efficiency_cos2d_eventswithmultcands = new TH1F("ctrlhist_selection_efficiency_cos2d_eventswithmultcands", "Efficiency of selection of candidate with largest vertex cos(angle) (events with multiple candidates)", 2, 0, 2);
    ctrlhist_selection_efficiency_kpt_allevents = new TH1F("ctrlhist_selection_efficiency_kpt_allevents", "Efficiency of selection of candidate with largest k pT (all events)", 2, 0, 2);
    ctrlhist_selection_efficiency_kpt_eventswithmultcands = new TH1F("ctrlhist_selection_efficiency_kpt_eventswithmultcands", "Efficiency of selection of candidate with largest k pT (events with multiple candidates)", 2, 0, 2);

    ctrlhist_ncand_wtriggeringmuon = new TH1F("ctrlhist_ncand_wtriggeringmuon", "ctrlhist_ncand_wtriggeringmuon", 2, 0, 2);

    my_file->cd();
  } // end define histograms
}


Bool_t BToKMuMuDumper::Process(Long64_t entry)
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
  UInt_t nCand_ctrl = *nBToKMuMu; 

  // branches common to signal and control channels
  the_event = *event; 
  the_run = *run;
  the_lumi = *luminosityBlock; 

  the_pv_npvs = *PV_npvs;

  the_hlt_mu7_ip4 = *HLT_Mu7_IP4;
  the_hlt_mu8_ip6 = *HLT_Mu8_IP6;
  the_hlt_mu8_ip5 = *HLT_Mu8_IP5;
  the_hlt_mu8_ip3 = *HLT_Mu8_IP3;
  the_hlt_mu8p5_ip3p5 = *HLT_Mu8p5_IP3p5;
  the_hlt_mu9_ip6 = *HLT_Mu9_IP6;
  the_hlt_mu9_ip5 = *HLT_Mu9_IP5;
  the_hlt_mu9_ip4 = *HLT_Mu9_IP4;
  the_hlt_mu10p5_ip3p5 = *HLT_Mu10p5_IP3p5;
  the_hlt_mu12_ip6 = *HLT_Mu12_IP6;


  //   ----- Control Channel -----  //

  int ncand_ctrl(0);
  int ncand_istriggering(0);

  if(nCand_ctrl > 0){ // at least one candidate per event
    ncand_ctrl = 1;

    // selecting the candidate as the one having the largest b pt
    // - create candIdx - b pt pair
    vector<pair<int,float>> pair_candIdx_desc_bpt_ctrl = createPairWithDesc(nCand_ctrl, BToKMuMu_fit_pt);
    // - sort pair with decreasing b pt
    stable_sort(pair_candIdx_desc_bpt_ctrl.begin(), pair_candIdx_desc_bpt_ctrl.end(), sortcansbydesc);
    // - for signal, priviledge matched candidates
    vector<pair<int,float>> pair_candIdx_desc_bpt_matched_ctrl = updatePairWithDesc(pair_candIdx_desc_bpt_ctrl, BToKMuMu_isMatched);
    stable_sort(pair_candIdx_desc_bpt_matched_ctrl.begin(), pair_candIdx_desc_bpt_matched_ctrl.end(), sortcansbydesc);
    // - fetch candIdx associated to the largest b pt and the best matching
    UInt_t selectedCandIdx_ctrl = pair_candIdx_desc_bpt_matched_ctrl[0].first;

    // temporary - we manually ask l1 to be the triggering muon
    if(Muon_isTriggeringBPark[BToKMuMu_l1Idx[selectedCandIdx_ctrl]] == 1){

      ncand_istriggering = 1; 

      // fill the control_tree
      the_ctrl_b_pt = BToKMuMu_fit_pt[selectedCandIdx_ctrl];
      the_ctrl_b_eta = BToKMuMu_fit_eta[selectedCandIdx_ctrl];
      the_ctrl_b_phi = BToKMuMu_fit_phi[selectedCandIdx_ctrl];
      the_ctrl_b_mass = BToKMuMu_fit_mass[selectedCandIdx_ctrl];
      the_ctrl_b_charge = BToKMuMu_charge[selectedCandIdx_ctrl];
      the_ctrl_b_pdgid = BToKMuMu_pdgId[selectedCandIdx_ctrl];
      the_ctrl_b_cos2d = BToKMuMu_fit_cos2D[selectedCandIdx_ctrl];
      the_ctrl_b_iso03 = BToKMuMu_b_iso03[selectedCandIdx_ctrl];
      the_ctrl_b_iso03_close = BToKMuMu_b_iso03_close[selectedCandIdx_ctrl];
      the_ctrl_b_iso04 = BToKMuMu_b_iso04[selectedCandIdx_ctrl];
      the_ctrl_b_iso04_close = BToKMuMu_b_iso04_close[selectedCandIdx_ctrl];

      the_ctrl_k_pt = BToKMuMu_fit_k_pt[selectedCandIdx_ctrl];
      the_ctrl_k_eta = BToKMuMu_fit_k_eta[selectedCandIdx_ctrl];
      the_ctrl_k_phi = BToKMuMu_fit_k_phi[selectedCandIdx_ctrl];
      the_ctrl_k_charge = ProbeTracks_charge[BToKMuMu_kIdx[selectedCandIdx_ctrl]];
      the_ctrl_k_iso03 = BToKMuMu_k_iso03[selectedCandIdx_ctrl];
      the_ctrl_k_iso03_close = BToKMuMu_k_iso03_close[selectedCandIdx_ctrl];
      the_ctrl_k_iso04 = BToKMuMu_k_iso04[selectedCandIdx_ctrl];
      the_ctrl_k_iso04_close = BToKMuMu_k_iso04_close[selectedCandIdx_ctrl];

      the_ctrl_l1_pt = BToKMuMu_fit_l1_pt[selectedCandIdx_ctrl];
      the_ctrl_l1_eta = BToKMuMu_fit_l1_eta[selectedCandIdx_ctrl];
      the_ctrl_l1_phi = BToKMuMu_fit_l1_phi[selectedCandIdx_ctrl];
      the_ctrl_l1_charge = Muon_charge[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_iso03 = BToKMuMu_l1_iso03[selectedCandIdx_ctrl];
      the_ctrl_l1_iso03_close = BToKMuMu_l1_iso03_close[selectedCandIdx_ctrl];
      the_ctrl_l1_iso04 = BToKMuMu_l1_iso04[selectedCandIdx_ctrl];
      the_ctrl_l1_iso04_close = BToKMuMu_l1_iso04_close[selectedCandIdx_ctrl];
      the_ctrl_l1_istriggering = Muon_isTriggeringBPark[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_fired_hlt_mu7_ip4 = Muon_fired_HLT_Mu7_IP4[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_fired_hlt_mu8_ip3 = Muon_fired_HLT_Mu8_IP3[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_fired_hlt_mu8_ip5 = Muon_fired_HLT_Mu8_IP5[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_fired_hlt_mu8_ip6 = Muon_fired_HLT_Mu8_IP6[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_fired_hlt_mu8p5_ip3p5 = Muon_fired_HLT_Mu8p5_IP3p5[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_fired_hlt_mu9_ip4 = Muon_fired_HLT_Mu9_IP4[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_fired_hlt_mu9_ip5 = Muon_fired_HLT_Mu9_IP5[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_fired_hlt_mu9_ip6 = Muon_fired_HLT_Mu9_IP6[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_fired_hlt_mu10p5_ip3p5 = Muon_fired_HLT_Mu10p5_IP3p5[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];
      the_ctrl_l1_fired_hlt_mu12_ip6 = Muon_fired_HLT_Mu12_IP6[BToKMuMu_l1Idx[selectedCandIdx_ctrl]];

      the_ctrl_l2_pt = BToKMuMu_fit_l2_pt[selectedCandIdx_ctrl];
      the_ctrl_l2_eta = BToKMuMu_fit_l2_eta[selectedCandIdx_ctrl];
      the_ctrl_l2_phi = BToKMuMu_fit_l2_phi[selectedCandIdx_ctrl];
      the_ctrl_l2_charge = Muon_charge[BToKMuMu_l2Idx[selectedCandIdx_ctrl]];
      the_ctrl_l2_iso03 = BToKMuMu_l2_iso03[selectedCandIdx_ctrl];
      the_ctrl_l2_iso03_close = BToKMuMu_l2_iso03_close[selectedCandIdx_ctrl];
      the_ctrl_l2_iso04 = BToKMuMu_l2_iso04[selectedCandIdx_ctrl];
      the_ctrl_l2_iso04_close = BToKMuMu_l2_iso04_close[selectedCandIdx_ctrl];
      the_ctrl_l2_istriggering = Muon_isTriggeringBPark[BToKMuMu_l2Idx[selectedCandIdx_ctrl]];

      the_ctrl_dimu_mass = BToKMuMu_mll_fullfit[selectedCandIdx_ctrl];
      the_ctrl_sv_x = BToKMuMu_vtx_x[selectedCandIdx_ctrl];
      the_ctrl_sv_y = BToKMuMu_vtx_x[selectedCandIdx_ctrl];
      the_ctrl_sv_z = BToKMuMu_vtx_x[selectedCandIdx_ctrl];
      the_ctrl_sv_lxy = BToKMuMu_l_xy[selectedCandIdx_ctrl];
      the_ctrl_sv_lxysig = BToKMuMu_l_xy[selectedCandIdx_ctrl]/BToKMuMu_l_xy_unc[selectedCandIdx_ctrl];
      the_ctrl_sv_prob = BToKMuMu_svprob[selectedCandIdx_ctrl];

      the_ctrl_ismatched = BToKMuMu_isMatched[selectedCandIdx_ctrl];

      // trigger scale factor
      the_ctrl_weight_hlt = isMC ? getTriggerScaleFactor(the_ctrl_l2_pt, fabs(the_ctrl_l2_eta)) : 1.;

      control_tree->Fill();
    } // l1 is triggering
  }// end at least one candidate in the event


  //   ----- Histograms -----  //

  if(do_fillhistograms){

    // number of matched candidates in the event (only saved if nCand non zero)
    UInt_t nMatchedCand_sig = 0;
    UInt_t nMatchedCand_ctrl = 0;

    ctrlhist_ncand_perevent->Fill(nCand_ctrl);

    if(nCand_ctrl > 0){
      for(unsigned int iCand(0); iCand < nCand_ctrl; ++iCand){
        //if(nCand_ctrl>1) cout << "cand " << iCand << " isMatched: " << BToKMuMu_isMatched[iCand] << " b pt: " << BToKMuMu_fit_pt[iCand]<< endl;
        if(BToKMuMu_isMatched[iCand] == 1) ++nMatchedCand_ctrl;
      }

      ctrlhist_ncand_matched_perevent->Fill(nMatchedCand_ctrl);

      ctrlhist_ncand_wtriggeringmuon->Fill(ncand_istriggering/ncand_ctrl);

      vector<pair<int,float>> pair_candIdx_desc_bpT_ctrl = createPairWithDesc(nCand_ctrl, BToKMuMu_pt);
      vector<pair<int,float>> pair_candIdx_desc_svprob_ctrl = createPairWithDesc(nCand_ctrl, BToKMuMu_svprob);
      vector<pair<int,float>> pair_candIdx_desc_cos2d_ctrl = createPairWithDesc(nCand_ctrl, BToKMuMu_fit_cos2D);
      vector<pair<int,float>> pair_candIdx_desc_kpt_ctrl = createPairWithDesc(nCand_ctrl, BToKMuMu_fit_k_pt);

      sort(pair_candIdx_desc_bpT_ctrl.begin(), pair_candIdx_desc_bpT_ctrl.end(), sortcansbydesc);
      sort(pair_candIdx_desc_svprob_ctrl.begin(), pair_candIdx_desc_svprob_ctrl.end(), sortcansbydesc);
      sort(pair_candIdx_desc_cos2d_ctrl.begin(), pair_candIdx_desc_cos2d_ctrl.end(), sortcansbydesc);
      sort(pair_candIdx_desc_kpt_ctrl.begin(), pair_candIdx_desc_kpt_ctrl.end(), sortcansbydesc);

      //if(nCand>1){
      //  for(unsigned int iPair(0); iPair < pair_candIdx_desc.size(); ++iPair){
      //    cout << pair_candIdx_desc[iPair].first << " " << pair_candIdx_desc[iPair].second << endl;
      //  }
      //}

      // number of matched selected candidate per event (0 or 1)
      UInt_t selEff_bpt_ctrl = BToKMuMu_isMatched[pair_candIdx_desc_bpT_ctrl[0].first];
      UInt_t selEff_svprob_ctrl = BToKMuMu_isMatched[pair_candIdx_desc_svprob_ctrl[0].first];
      UInt_t selEff_cos2d_ctrl = BToKMuMu_isMatched[pair_candIdx_desc_cos2d_ctrl[0].first];
      UInt_t selEff_kpt_ctrl = BToKMuMu_isMatched[pair_candIdx_desc_kpt_ctrl[0].first];
      //cout << selEff << endl;

      if(nMatchedCand_ctrl != 0 && Muon_isTriggeringBPark[BToKMuMu_l1Idx[pair_candIdx_desc_bpT_ctrl[0].first]]==1){
        ctrlhist_selection_efficiency_bpt_allevents->Fill(selEff_bpt_ctrl);
        ctrlhist_selection_efficiency_svprob_allevents->Fill(selEff_svprob_ctrl);
        ctrlhist_selection_efficiency_cos2d_allevents->Fill(selEff_cos2d_ctrl);
        ctrlhist_selection_efficiency_kpt_allevents->Fill(selEff_kpt_ctrl);

        if(nCand_ctrl > 1) ctrlhist_selection_efficiency_bpt_eventswithmultcands->Fill(selEff_bpt_ctrl);
        if(nCand_ctrl > 1) ctrlhist_selection_efficiency_svprob_eventswithmultcands->Fill(selEff_svprob_ctrl);
        if(nCand_ctrl > 1) ctrlhist_selection_efficiency_cos2d_eventswithmultcands->Fill(selEff_cos2d_ctrl);
        if(nCand_ctrl > 1) ctrlhist_selection_efficiency_kpt_eventswithmultcands->Fill(selEff_kpt_ctrl);
      }
    }
  } // end fill histograms

  return kTRUE;
}


void BToKMuMuDumper::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}


void BToKMuMuDumper::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  my_file->cd();
  if(do_fillhistograms) my_file->Write();

  control_tree->Write("", TObject::kOverwrite);

  my_file->Close();

  cout << "- End B->KMuMu Dumper -" << endl;
}

