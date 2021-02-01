#define NanoDumper_cxx
// The class definition in NanoDumper.h has been generated automatically
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
// root> T->Process("NanoDumper.C")
// root> T->Process("NanoDumper.C","some options")
// root> T->Process("NanoDumper.C+")
//


#include "NanoDumper.h"
#include <TMath.h>
#include <cmath>
#include <TH2.h>
#include <TStyle.h>

using namespace std;

// extra functions
bool sortcansbydesc(const pair<int, float> &a1, const pair<int, float> &a2){
  return a1.second > a2.second;
}

bool sortcansbydesc_opp(const pair<int, float> &a1, const pair<int, float> &a2){
  return a1.second < a2.second;
}

// class functions
vector<pair<int,float>> NanoDumper::createPairWithDesc(const UInt_t& nCand, const TTreeReaderArray<Float_t>& desc){
  vector<pair<int,float>> pair_candIdx_desc;

  for(unsigned int iCand(0); iCand < nCand; ++iCand){
    pair<int, float> pair_candIdx_desc_tmp;
    pair_candIdx_desc_tmp.first  = iCand;
    pair_candIdx_desc_tmp.second = desc[iCand] ;
    pair_candIdx_desc.push_back(pair_candIdx_desc_tmp);
  }

  return pair_candIdx_desc;
}


vector<pair<int,float>> NanoDumper::updatePairWithDesc(vector<pair<int,float>> the_ini_pair, const TTreeReaderArray<Int_t>& charge){
  vector<pair<int,float>> pair_candIdx_desc;

  for(unsigned int iCand(0); iCand < the_ini_pair.size(); ++iCand){
    pair<int, float> pair_candIdx_desc_tmp;
    pair_candIdx_desc_tmp.first = the_ini_pair[iCand].first;
    pair_candIdx_desc_tmp.second = fabs(charge[the_ini_pair[iCand].first]);
    pair_candIdx_desc.push_back(pair_candIdx_desc_tmp);
  }
  return pair_candIdx_desc;
}


void NanoDumper::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  cout << " --------------------------" << endl;
  cout << "        Nano Dumper        " << endl;
  cout << " --------------------------" << endl;
}


void NanoDumper::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  TString outFileName = option;

  // this option is intentionally hardcoded
  do_fillhistograms = false;

  my_file = new TFile(outFileName, "RECREATE");  
  my_file->cd();

  // getting the control tree ready
  control_tree = new TTree("control_tree", "control_tree");

  control_tree->Branch("event", &the_event, "the_event/F");
  control_tree->Branch("run", &the_run, "the_run/I");
  control_tree->Branch("lumi", &the_lumi, "the_lumi/I");

  control_tree->Branch("b_pt", &the_ctrl_b_pt, "the_ctrl_b_pt/F");
  control_tree->Branch("b_eta", &the_ctrl_b_eta, "the_ctrl_b_eta/F");
  control_tree->Branch("b_phi", &the_ctrl_b_phi, "the_ctrl_b_phi/F");
  control_tree->Branch("b_mass", &the_ctrl_b_mass, "the_ctrl_b_mass/F");
  control_tree->Branch("b_charge", &the_ctrl_b_charge, "the_ctrl_b_charge/F");
  control_tree->Branch("b_pdgid", &the_ctrl_b_pdgid, "the_ctrl_b_pdgid/I");
  control_tree->Branch("b_cos2d", &the_ctrl_b_cos2d, "the_ctrl_b_cos2d/F");
  control_tree->Branch("b_iso03", &the_ctrl_b_iso03, "the_ctrl_b_iso03/F");
  control_tree->Branch("b_iso03_close", &the_ctrl_b_iso03_close, "the_ctrl_b_iso03_close/F");
  control_tree->Branch("b_iso04", &the_ctrl_b_iso04, "the_ctrl_b_iso04/F");
  control_tree->Branch("b_iso04_close", &the_ctrl_b_iso04_close, "the_ctrl_b_iso04_close/F");

  control_tree->Branch("k_pt", &the_ctrl_k_pt, "the_ctrl_k_pt/F");
  control_tree->Branch("k_eta", &the_ctrl_k_eta, "the_ctrl_k_eta/F");
  control_tree->Branch("k_phi", &the_ctrl_k_phi, "the_ctrl_k_phi/F");
  control_tree->Branch("k_iso03", &the_ctrl_k_iso03, "the_ctrl_k_iso03/F");
  control_tree->Branch("k_iso03_close", &the_ctrl_k_iso03_close, "the_ctrl_k_iso03_close/F");
  control_tree->Branch("k_iso04", &the_ctrl_k_iso04, "the_ctrl_k_iso04/F");
  control_tree->Branch("k_iso04_close", &the_ctrl_k_iso04_close, "the_ctrl_k_iso04_close/F");

  control_tree->Branch("l1_pt", &the_ctrl_l1_pt, "the_ctrl_l1_pt/F");
  control_tree->Branch("l1_eta", &the_ctrl_l1_eta, "the_ctrl_l1_eta/F");
  control_tree->Branch("l1_phi", &the_ctrl_l1_phi, "the_ctrl_l1_phi/F");
  control_tree->Branch("l1_iso03", &the_ctrl_l1_iso03, "the_ctrl_l1_iso03/F");
  control_tree->Branch("l1_iso03_close", &the_ctrl_l1_iso03_close, "the_ctrl_l1_iso03_close/F");
  control_tree->Branch("l1_iso04", &the_ctrl_l1_iso04, "the_ctrl_l1_iso04/F");
  control_tree->Branch("l1_iso04_close", &the_ctrl_l1_iso04_close, "the_ctrl_l1_iso04_close/F");

  control_tree->Branch("l2_pt", &the_ctrl_l2_pt, "the_ctrl_l2_pt/F");
  control_tree->Branch("l2_eta", &the_ctrl_l2_eta, "the_ctrl_l2_eta/F");
  control_tree->Branch("l2_phi", &the_ctrl_l2_phi, "the_ctrl_l2_phi/F");
  control_tree->Branch("l2_iso03", &the_ctrl_l2_iso03, "the_ctrl_l2_iso03/F");
  control_tree->Branch("l2_iso03_close", &the_ctrl_l2_iso03_close, "the_ctrl_l2_iso03_close/F");
  control_tree->Branch("l2_iso04", &the_ctrl_l2_iso04, "the_ctrl_l2_iso04/F");
  control_tree->Branch("l2_iso04_close", &the_ctrl_l2_iso04_close, "the_ctrl_l2_iso04_close/F");

  control_tree->Branch("dimu_mass", &the_ctrl_dimu_mass, "the_ctrl_dimu_mass/F");
  control_tree->Branch("vtx_x", &the_ctrl_vtx_x, "the_ctrl_vtx_x/F");
  control_tree->Branch("vtx_y", &the_ctrl_vtx_y, "the_ctrl_vtx_y/F");
  control_tree->Branch("vtx_z", &the_ctrl_vtx_z, "the_ctrl_vtx_z/F");
  control_tree->Branch("vtx_lxy", &the_ctrl_vtx_lxy, "the_ctrl_vtx_lxy/F");
  control_tree->Branch("sv_prob", &the_ctrl_sv_prob, "the_ctrl_sv_prob/F");

  control_tree->Branch("pv_npvs", &the_pv_npvs, "the_pv_npvs/I");

  control_tree->Branch("hlt_mu7_ip4", &the_hlt_mu7_ip4, "the_hlt_mu7_ip4/I");
  control_tree->Branch("hlt_mu8_ip6", &the_hlt_mu8_ip6, "the_hlt_mu8_ip6/I");
  control_tree->Branch("hlt_mu8_ip5", &the_hlt_mu8_ip5, "the_hlt_mu8_ip5/I");
  control_tree->Branch("hlt_mu8_ip3", &the_hlt_mu8_ip3, "the_hlt_mu8_ip3/I");
  control_tree->Branch("hlt_mu8p5_ip3p5", &the_hlt_mu8p5_ip3p5, "the_hlt_mu8p5_ip3p5/I");
  control_tree->Branch("hlt_mu9_ip6", &the_hlt_mu9_ip6, "the_hlt_mu9_ip6/I");
  control_tree->Branch("hlt_mu9_ip5", &the_hlt_mu9_ip5, "the_hlt_mu9_ip5/I");
  control_tree->Branch("hlt_mu9_ip4", &the_hlt_mu9_ip4, "the_hlt_mu9_ip4/I");
  control_tree->Branch("hlt_mu10p5_ip3p5", &the_hlt_mu10p5_ip3p5, "the_hlt_mu10p5_ip3p5/I");
  control_tree->Branch("hlt_mu12_ip6", &the_hlt_mu12_ip6, "the_hlt_mu12_ip6/I");


  // getting the signal tree ready
  signal_tree = new TTree("signal_tree", "signal_tree");

  signal_tree->Branch("event", &the_event, "the_event/I");
  signal_tree->Branch("run", &the_run, "the_run/I");
  signal_tree->Branch("lumi", &the_lumi, "the_lumi/I");
  
  signal_tree->Branch("b_pt", &the_sig_b_pt, "the_sig_b_pt/F");
  signal_tree->Branch("b_eta", &the_sig_b_eta, "the_sig_b_eta/F");
  signal_tree->Branch("b_phi", &the_sig_b_phi, "the_sig_b_phi/F");
  signal_tree->Branch("b_mass", &the_sig_b_mass, "the_sig_b_mass/F");
  signal_tree->Branch("b_charge", &the_sig_b_charge, "the_sig_b_charge/F");
  signal_tree->Branch("b_pdgid", &the_sig_b_pdgid, "the_sig_b_pdgid/I");

  signal_tree->Branch("hnl_pt", &the_sig_hnl_pt, "the_sig_hnl_pt/F");
  signal_tree->Branch("hnl_eta", &the_sig_hnl_eta, "the_sig_hnl_eta/F");
  signal_tree->Branch("hnl_phi", &the_sig_hnl_phi, "the_sig_hnl_phi/F");
  signal_tree->Branch("hnl_mass", &the_sig_hnl_mass, "the_sig_hnl_mass/F");
  signal_tree->Branch("hnl_charge", &the_sig_hnl_charge, "the_sig_hnl_charge/F");
  signal_tree->Branch("hnl_cos2d", &the_sig_hnl_cos2d, "the_sig_hnl_cos2d/F");
  signal_tree->Branch("hnl_iso03", &the_sig_hnl_iso03, "the_sig_hnl_iso03/F");
  signal_tree->Branch("hnl_iso03_close", &the_sig_hnl_iso03_close, "the_sig_hnl_iso03_close/F");
  signal_tree->Branch("hnl_iso04", &the_sig_hnl_iso04, "the_sig_hnl_iso04/F");
  signal_tree->Branch("hnl_iso04_close", &the_sig_hnl_iso04_close, "the_sig_hnl_iso04_close/F");

  signal_tree->Branch("trgmu_pt", &the_sig_trgmu_pt, "the_sig_trgmu_pt/F");
  signal_tree->Branch("trgmu_eta", &the_sig_trgmu_eta, "the_sig_trgmu_eta/F");
  signal_tree->Branch("trgmu_phi", &the_sig_trgmu_phi, "the_sig_trgmu_phi/F");
  signal_tree->Branch("trgmu_dxy", &the_sig_trgmu_dxy, "the_sig_trgmu_dxy/F");
  signal_tree->Branch("trgmu_dz", &the_sig_trgmu_dz, "the_sig_trgmu_dz/F");
  signal_tree->Branch("trgmu_ip3d", &the_sig_trgmu_ip3d, "the_sig_trgmu_ip3d/F");
  signal_tree->Branch("trgmu_ip3dsig", &the_sig_trgmu_ip3dsig, "the_sig_trgmu_ip3dsig/F");
  signal_tree->Branch("trgmu_iso03", &the_sig_trgmu_iso03, "the_sig_trgmu_iso03/F");
  signal_tree->Branch("trgmu_iso03_close", &the_sig_trgmu_iso03_close, "the_sig_trgmu_iso03_close/F");
  signal_tree->Branch("trgmu_iso04", &the_sig_trgmu_iso04, "the_sig_trgmu_iso04/F");
  signal_tree->Branch("trgmu_iso04_close", &the_sig_trgmu_iso04_close, "the_sig_trgmu_iso04_close/F");

  signal_tree->Branch("mu_pt", &the_sig_mu_pt, "the_sig_mu_pt/F");
  signal_tree->Branch("mu_eta", &the_sig_mu_eta, "the_sig_mu_eta/F");
  signal_tree->Branch("mu_phi", &the_sig_mu_phi, "the_sig_mu_phi/F");
  signal_tree->Branch("mu_dxy", &the_sig_mu_dxy, "the_sig_mu_dxy/F");
  signal_tree->Branch("mu_dz", &the_sig_mu_dz, "the_sig_mu_dz/F");
  signal_tree->Branch("mu_ip3d", &the_sig_mu_ip3d, "the_sig_mu_ip3d/F");
  signal_tree->Branch("mu_ip3dsig", &the_sig_mu_ip3dsig, "the_sig_mu_ip3dsig/F");
  signal_tree->Branch("mu_iso03", &the_sig_mu_iso03, "the_sig_mu_iso03/F");
  signal_tree->Branch("mu_iso03_close", &the_sig_mu_iso03_close, "the_sig_mu_iso03_close/F");
  signal_tree->Branch("mu_iso04", &the_sig_mu_iso04, "the_sig_mu_iso04/F");
  signal_tree->Branch("mu_iso04_close", &the_sig_mu_iso04_close, "the_sig_mu_iso04_close/F");
  signal_tree->Branch("mu_isloose", &the_sig_mu_isloose, "the_sig_mu_isloose/F");
  signal_tree->Branch("mu_ismedium", &the_sig_mu_ismedium, "the_sig_mu_ismedium/F");
  signal_tree->Branch("mu_istight", &the_sig_mu_istight, "the_sig_mu_istight/F");
  signal_tree->Branch("mu_issoft", &the_sig_mu_issoft, "the_sig_mu_issoft/F");

  signal_tree->Branch("pi_pt", &the_sig_pi_pt, "the_sig_pi_pt/F");
  signal_tree->Branch("pi_eta", &the_sig_pi_eta, "the_sig_pi_eta/F");
  signal_tree->Branch("pi_phi", &the_sig_pi_phi, "the_sig_pi_phi/F");
  signal_tree->Branch("pi_dcasig", &the_sig_pi_dcasig, "the_sig_pi_dcasig/F");
  signal_tree->Branch("pi_dxy", &the_sig_pi_dxy, "the_sig_pi_dxy/F");
  signal_tree->Branch("pi_dz", &the_sig_pi_dz, "the_sig_pi_dz/F");
  signal_tree->Branch("pi_dxysig", &the_sig_pi_dxysig, "the_sig_pi_dxysig/F");
  signal_tree->Branch("pi_dzsig", &the_sig_pi_dzsig, "the_sig_pi_dzsig/F");
  signal_tree->Branch("pi_iso03", &the_sig_pi_iso03, "the_sig_pi_iso03/F");
  signal_tree->Branch("pi_iso03_close", &the_sig_pi_iso03_close, "the_sig_pi_iso03_close/F");
  signal_tree->Branch("pi_iso04", &the_sig_pi_iso04, "the_sig_pi_iso04/F");
  signal_tree->Branch("pi_iso04_close", &the_sig_pi_iso04_close, "the_sig_pi_iso04_close/F");

  signal_tree->Branch("dimu_mass", &the_sig_dimu_mass, "the_sig_dimu_mass/F");
  signal_tree->Branch("dimu_pt", &the_sig_dimu_pt, "the_sig_dimu_pt/F");
  signal_tree->Branch("dimu_lxy", &the_sig_dimu_lxy, "the_sig_dimu_lxy/F");
  signal_tree->Branch("dimu_lxyz", &the_sig_dimu_lxyz, "the_sig_dimu_lxyz/F");
  signal_tree->Branch("dimu_vxdiff", &the_sig_dimu_vxdiff, "the_sig_dimu_vxdiff/F");
  signal_tree->Branch("dimu_vydiff", &the_sig_dimu_vydiff, "the_sig_dimu_vydiff/F");
  signal_tree->Branch("dimu_vzdiff", &the_sig_dimu_vzdiff, "the_sig_dimu_vzdiff/F");

  signal_tree->Branch("deltar_mu_pi", &the_sig_deltar_mu_pi, "the_sig_deltar_mu_pi/F");
  signal_tree->Branch("deltar_trgmu_hnl", &the_sig_deltar_trgmu_hnl, "the_sig_deltar_trgmu_hnl/F");
  signal_tree->Branch("sv_chi2", &the_sig_sv_chi2, "the_sig_sv_chi2/F");
  signal_tree->Branch("sv_lxy", &the_sig_sv_lxy, "the_sig_sv_lxy/F");
  signal_tree->Branch("sv_lxysig", &the_sig_sv_lxysig, "the_sig_sv_lxysig/F");
  signal_tree->Branch("sv_prob", &the_sig_sv_prob, "the_sig_sv_prob/F");
  signal_tree->Branch("sv_x", &the_sig_sv_x, "the_sig_sv_x/F");
  signal_tree->Branch("sv_y", &the_sig_sv_y, "the_sig_sv_y/F");
  signal_tree->Branch("sv_z", &the_sig_sv_z, "the_sig_sv_z/F");

  signal_tree->Branch("pi_mu_vzdiff", &the_sig_pi_mu_vzdiff, "the_sig_pi_mu_vzdiff/F");

  signal_tree->Branch("pv_npvs", &the_pv_npvs, "the_pv_npvs/I");

  signal_tree->Branch("hlt_mu7_ip4", &the_hlt_mu7_ip4, "the_hlt_mu7_ip4/C");
  signal_tree->Branch("hlt_mu8_ip6", &the_hlt_mu8_ip6, "the_hlt_mu8_ip6/C");
  signal_tree->Branch("hlt_mu8_ip5", &the_hlt_mu8_ip5, "the_hlt_mu8_ip5/C");
  signal_tree->Branch("hlt_mu8_ip3", &the_hlt_mu8_ip3, "the_hlt_mu8_ip3/C");
  signal_tree->Branch("hlt_mu8p5_ip3p5", &the_hlt_mu8p5_ip3p5, "the_hlt_mu8p5_ip3p5/C");
  signal_tree->Branch("hlt_mu9_ip6", &the_hlt_mu9_ip6, "the_hlt_mu9_ip6/C");
  signal_tree->Branch("hlt_mu9_ip5", &the_hlt_mu9_ip5, "the_hlt_mu9_ip5/C");
  signal_tree->Branch("hlt_mu9_ip4", &the_hlt_mu9_ip4, "the_hlt_mu9_ip4/C");
  signal_tree->Branch("hlt_mu10p5_ip3p5", &the_hlt_mu10p5_ip3p5, "the_hlt_mu10p5_ip3p5/C");
  signal_tree->Branch("hlt_mu12_ip6", &the_hlt_mu12_ip6, "the_hlt_mu12_ip6/C");


  // defining histograms
  if(do_fillhistograms){
    my_file->mkdir("control_channel");
    my_file->mkdir("signal_channel");
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

    my_file->cd("signal_channel");

    sighist_ncand_perevent = new TH1F("sighist_ncand_perevent", "sighist_ncand_perevent", 10, 0, 10);
    sighist_ncand_matched_perevent = new TH1F("sighist_ncand_matched_perevent", "sighist_ncand_matched__perevent", 10, 0, 10);
    sighist_selection_efficiency_hnlpt_allevents = new TH1F("sighist_selection_efficiency_hnlpt_allevents", "Efficiency of selection of candidate with largest hnl pT (all events)", 2, 0, 2);
    sighist_selection_efficiency_hnlpt_eventswithmultcands = new TH1F("sighist_selection_efficiency_hnlpt_eventswithmultcands", "Efficiency of selection of candidate with largest hnl pT (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_bpt_allevents = new TH1F("sighist_selection_efficiency_bpt_allevents", "Efficiency of selection of candidate with largest b pT (all events)", 2, 0, 2);
    sighist_selection_efficiency_bpt_eventswithmultcands = new TH1F("sighist_selection_efficiency_bpt_eventswithmultcands", "Efficiency of selection of candidate with largest b pT (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_trgmupt_allevents = new TH1F("sighist_selection_efficiency_trgmupt_allevents", "Efficiency of selection of candidate with largest trigger muon pT (all events)", 2, 0, 2);
    sighist_selection_efficiency_trgmupt_eventswithmultcands = new TH1F("sighist_selection_efficiency_trgmupt_eventswithmultcands", "Efficiency of selection of candidate with largest trigger muon pT (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_pipt_allevents = new TH1F("sighist_selection_efficiency_pipt_allevents", "Efficiency of selection of candidate with largest pi pT (all events)", 2, 0, 2);
    sighist_selection_efficiency_pipt_eventswithmultcands = new TH1F("sighist_selection_efficiency_pipt_eventswithmultcands", "Efficiency of selection of candidate with largest pi pT (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_svprob_allevents = new TH1F("sighist_selection_efficiency_svprob_allevents", "Efficiency of selection of candidate with largest vertex prob. (all events)", 2, 0, 2);
    sighist_selection_efficiency_svprob_eventswithmultcands = new TH1F("sighist_selection_efficiency_svprob_eventswithmultcands", "Efficiency of selection of candidate with largest vertex prob. (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_svchi2_allevents = new TH1F("sighist_selection_efficiency_svchi2_allevents", "Efficiency of selection of candidate with largest vertex chi2 (all events)", 2, 0, 2);
    sighist_selection_efficiency_svchi2_eventswithmultcands = new TH1F("sighist_selection_efficiency_svchi2_eventswithmultcands", "Efficiency of selection of candidate with largest vertex chi2 (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_cos2d_allevents = new TH1F("sighist_selection_efficiency_cos2d_allevents", "Efficiency of selection of candidate with largest vertex cos(angle) (all events)", 2, 0, 2);
    sighist_selection_efficiency_cos2d_eventswithmultcands = new TH1F("sighist_selection_efficiency_cos2d_eventswithmultcands", "Efficiency of selection of candidate with largest vertex cos(angle) (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_hnliso4_allevents = new TH1F("sighist_selection_efficiency_hnliso4_allevents", "Efficiency of selection of candidate with smallest hnl isolation (all events)", 2, 0, 2);
    sighist_selection_efficiency_hnliso4_eventswithmultcands = new TH1F("sighist_selection_efficiency_hnliso4_eventswithmultcands", "Efficiency of selection of candidate with smallest hnl isolation (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_dr_allevents = new TH1F("sighist_selection_efficiency_dr_allevents", "Efficiency of selection of candidate with smallest dR(trgmu, hnl) (all events)", 2, 0, 2);
    sighist_selection_efficiency_dr_eventswithmultcands = new TH1F("sighist_selection_efficiency_dr_eventswithmultcands", "Efficiency of selection of candidate with smallest dR(trgmu, hnl) (events with multiple candidates)", 2, 0, 2);
  } // end define histograms
}


Bool_t NanoDumper::Process(Long64_t entry)
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

  fReader.SetEntry(entry);
  //cout << endl << "--- Entry " << entry << " ---" << endl;

  // number of candidates in the event
  UInt_t nCand_ctrl = *nBToKMuMu; 
  UInt_t nCand_sig = *nBToMuMuPi; 

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
  
  if(nCand_ctrl > 0){ // at least one candidate per event

    // selecting the candidate as the one having the largest b pt
    // - create candIdx - b pt pair
    vector<pair<int,float>> pair_candIdx_desc_bpt_ctrl = createPairWithDesc(nCand_ctrl, BToKMuMu_fit_pt);
    // - sort pair with decreasing b pt
    sort(pair_candIdx_desc_bpt_ctrl.begin(), pair_candIdx_desc_bpt_ctrl.end(), sortcansbydesc);
    // - fetch candIdx associated to the largest b pt
    UInt_t selectedCandIdx_ctrl = pair_candIdx_desc_bpt_ctrl[0].first;

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
    the_ctrl_k_iso03 = BToKMuMu_k_iso03[selectedCandIdx_ctrl];
    the_ctrl_k_iso03_close = BToKMuMu_k_iso03_close[selectedCandIdx_ctrl];
    the_ctrl_k_iso04 = BToKMuMu_k_iso04[selectedCandIdx_ctrl];
    the_ctrl_k_iso04_close = BToKMuMu_k_iso04_close[selectedCandIdx_ctrl];

    the_ctrl_l1_pt = BToKMuMu_fit_l1_pt[selectedCandIdx_ctrl];
    the_ctrl_l1_eta = BToKMuMu_fit_l1_eta[selectedCandIdx_ctrl];
    the_ctrl_l1_phi = BToKMuMu_fit_l1_phi[selectedCandIdx_ctrl];
    the_ctrl_l1_iso03 = BToKMuMu_l1_iso03[selectedCandIdx_ctrl];
    the_ctrl_l1_iso03_close = BToKMuMu_l1_iso03_close[selectedCandIdx_ctrl];
    the_ctrl_l1_iso04 = BToKMuMu_l1_iso04[selectedCandIdx_ctrl];
    the_ctrl_l1_iso04_close = BToKMuMu_l1_iso04_close[selectedCandIdx_ctrl];

    the_ctrl_l2_pt = BToKMuMu_fit_l2_pt[selectedCandIdx_ctrl];
    the_ctrl_l2_eta = BToKMuMu_fit_l2_eta[selectedCandIdx_ctrl];
    the_ctrl_l2_phi = BToKMuMu_fit_l2_phi[selectedCandIdx_ctrl];
    the_ctrl_l2_iso03 = BToKMuMu_l2_iso03[selectedCandIdx_ctrl];
    the_ctrl_l2_iso03_close = BToKMuMu_l2_iso03_close[selectedCandIdx_ctrl];
    the_ctrl_l2_iso04 = BToKMuMu_l2_iso04[selectedCandIdx_ctrl];
    the_ctrl_l2_iso04_close = BToKMuMu_l2_iso04_close[selectedCandIdx_ctrl];

    the_ctrl_dimu_mass = BToKMuMu_mll_fullfit[selectedCandIdx_ctrl];
    the_ctrl_vtx_x = BToKMuMu_vtx_x[selectedCandIdx_ctrl];
    the_ctrl_vtx_y = BToKMuMu_vtx_x[selectedCandIdx_ctrl];
    the_ctrl_vtx_z = BToKMuMu_vtx_x[selectedCandIdx_ctrl];
    the_ctrl_vtx_lxy = BToKMuMu_vtx_x[selectedCandIdx_ctrl];
    the_ctrl_sv_prob = BToKMuMu_svprob[selectedCandIdx_ctrl];

    control_tree->Fill();

  }// end at least one candidate in the event


  //   ----- Signal Channel -----  //

  if(nCand_sig > 0){ // at least one candidate per event

    // selecting the candidate as the one having the largest hnl pt
    // - create candIdx - hnl pt pairs
    vector<pair<int,float>> pair_candIdx_desc_hnlpt_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_pt);
    // - sort it in decreasing hnl pt
    sort(pair_candIdx_desc_hnlpt_sig.begin(), pair_candIdx_desc_hnlpt_sig.end(), sortcansbydesc);
    // - then privilege OS cand over SS ones
    vector<pair<int,float>> pair_candIdx_desc_hnlpt_sig_up = updatePairWithDesc(pair_candIdx_desc_hnlpt_sig, BToMuMuPi_hnl_charge);
    sort(pair_candIdx_desc_hnlpt_sig_up.begin(), pair_candIdx_desc_hnlpt_sig_up.end(), sortcansbydesc_opp);
    // - and select the OS cand with the largest hnl pt
    UInt_t selectedCandIdx_sig = pair_candIdx_desc_hnlpt_sig_up[0].first;

    // fill the signal_tree
    the_sig_b_pt = BToMuMuPi_pt[selectedCandIdx_sig];
    the_sig_b_eta = BToMuMuPi_eta[selectedCandIdx_sig];
    the_sig_b_phi = BToMuMuPi_phi[selectedCandIdx_sig];
    the_sig_b_mass = BToMuMuPi_mass[selectedCandIdx_sig];
    the_sig_b_charge = BToMuMuPi_charge[selectedCandIdx_sig];
    the_sig_b_pdgid = BToMuMuPi_pdgId[selectedCandIdx_sig];

    the_sig_hnl_pt = BToMuMuPi_hnl_pt[selectedCandIdx_sig];
    the_sig_hnl_eta = BToMuMuPi_hnl_eta[selectedCandIdx_sig];
    the_sig_hnl_phi = BToMuMuPi_hnl_phi[selectedCandIdx_sig];
    the_sig_hnl_mass = BToMuMuPi_hnl_mass[selectedCandIdx_sig];
    the_sig_hnl_charge = BToMuMuPi_hnl_charge[selectedCandIdx_sig];
    the_sig_hnl_cos2d = BToMuMuPi_hnl_cos2D[selectedCandIdx_sig];
    the_sig_hnl_iso03 = BToMuMuPi_hnl_iso03[selectedCandIdx_sig];
    the_sig_hnl_iso03_close = BToMuMuPi_hnl_iso03_close[selectedCandIdx_sig];
    the_sig_hnl_iso04 = BToMuMuPi_hnl_iso04[selectedCandIdx_sig];
    the_sig_hnl_iso04_close = BToMuMuPi_hnl_iso04_close[selectedCandIdx_sig];

    the_sig_trgmu_pt = BToMuMuPi_trg_mu_pt[selectedCandIdx_sig];
    the_sig_trgmu_eta = BToMuMuPi_trg_mu_eta[selectedCandIdx_sig];
    the_sig_trgmu_phi = BToMuMuPi_trg_mu_phi[selectedCandIdx_sig];
    the_sig_trgmu_dxy = BToMuMuPi_trg_mu_dxy[selectedCandIdx_sig];
    the_sig_trgmu_dz = BToMuMuPi_trg_mu_dz[selectedCandIdx_sig];
    the_sig_trgmu_ip3d = BToMuMuPi_trg_mu_ip3d[selectedCandIdx_sig];
    the_sig_trgmu_ip3dsig = BToMuMuPi_trg_mu_sip3d[selectedCandIdx_sig];
    the_sig_trgmu_iso03 = BToMuMuPi_trg_mu_iso03[selectedCandIdx_sig];
    the_sig_trgmu_iso03_close = BToMuMuPi_trg_mu_iso03_close[selectedCandIdx_sig];
    the_sig_trgmu_iso04 = BToMuMuPi_trg_mu_iso04[selectedCandIdx_sig];
    the_sig_trgmu_iso04_close = BToMuMuPi_trg_mu_iso04_close[selectedCandIdx_sig];

    the_sig_mu_pt = BToMuMuPi_fit_mu_pt[selectedCandIdx_sig];
    the_sig_mu_eta = BToMuMuPi_fit_mu_eta[selectedCandIdx_sig];
    the_sig_mu_phi = BToMuMuPi_fit_mu_phi[selectedCandIdx_sig];
    the_sig_mu_dxy = BToMuMuPi_sel_mu_dxy[selectedCandIdx_sig];
    the_sig_mu_dz = BToMuMuPi_sel_mu_dz[selectedCandIdx_sig];
    the_sig_mu_ip3d = BToMuMuPi_sel_mu_ip3d[selectedCandIdx_sig];
    the_sig_mu_ip3dsig = BToMuMuPi_sel_mu_sip3d[selectedCandIdx_sig];
    the_sig_mu_iso03 = BToMuMuPi_sel_mu_iso03[selectedCandIdx_sig];
    the_sig_mu_iso03_close = BToMuMuPi_sel_mu_iso03_close[selectedCandIdx_sig];
    the_sig_mu_iso04 = BToMuMuPi_sel_mu_iso04[selectedCandIdx_sig];
    the_sig_mu_iso04_close = BToMuMuPi_sel_mu_iso04_close[selectedCandIdx_sig];
    the_sig_mu_isloose = BToMuMuPi_sel_mu_isLoose[selectedCandIdx_sig];
    the_sig_mu_ismedium = BToMuMuPi_sel_mu_isMedium[selectedCandIdx_sig];
    the_sig_mu_istight = BToMuMuPi_sel_mu_isTight[selectedCandIdx_sig];
    the_sig_mu_issoft = BToMuMuPi_sel_mu_isSoft[selectedCandIdx_sig];

    the_sig_pi_pt = BToMuMuPi_fit_pi_pt[selectedCandIdx_sig];
    the_sig_pi_eta = BToMuMuPi_fit_pi_eta[selectedCandIdx_sig];
    the_sig_pi_phi = BToMuMuPi_fit_pi_phi[selectedCandIdx_sig];
    the_sig_pi_dcasig = BToMuMuPi_pi_DCASig[selectedCandIdx_sig];
    the_sig_pi_dxy = BToMuMuPi_pi_dxy[selectedCandIdx_sig];
    the_sig_pi_dz = BToMuMuPi_pi_dz[selectedCandIdx_sig];
    the_sig_pi_dxysig = BToMuMuPi_pi_dxyS[selectedCandIdx_sig];
    the_sig_pi_dzsig = BToMuMuPi_pi_dzS[selectedCandIdx_sig];
    the_sig_pi_iso03 = BToMuMuPi_pi_iso03[selectedCandIdx_sig];
    the_sig_pi_iso03_close = BToMuMuPi_pi_iso03_close[selectedCandIdx_sig];
    the_sig_pi_iso04 = BToMuMuPi_pi_iso04[selectedCandIdx_sig];
    the_sig_pi_iso04_close = BToMuMuPi_pi_iso04_close[selectedCandIdx_sig];

    the_sig_dimu_mass = BToMuMuPi_dilepton_mass[selectedCandIdx_sig];
    the_sig_dimu_pt = BToMuMuPi_dilepton_pt[selectedCandIdx_sig];
    the_sig_dimu_lxy = BToMuMuPi_dimu_Lxy[selectedCandIdx_sig];
    the_sig_dimu_lxyz = BToMuMuPi_dimu_Lxyz[selectedCandIdx_sig];
    the_sig_dimu_vxdiff = BToMuMuPi_dimu_vxdiff[selectedCandIdx_sig];
    the_sig_dimu_vydiff = BToMuMuPi_dimu_vydiff[selectedCandIdx_sig];
    the_sig_dimu_vzdiff = BToMuMuPi_dimu_vzdiff[selectedCandIdx_sig];

    the_sig_deltar_mu_pi = BToMuMuPi_dr_mu_pi[selectedCandIdx_sig];
    the_sig_deltar_trgmu_hnl = BToMuMuPi_dr_trgmu_hnl[selectedCandIdx_sig];

    the_sig_sv_chi2 = BToMuMuPi_sv_chi2[selectedCandIdx_sig];
    the_sig_sv_lxy = BToMuMuPi_sv_lxy[selectedCandIdx_sig];
    the_sig_sv_lxysig = BToMuMuPi_sv_lxy_sig[selectedCandIdx_sig];
    the_sig_sv_prob = BToMuMuPi_sv_prob[selectedCandIdx_sig];
    the_sig_sv_x = BToMuMuPi_sv_x[selectedCandIdx_sig];
    the_sig_sv_y = BToMuMuPi_sv_y[selectedCandIdx_sig];
    the_sig_sv_z = BToMuMuPi_sv_z[selectedCandIdx_sig];

    the_sig_pi_mu_vzdiff = BToMuMuPi_pi_mu_vzdiff[selectedCandIdx_sig];

    signal_tree->Fill();

  }// end at least one candidate in the event

  
  //   ----- Histograms -----  //
  
  if(do_fillhistograms){
    
    // number of matched candidates in the event (only saved if nCand non zero)
    UInt_t nMatchedCand_sig = 0;
    UInt_t nMatchedCand_ctrl = 0;
    
    sighist_ncand_perevent->Fill(nCand_sig);
    ctrlhist_ncand_perevent->Fill(nCand_ctrl);

    if(nCand_ctrl > 0){
      for(unsigned int iCand(0); iCand < nCand_ctrl; ++iCand){
        //if(nCand_ctrl>1) cout << "cand " << iCand << " isMatched: " << BToKMuMu_isMatched[iCand] << " b pt: " << BToKMuMu_fit_pt[iCand]<< endl;
        if(BToKMuMu_isMatched[iCand] == 1) ++nMatchedCand_ctrl;
      }

      ctrlhist_ncand_matched_perevent->Fill(nMatchedCand_ctrl);

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

      if(nMatchedCand_ctrl != 0){
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

    if(nCand_sig > 0){
      for(unsigned int iCand(0); iCand < nCand_sig; ++iCand){
        //cout << "cand " << iCand << " isMatched: " << BToMuMuPi_isMatched[iCand] << " b pt: " <<  BToMuMuPi_pt[iCand] << " hnl charge " << BToMuMuPi_hnl_charge[iCand] << endl;
        if(BToMuMuPi_isMatched[iCand] == 1) ++nMatchedCand_sig;
      }

      sighist_ncand_matched_perevent->Fill(nMatchedCand_sig);

      vector<pair<int,float>> pair_candIdx_desc_hnlpT_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_pt);
      vector<pair<int,float>> pair_candIdx_desc_bpt_sig = createPairWithDesc(nCand_sig, BToMuMuPi_pt);
      vector<pair<int,float>> pair_candIdx_desc_trgmupt_sig = createPairWithDesc(nCand_sig, BToMuMuPi_trg_mu_pt);
      vector<pair<int,float>> pair_candIdx_desc_pipt_sig = createPairWithDesc(nCand_sig, BToMuMuPi_fit_pi_pt);
      vector<pair<int,float>> pair_candIdx_desc_svprob_sig = createPairWithDesc(nCand_sig, BToMuMuPi_sv_prob);
      vector<pair<int,float>> pair_candIdx_desc_svchi2_sig = createPairWithDesc(nCand_sig, BToMuMuPi_sv_chi2);
      vector<pair<int,float>> pair_candIdx_desc_cos2d_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_cos2D);
      vector<pair<int,float>> pair_candIdx_desc_dr_sig = createPairWithDesc(nCand_sig, BToMuMuPi_dr_trgmu_hnl);
      vector<pair<int,float>> pair_candIdx_desc_hnliso4_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_iso04_close);

      sort(pair_candIdx_desc_hnlpT_sig.begin(), pair_candIdx_desc_hnlpT_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_bpt_sig.begin(), pair_candIdx_desc_bpt_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_trgmupt_sig.begin(), pair_candIdx_desc_trgmupt_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_pipt_sig.begin(), pair_candIdx_desc_pipt_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_svprob_sig.begin(), pair_candIdx_desc_svprob_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_svchi2_sig.begin(), pair_candIdx_desc_svchi2_sig.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_cos2d_sig.begin(), pair_candIdx_desc_cos2d_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_hnliso4_sig.begin(), pair_candIdx_desc_hnliso4_sig.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_dr_sig.begin(), pair_candIdx_desc_dr_sig.end(), sortcansbydesc_opp);
    
      // then privilege OS cand over SS ones
      vector<pair<int,float>> pair_candIdx_desc_hnlpT_sig_up = updatePairWithDesc(pair_candIdx_desc_hnlpT_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_bpt_sig_up = updatePairWithDesc(pair_candIdx_desc_bpt_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_trgmupt_sig_up = updatePairWithDesc(pair_candIdx_desc_trgmupt_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_pipt_sig_up = updatePairWithDesc(pair_candIdx_desc_pipt_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_svprob_sig_up = updatePairWithDesc(pair_candIdx_desc_svprob_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_svchi2_sig_up = updatePairWithDesc(pair_candIdx_desc_svchi2_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_cos2d_sig_up = updatePairWithDesc(pair_candIdx_desc_cos2d_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_hnliso4_sig_up = updatePairWithDesc(pair_candIdx_desc_hnliso4_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_dr_sig_up = updatePairWithDesc(pair_candIdx_desc_dr_sig, BToMuMuPi_hnl_charge);

      sort(pair_candIdx_desc_hnlpT_sig_up.begin(), pair_candIdx_desc_hnlpT_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_bpt_sig_up.begin(), pair_candIdx_desc_bpt_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_trgmupt_sig_up.begin(), pair_candIdx_desc_trgmupt_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_pipt_sig_up.begin(), pair_candIdx_desc_pipt_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_svprob_sig_up.begin(), pair_candIdx_desc_svprob_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_svchi2_sig_up.begin(), pair_candIdx_desc_svchi2_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_cos2d_sig_up.begin(), pair_candIdx_desc_cos2d_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_hnliso4_sig_up.begin(), pair_candIdx_desc_hnliso4_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_dr_sig_up.begin(), pair_candIdx_desc_dr_sig_up.end(), sortcansbydesc_opp);

      //if(nCand_sig>1){
      //for(unsigned int iPair(0); iPair < pair_candIdx_desc_bpt_sig.size(); ++iPair){
      //  cout << "idx: " << pair_candIdx_desc_bpt_sig[iPair].first << " b pt: " << pair_candIdx_desc_bpt_sig[iPair].second << endl;
      //}
      //for(unsigned int iPair(0); iPair < pair_candIdx_desc_bpt_sig_up.size(); ++iPair){
      //  cout << "idx: " << pair_candIdx_desc_bpt_sig_up[iPair].first << " charge: " << pair_candIdx_desc_bpt_sig_up[iPair].second << endl;
      //}
      //}
      //cout << "selected cand idx: " << selectedCandIdx_sig << " is matched: " << BToMuMuPi_isMatched[selectedCandIdx_sig] << endl;

      // number of matched selected candidate per event (0 or 1)
      UInt_t selEff_hnlpt_sig = BToMuMuPi_isMatched[pair_candIdx_desc_hnlpT_sig_up[0].first];
      UInt_t selEff_bpt_sig = BToMuMuPi_isMatched[pair_candIdx_desc_bpt_sig_up[0].first];
      UInt_t selEff_trgmupt_sig = BToMuMuPi_isMatched[pair_candIdx_desc_trgmupt_sig_up[0].first];
      UInt_t selEff_pipt_sig = BToMuMuPi_isMatched[pair_candIdx_desc_pipt_sig_up[0].first];
      UInt_t selEff_svprob_sig = BToMuMuPi_isMatched[pair_candIdx_desc_svprob_sig_up[0].first];
      UInt_t selEff_svchi2_sig = BToMuMuPi_isMatched[pair_candIdx_desc_svchi2_sig_up[0].first];
      UInt_t selEff_cos2d_sig = BToMuMuPi_isMatched[pair_candIdx_desc_cos2d_sig_up[0].first];
      UInt_t selEff_hnliso4_sig = BToMuMuPi_isMatched[pair_candIdx_desc_hnliso4_sig_up[0].first];
      UInt_t selEff_dr_sig = BToMuMuPi_isMatched[pair_candIdx_desc_dr_sig_up[0].first];

      if(nMatchedCand_sig != 0){
        sighist_selection_efficiency_hnlpt_allevents->Fill(selEff_hnlpt_sig);
        sighist_selection_efficiency_bpt_allevents->Fill(selEff_bpt_sig);
        sighist_selection_efficiency_trgmupt_allevents->Fill(selEff_trgmupt_sig);
        sighist_selection_efficiency_pipt_allevents->Fill(selEff_pipt_sig);
        sighist_selection_efficiency_svprob_allevents->Fill(selEff_svprob_sig);
        sighist_selection_efficiency_svchi2_allevents->Fill(selEff_svchi2_sig);
        sighist_selection_efficiency_cos2d_allevents->Fill(selEff_cos2d_sig);
        sighist_selection_efficiency_hnliso4_allevents->Fill(selEff_hnliso4_sig);
        sighist_selection_efficiency_dr_allevents->Fill(selEff_dr_sig);

        if(nCand_sig > 1) sighist_selection_efficiency_hnlpt_eventswithmultcands->Fill(selEff_hnlpt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_bpt_eventswithmultcands->Fill(selEff_bpt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_trgmupt_eventswithmultcands->Fill(selEff_trgmupt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_pipt_eventswithmultcands->Fill(selEff_pipt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_svprob_eventswithmultcands->Fill(selEff_svprob_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_svchi2_eventswithmultcands->Fill(selEff_svchi2_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_cos2d_eventswithmultcands->Fill(selEff_cos2d_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_hnliso4_eventswithmultcands->Fill(selEff_hnliso4_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_dr_eventswithmultcands->Fill(selEff_dr_sig);
      }
    } // end at least one candidate
  } // end fill histograms
   
  return kTRUE;
}


void NanoDumper::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}


void NanoDumper::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  my_file->Write();
  my_file->Close();

  cout << "- End Nano Dumper -" << endl;
}
