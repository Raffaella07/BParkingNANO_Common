//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 23 15:15:57 2021 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: bparknano_nj89.root
//////////////////////////////////////////////////////////

#ifndef NanoDumper_h
#define NanoDumper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class NanoDumper : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> nBToKMuMu = {fReader, "nBToKMuMu"};
   TTreeReaderArray<Float_t> BToKMuMu_b_iso03 = {fReader, "BToKMuMu_b_iso03"};
   TTreeReaderArray<Float_t> BToKMuMu_b_iso03_close = {fReader, "BToKMuMu_b_iso03_close"};
   TTreeReaderArray<Float_t> BToKMuMu_b_iso04 = {fReader, "BToKMuMu_b_iso04"};
   TTreeReaderArray<Float_t> BToKMuMu_b_iso04_close = {fReader, "BToKMuMu_b_iso04_close"};
   TTreeReaderArray<Float_t> BToKMuMu_cos2D = {fReader, "BToKMuMu_cos2D"};
   TTreeReaderArray<Float_t> BToKMuMu_eta = {fReader, "BToKMuMu_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_cos2D = {fReader, "BToKMuMu_fit_cos2D"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_eta = {fReader, "BToKMuMu_fit_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_k_eta = {fReader, "BToKMuMu_fit_k_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_k_phi = {fReader, "BToKMuMu_fit_k_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_k_pt = {fReader, "BToKMuMu_fit_k_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l1_eta = {fReader, "BToKMuMu_fit_l1_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l1_phi = {fReader, "BToKMuMu_fit_l1_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l1_pt = {fReader, "BToKMuMu_fit_l1_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l2_eta = {fReader, "BToKMuMu_fit_l2_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l2_phi = {fReader, "BToKMuMu_fit_l2_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l2_pt = {fReader, "BToKMuMu_fit_l2_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_mass = {fReader, "BToKMuMu_fit_mass"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_massErr = {fReader, "BToKMuMu_fit_massErr"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_phi = {fReader, "BToKMuMu_fit_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_pt = {fReader, "BToKMuMu_fit_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_k_iso03 = {fReader, "BToKMuMu_k_iso03"};
   TTreeReaderArray<Float_t> BToKMuMu_k_iso03_close = {fReader, "BToKMuMu_k_iso03_close"};
   TTreeReaderArray<Float_t> BToKMuMu_k_iso04 = {fReader, "BToKMuMu_k_iso04"};
   TTreeReaderArray<Float_t> BToKMuMu_k_iso04_close = {fReader, "BToKMuMu_k_iso04_close"};
   TTreeReaderArray<Float_t> BToKMuMu_l1_iso03 = {fReader, "BToKMuMu_l1_iso03"};
   TTreeReaderArray<Float_t> BToKMuMu_l1_iso03_close = {fReader, "BToKMuMu_l1_iso03_close"};
   TTreeReaderArray<Float_t> BToKMuMu_l1_iso04 = {fReader, "BToKMuMu_l1_iso04"};
   TTreeReaderArray<Float_t> BToKMuMu_l1_iso04_close = {fReader, "BToKMuMu_l1_iso04_close"};
   TTreeReaderArray<Float_t> BToKMuMu_l2_iso03 = {fReader, "BToKMuMu_l2_iso03"};
   TTreeReaderArray<Float_t> BToKMuMu_l2_iso03_close = {fReader, "BToKMuMu_l2_iso03_close"};
   TTreeReaderArray<Float_t> BToKMuMu_l2_iso04 = {fReader, "BToKMuMu_l2_iso04"};
   TTreeReaderArray<Float_t> BToKMuMu_l2_iso04_close = {fReader, "BToKMuMu_l2_iso04_close"};
   TTreeReaderArray<Float_t> BToKMuMu_l_xy = {fReader, "BToKMuMu_l_xy"};
   TTreeReaderArray<Float_t> BToKMuMu_l_xy_unc = {fReader, "BToKMuMu_l_xy_unc"};
   TTreeReaderArray<Float_t> BToKMuMu_mass = {fReader, "BToKMuMu_mass"};
   TTreeReaderArray<Float_t> BToKMuMu_maxDR = {fReader, "BToKMuMu_maxDR"};
   TTreeReaderArray<Float_t> BToKMuMu_minDR = {fReader, "BToKMuMu_minDR"};
   TTreeReaderArray<Float_t> BToKMuMu_mllErr_llfit = {fReader, "BToKMuMu_mllErr_llfit"};
   TTreeReaderArray<Float_t> BToKMuMu_mll_fullfit = {fReader, "BToKMuMu_mll_fullfit"};
   TTreeReaderArray<Float_t> BToKMuMu_mll_llfit = {fReader, "BToKMuMu_mll_llfit"};
   TTreeReaderArray<Float_t> BToKMuMu_mll_raw = {fReader, "BToKMuMu_mll_raw"};
   TTreeReaderArray<Float_t> BToKMuMu_phi = {fReader, "BToKMuMu_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_pt = {fReader, "BToKMuMu_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_svprob = {fReader, "BToKMuMu_svprob"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_ex = {fReader, "BToKMuMu_vtx_ex"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_ey = {fReader, "BToKMuMu_vtx_ey"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_ez = {fReader, "BToKMuMu_vtx_ez"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_x = {fReader, "BToKMuMu_vtx_x"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_y = {fReader, "BToKMuMu_vtx_y"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_z = {fReader, "BToKMuMu_vtx_z"};
   TTreeReaderArray<Int_t> BToKMuMu_charge = {fReader, "BToKMuMu_charge"};
   TTreeReaderArray<Int_t> BToKMuMu_isMatched = {fReader, "BToKMuMu_isMatched"};
   TTreeReaderArray<Int_t> BToKMuMu_kIdx = {fReader, "BToKMuMu_kIdx"};
   TTreeReaderArray<Int_t> BToKMuMu_l1Idx = {fReader, "BToKMuMu_l1Idx"};
   TTreeReaderArray<Int_t> BToKMuMu_l2Idx = {fReader, "BToKMuMu_l2Idx"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_k_genIdx = {fReader, "BToKMuMu_matching_k_genIdx"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_k_motherPdgId = {fReader, "BToKMuMu_matching_k_motherPdgId"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_l1_genIdx = {fReader, "BToKMuMu_matching_l1_genIdx"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_l1_motherPdgId = {fReader, "BToKMuMu_matching_l1_motherPdgId"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_l2_genIdx = {fReader, "BToKMuMu_matching_l2_genIdx"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_l2_motherPdgId = {fReader, "BToKMuMu_matching_l2_motherPdgId"};
   TTreeReaderArray<Int_t> BToKMuMu_n_k_used = {fReader, "BToKMuMu_n_k_used"};
   TTreeReaderArray<Int_t> BToKMuMu_n_l1_used = {fReader, "BToKMuMu_n_l1_used"};
   TTreeReaderArray<Int_t> BToKMuMu_n_l2_used = {fReader, "BToKMuMu_n_l2_used"};
   TTreeReaderArray<Int_t> BToKMuMu_pdgId = {fReader, "BToKMuMu_pdgId"};
   TTreeReaderValue<UInt_t> nBToMuMuPi = {fReader, "nBToMuMuPi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dilepton_mass = {fReader, "BToMuMuPi_dilepton_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_dilepton_pt = {fReader, "BToMuMuPi_dilepton_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_dimu_Lxy = {fReader, "BToMuMuPi_dimu_Lxy"};
   TTreeReaderArray<Float_t> BToMuMuPi_dimu_Lxyz = {fReader, "BToMuMuPi_dimu_Lxyz"};
   TTreeReaderArray<Float_t> BToMuMuPi_dimu_vxdiff = {fReader, "BToMuMuPi_dimu_vxdiff"};
   TTreeReaderArray<Float_t> BToMuMuPi_dimu_vydiff = {fReader, "BToMuMuPi_dimu_vydiff"};
   TTreeReaderArray<Float_t> BToMuMuPi_dimu_vzdiff = {fReader, "BToMuMuPi_dimu_vzdiff"};
   TTreeReaderArray<Float_t> BToMuMuPi_dr_mu_pi = {fReader, "BToMuMuPi_dr_mu_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dr_trgmu_hnl = {fReader, "BToMuMuPi_dr_trgmu_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_eta = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_mu_eta = {fReader, "BToMuMuPi_fit_mu_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_mu_mass = {fReader, "BToMuMuPi_fit_mu_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_mu_phi = {fReader, "BToMuMuPi_fit_mu_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_mu_pt = {fReader, "BToMuMuPi_fit_mu_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_pi_eta = {fReader, "BToMuMuPi_fit_pi_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_pi_mass = {fReader, "BToMuMuPi_fit_pi_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_pi_phi = {fReader, "BToMuMuPi_fit_pi_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_pi_pt = {fReader, "BToMuMuPi_fit_pi_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_cos2D = {fReader, "BToMuMuPi_hnl_cos2D"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_eta = {fReader, "BToMuMuPi_hnl_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_iso03 = {fReader, "BToMuMuPi_hnl_iso03"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_iso03_close = {fReader, "BToMuMuPi_hnl_iso03_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_iso04 = {fReader, "BToMuMuPi_hnl_iso04"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_iso04_close = {fReader, "BToMuMuPi_hnl_iso04_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_mass = {fReader, "BToMuMuPi_hnl_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_masserr = {fReader, "BToMuMuPi_hnl_masserr"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_phi = {fReader, "BToMuMuPi_hnl_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_pt = {fReader, "BToMuMuPi_hnl_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_mass = {fReader, "BToMuMuPi_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_phi = {fReader, "BToMuMuPi_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_DCASig = {fReader, "BToMuMuPi_pi_DCASig"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_dxy = {fReader, "BToMuMuPi_pi_dxy"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_dxyS = {fReader, "BToMuMuPi_pi_dxyS"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_dz = {fReader, "BToMuMuPi_pi_dz"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_dzS = {fReader, "BToMuMuPi_pi_dzS"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_iso03 = {fReader, "BToMuMuPi_pi_iso03"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_iso03_close = {fReader, "BToMuMuPi_pi_iso03_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_iso04 = {fReader, "BToMuMuPi_pi_iso04"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_iso04_close = {fReader, "BToMuMuPi_pi_iso04_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_mu_vzdiff = {fReader, "BToMuMuPi_pi_mu_vzdiff"};
   TTreeReaderArray<Float_t> BToMuMuPi_pt = {fReader, "BToMuMuPi_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_dxy = {fReader, "BToMuMuPi_sel_mu_dxy"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_dz = {fReader, "BToMuMuPi_sel_mu_dz"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_ip3d = {fReader, "BToMuMuPi_sel_mu_ip3d"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_isLoose = {fReader, "BToMuMuPi_sel_mu_isLoose"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_isMedium = {fReader, "BToMuMuPi_sel_mu_isMedium"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_isSoft = {fReader, "BToMuMuPi_sel_mu_isSoft"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_isTight = {fReader, "BToMuMuPi_sel_mu_isTight"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso03 = {fReader, "BToMuMuPi_sel_mu_iso03"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso03_close = {fReader, "BToMuMuPi_sel_mu_iso03_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso04 = {fReader, "BToMuMuPi_sel_mu_iso04"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso04_close = {fReader, "BToMuMuPi_sel_mu_iso04_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_sip3d = {fReader, "BToMuMuPi_sel_mu_sip3d"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_chi2 = {fReader, "BToMuMuPi_sv_chi2"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_lxy = {fReader, "BToMuMuPi_sv_lxy"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_lxy_sig = {fReader, "BToMuMuPi_sv_lxy_sig"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_lxye = {fReader, "BToMuMuPi_sv_lxye"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_prob = {fReader, "BToMuMuPi_sv_prob"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_x = {fReader, "BToMuMuPi_sv_x"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_xe = {fReader, "BToMuMuPi_sv_xe"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_y = {fReader, "BToMuMuPi_sv_y"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_ye = {fReader, "BToMuMuPi_sv_ye"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_z = {fReader, "BToMuMuPi_sv_z"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_ze = {fReader, "BToMuMuPi_sv_ze"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_dxy = {fReader, "BToMuMuPi_trg_mu_dxy"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_dz = {fReader, "BToMuMuPi_trg_mu_dz"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_eta = {fReader, "BToMuMuPi_trg_mu_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_ip3d = {fReader, "BToMuMuPi_trg_mu_ip3d"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso03 = {fReader, "BToMuMuPi_trg_mu_iso03"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso03_close = {fReader, "BToMuMuPi_trg_mu_iso03_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso04 = {fReader, "BToMuMuPi_trg_mu_iso04"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso04_close = {fReader, "BToMuMuPi_trg_mu_iso04_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_phi = {fReader, "BToMuMuPi_trg_mu_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_pt = {fReader, "BToMuMuPi_trg_mu_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_sip3d = {fReader, "BToMuMuPi_trg_mu_sip3d"};
   TTreeReaderArray<Int_t> BToMuMuPi_charge = {fReader, "BToMuMuPi_charge"};
   TTreeReaderArray<Int_t> BToMuMuPi_hnl_charge = {fReader, "BToMuMuPi_hnl_charge"};
   TTreeReaderArray<Int_t> BToMuMuPi_isMatched = {fReader, "BToMuMuPi_isMatched"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_pi_genIdx = {fReader, "BToMuMuPi_matching_pi_genIdx"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_pi_motherPdgId = {fReader, "BToMuMuPi_matching_pi_motherPdgId"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_sel_mu_genIdx = {fReader, "BToMuMuPi_matching_sel_mu_genIdx"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_sel_mu_motherPdgId = {fReader, "BToMuMuPi_matching_sel_mu_motherPdgId"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_trg_mu_genIdx = {fReader, "BToMuMuPi_matching_trg_mu_genIdx"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_trg_mu_motherPdgId = {fReader, "BToMuMuPi_matching_trg_mu_motherPdgId"};
   TTreeReaderArray<Int_t> BToMuMuPi_pdgId = {fReader, "BToMuMuPi_pdgId"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_idx = {fReader, "BToMuMuPi_pi_idx"};
   TTreeReaderArray<Int_t> BToMuMuPi_sel_mu_idx = {fReader, "BToMuMuPi_sel_mu_idx"};
   TTreeReaderArray<Int_t> BToMuMuPi_trg_mu_idx = {fReader, "BToMuMuPi_trg_mu_idx"};
   TTreeReaderValue<UInt_t> nMuon = {fReader, "nMuon"};
   TTreeReaderArray<Float_t> Muon_dxy = {fReader, "Muon_dxy"};
   TTreeReaderArray<Float_t> Muon_dxyErr = {fReader, "Muon_dxyErr"};
   TTreeReaderArray<Float_t> Muon_dz = {fReader, "Muon_dz"};
   TTreeReaderArray<Float_t> Muon_dzErr = {fReader, "Muon_dzErr"};
   TTreeReaderArray<Float_t> Muon_eta = {fReader, "Muon_eta"};
   TTreeReaderArray<Float_t> Muon_ip3d = {fReader, "Muon_ip3d"};
   TTreeReaderArray<Float_t> Muon_mass = {fReader, "Muon_mass"};
   TTreeReaderArray<Float_t> Muon_pfRelIso03_all = {fReader, "Muon_pfRelIso03_all"};
   TTreeReaderArray<Float_t> Muon_pfRelIso04_all = {fReader, "Muon_pfRelIso04_all"};
   TTreeReaderArray<Float_t> Muon_phi = {fReader, "Muon_phi"};
   TTreeReaderArray<Float_t> Muon_pt = {fReader, "Muon_pt"};
   TTreeReaderArray<Float_t> Muon_ptErr = {fReader, "Muon_ptErr"};
   TTreeReaderArray<Float_t> Muon_sip3d = {fReader, "Muon_sip3d"};
   TTreeReaderArray<Float_t> Muon_vx = {fReader, "Muon_vx"};
   TTreeReaderArray<Float_t> Muon_vy = {fReader, "Muon_vy"};
   TTreeReaderArray<Float_t> Muon_vz = {fReader, "Muon_vz"};
   TTreeReaderArray<Int_t> Muon_charge = {fReader, "Muon_charge"};
   TTreeReaderArray<Int_t> Muon_isTriggering = {fReader, "Muon_isTriggering"};
   TTreeReaderArray<Int_t> Muon_pdgId = {fReader, "Muon_pdgId"};
   TTreeReaderArray<Bool_t> Muon_isGlobal = {fReader, "Muon_isGlobal"};
   TTreeReaderArray<Bool_t> Muon_isPFcand = {fReader, "Muon_isPFcand"};
   TTreeReaderArray<Bool_t> Muon_isTracker = {fReader, "Muon_isTracker"};
   TTreeReaderArray<Bool_t> Muon_mediumId = {fReader, "Muon_mediumId"};
   TTreeReaderArray<UChar_t> Muon_pfIsoId = {fReader, "Muon_pfIsoId"};
   TTreeReaderArray<Bool_t> Muon_softId = {fReader, "Muon_softId"};
   TTreeReaderArray<Bool_t> Muon_tightId = {fReader, "Muon_tightId"};
   TTreeReaderArray<UChar_t> Muon_tkIsoId = {fReader, "Muon_tkIsoId"};
   TTreeReaderArray<Bool_t> Muon_triggerIdLoose = {fReader, "Muon_triggerIdLoose"};
   TTreeReaderValue<UInt_t> nTriggerMuon = {fReader, "nTriggerMuon"};
   TTreeReaderArray<Float_t> TriggerMuon_eta = {fReader, "TriggerMuon_eta"};
   TTreeReaderArray<Float_t> TriggerMuon_mass = {fReader, "TriggerMuon_mass"};
   TTreeReaderArray<Float_t> TriggerMuon_phi = {fReader, "TriggerMuon_phi"};
   TTreeReaderArray<Float_t> TriggerMuon_pt = {fReader, "TriggerMuon_pt"};
   TTreeReaderArray<Float_t> TriggerMuon_vx = {fReader, "TriggerMuon_vx"};
   TTreeReaderArray<Float_t> TriggerMuon_vy = {fReader, "TriggerMuon_vy"};
   TTreeReaderArray<Float_t> TriggerMuon_vz = {fReader, "TriggerMuon_vz"};
   TTreeReaderArray<Int_t> TriggerMuon_charge = {fReader, "TriggerMuon_charge"};
   TTreeReaderArray<Int_t> TriggerMuon_pdgId = {fReader, "TriggerMuon_pdgId"};
   //TTreeReaderArray<Int_t> TriggerMuon_trgMuonIndex = {fReader, "TriggerMuon_trgMuonIndex"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetAll = {fReader, "fixedGridRhoFastjetAll"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentral = {fReader, "fixedGridRhoFastjetCentral"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralCalo = {fReader, "fixedGridRhoFastjetCentralCalo"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralChargedPileUp = {fReader, "fixedGridRhoFastjetCentralChargedPileUp"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralNeutral = {fReader, "fixedGridRhoFastjetCentralNeutral"};
   TTreeReaderValue<UInt_t> nProbeTracks = {fReader, "nProbeTracks"};
   TTreeReaderArray<Float_t> ProbeTracks_DCASig = {fReader, "ProbeTracks_DCASig"};
   TTreeReaderArray<Float_t> ProbeTracks_dxy = {fReader, "ProbeTracks_dxy"};
   TTreeReaderArray<Float_t> ProbeTracks_dxyS = {fReader, "ProbeTracks_dxyS"};
   TTreeReaderArray<Float_t> ProbeTracks_dz = {fReader, "ProbeTracks_dz"};
   TTreeReaderArray<Float_t> ProbeTracks_dzS = {fReader, "ProbeTracks_dzS"};
   TTreeReaderArray<Float_t> ProbeTracks_eta = {fReader, "ProbeTracks_eta"};
   TTreeReaderArray<Float_t> ProbeTracks_mass = {fReader, "ProbeTracks_mass"};
   TTreeReaderArray<Float_t> ProbeTracks_phi = {fReader, "ProbeTracks_phi"};
   TTreeReaderArray<Float_t> ProbeTracks_pt = {fReader, "ProbeTracks_pt"};
   TTreeReaderArray<Float_t> ProbeTracks_vx = {fReader, "ProbeTracks_vx"};
   TTreeReaderArray<Float_t> ProbeTracks_vy = {fReader, "ProbeTracks_vy"};
   TTreeReaderArray<Float_t> ProbeTracks_vz = {fReader, "ProbeTracks_vz"};
   TTreeReaderArray<Int_t> ProbeTracks_charge = {fReader, "ProbeTracks_charge"};
   TTreeReaderArray<Int_t> ProbeTracks_isLostTrk = {fReader, "ProbeTracks_isLostTrk"};
   TTreeReaderArray<Int_t> ProbeTracks_isPacked = {fReader, "ProbeTracks_isPacked"};
   TTreeReaderArray<Int_t> ProbeTracks_nValidHits = {fReader, "ProbeTracks_nValidHits"};
   TTreeReaderArray<Int_t> ProbeTracks_pdgId = {fReader, "ProbeTracks_pdgId"};
   TTreeReaderArray<Bool_t> ProbeTracks_isMatchedToEle = {fReader, "ProbeTracks_isMatchedToEle"};
   TTreeReaderArray<Bool_t> ProbeTracks_isMatchedToLooseMuon = {fReader, "ProbeTracks_isMatchedToLooseMuon"};
   TTreeReaderArray<Bool_t> ProbeTracks_isMatchedToMediumMuon = {fReader, "ProbeTracks_isMatchedToMediumMuon"};
   TTreeReaderArray<Bool_t> ProbeTracks_isMatchedToMuon = {fReader, "ProbeTracks_isMatchedToMuon"};
   TTreeReaderArray<Bool_t> ProbeTracks_isMatchedToSoftMuon = {fReader, "ProbeTracks_isMatchedToSoftMuon"};
   TTreeReaderValue<UChar_t> HLT_Mu7_IP4 = {fReader, "HLT_Mu7_IP4"};
   TTreeReaderValue<UChar_t> HLT_Mu8_IP6 = {fReader, "HLT_Mu8_IP6"};
   TTreeReaderValue<UChar_t> HLT_Mu8_IP5 = {fReader, "HLT_Mu8_IP5"};
   TTreeReaderValue<UChar_t> HLT_Mu8_IP3 = {fReader, "HLT_Mu8_IP3"};
   TTreeReaderValue<UChar_t> HLT_Mu8p5_IP3p5 = {fReader, "HLT_Mu8p5_IP3p5"};
   TTreeReaderValue<UChar_t> HLT_Mu9_IP6 = {fReader, "HLT_Mu9_IP6"};
   TTreeReaderValue<UChar_t> HLT_Mu9_IP5 = {fReader, "HLT_Mu9_IP5"};
   TTreeReaderValue<UChar_t> HLT_Mu9_IP4 = {fReader, "HLT_Mu9_IP4"};
   TTreeReaderValue<UChar_t> HLT_Mu10p5_IP3p5 = {fReader, "HLT_Mu10p5_IP3p5"};
   TTreeReaderValue<UChar_t> HLT_Mu12_IP6 = {fReader, "HLT_Mu12_IP6"};
   TTreeReaderValue<UChar_t> L1_SingleMu7er1p5 = {fReader, "L1_SingleMu7er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu8er1p5 = {fReader, "L1_SingleMu8er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu9er1p5 = {fReader, "L1_SingleMu9er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu10er1p5 = {fReader, "L1_SingleMu10er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu12er1p5 = {fReader, "L1_SingleMu12er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu22 = {fReader, "L1_SingleMu22"};
   TTreeReaderValue<UInt_t> nTrigObj = {fReader, "nTrigObj"};
   TTreeReaderArray<Float_t> TrigObj_pt = {fReader, "TrigObj_pt"};
   TTreeReaderArray<Float_t> TrigObj_eta = {fReader, "TrigObj_eta"};
   TTreeReaderArray<Float_t> TrigObj_phi = {fReader, "TrigObj_phi"};
   TTreeReaderArray<Float_t> TrigObj_l1pt = {fReader, "TrigObj_l1pt"};
   TTreeReaderArray<Float_t> TrigObj_l1pt_2 = {fReader, "TrigObj_l1pt_2"};
   TTreeReaderArray<Float_t> TrigObj_l2pt = {fReader, "TrigObj_l2pt"};
   TTreeReaderArray<Int_t> TrigObj_id = {fReader, "TrigObj_id"};
   TTreeReaderArray<Int_t> TrigObj_l1iso = {fReader, "TrigObj_l1iso"};
   TTreeReaderArray<Int_t> TrigObj_l1charge = {fReader, "TrigObj_l1charge"};
   TTreeReaderArray<Int_t> TrigObj_filterBits = {fReader, "TrigObj_filterBits"};
   TTreeReaderValue<UInt_t> nOtherPV = {fReader, "nOtherPV"};
   TTreeReaderArray<Float_t> OtherPV_z = {fReader, "OtherPV_z"};
   TTreeReaderValue<Float_t> PV_ndof = {fReader, "PV_ndof"};
   TTreeReaderValue<Float_t> PV_x = {fReader, "PV_x"};
   TTreeReaderValue<Float_t> PV_y = {fReader, "PV_y"};
   TTreeReaderValue<Float_t> PV_z = {fReader, "PV_z"};
   TTreeReaderValue<Float_t> PV_chi2 = {fReader, "PV_chi2"};
   TTreeReaderValue<Float_t> PV_score = {fReader, "PV_score"};
   TTreeReaderValue<Int_t> PV_npvs = {fReader, "PV_npvs"};
   TTreeReaderValue<Int_t> PV_npvsGood = {fReader, "PV_npvsGood"};
   TTreeReaderValue<UInt_t> nSV = {fReader, "nSV"};
   TTreeReaderArray<Float_t> SV_dlen = {fReader, "SV_dlen"};
   TTreeReaderArray<Float_t> SV_dlenSig = {fReader, "SV_dlenSig"};
   TTreeReaderArray<Float_t> SV_pAngle = {fReader, "SV_pAngle"};
   TTreeReaderArray<Float_t> SV_chi2 = {fReader, "SV_chi2"};
   TTreeReaderArray<Float_t> SV_eta = {fReader, "SV_eta"};
   TTreeReaderArray<Float_t> SV_mass = {fReader, "SV_mass"};
   TTreeReaderArray<Float_t> SV_ndof = {fReader, "SV_ndof"};
   TTreeReaderArray<Float_t> SV_phi = {fReader,"SV_phi"};
   TTreeReaderArray<Float_t> SV_pt = {fReader, "SV_pt"};
   TTreeReaderArray<Float_t> SV_x = {fReader, "SV_x"};
   TTreeReaderArray<Float_t> SV_y = {fReader, "SV_y"};
   TTreeReaderArray<Float_t> SV_z = {fReader, "SV_z"};
   // branches intentionally pointed to wrong variable, to avoid crash when running on data
   TTreeReaderValue<UInt_t> nGenPart = {fReader, "nBToMuMuPi"};
   TTreeReaderArray<Float_t> GenPart_eta = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_mass = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_phi = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_pt = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_vx = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_vy = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_vz = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {fReader, "BToMuMuPi_charge"};
   TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "BToMuMuPi_charge"};
   TTreeReaderArray<Int_t> GenPart_status = {fReader, "BToMuMuPi_charge"};
   TTreeReaderArray<Int_t> GenPart_statusFlags = {fReader, "BToMuMuPi_charge"};
   TTreeReaderArray<Int_t> Muon_genPartIdx = {fReader, "BToMuMuPi_charge"};


   NanoDumper(TTree * /*tree*/ =0) { }
   virtual ~NanoDumper() { }
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

   // additional functions
   virtual vector<pair<int,float>> createPairWithDesc(const UInt_t&, const TTreeReaderArray<Float_t>&);
   virtual vector<pair<int,float>> updatePairWithDesc(vector<pair<int,float>> the_ini_pair, const TTreeReaderArray<Int_t>& charge);

   // output file
   TFile* my_file;  

   // some branches will be added only if sample is MC
   Bool_t isMC;

   // some option
   // this option is intentionally hardcoded
   Bool_t do_fillhistograms = false;

   // trees to fill
   TTree* signal_tree;
   TTree* control_tree;
   
   // filling variables
   ULong64_t the_event;
   UInt_t the_run;
   UInt_t the_lumi;

   Int_t the_pv_npvs;

   UChar_t the_hlt_mu7_ip4;
   UChar_t the_hlt_mu8_ip6;
   UChar_t the_hlt_mu8_ip5;
   UChar_t the_hlt_mu8_ip3;
   UChar_t the_hlt_mu8p5_ip3p5;
   UChar_t the_hlt_mu9_ip6;
   UChar_t the_hlt_mu9_ip5;
   UChar_t the_hlt_mu9_ip4;
   UChar_t the_hlt_mu10p5_ip3p5;
   UChar_t the_hlt_mu12_ip6;

   Float_t the_sig_b_pt;
   Float_t the_sig_b_eta;
   Float_t the_sig_b_phi;
   Float_t the_sig_b_mass;
   Int_t the_sig_b_charge;
   Int_t the_sig_b_pdgid;

   Float_t the_sig_hnl_pt;
   Float_t the_sig_hnl_eta;
   Float_t the_sig_hnl_phi;
   Float_t the_sig_hnl_mass;
   Int_t the_sig_hnl_charge;
   Float_t the_sig_hnl_cos2d;
   Float_t the_sig_hnl_iso03;
   Float_t the_sig_hnl_iso03_close;
   Float_t the_sig_hnl_iso04;
   Float_t the_sig_hnl_iso04_close;

   Float_t the_sig_trgmu_pt;
   Float_t the_sig_trgmu_eta;
   Float_t the_sig_trgmu_phi;
   Int_t the_sig_trgmu_charge;
   Float_t the_sig_trgmu_dxy;
   Float_t the_sig_trgmu_dz;
   Float_t the_sig_trgmu_ip3d;
   Float_t the_sig_trgmu_ip3dsig;
   Float_t the_sig_trgmu_iso03;
   Float_t the_sig_trgmu_iso03_close;
   Float_t the_sig_trgmu_iso04;
   Float_t the_sig_trgmu_iso04_close;
   Float_t the_sig_trgmu_istriggering;

   Float_t the_sig_mu_pt;
   Float_t the_sig_mu_eta;
   Float_t the_sig_mu_phi;
   Int_t the_sig_mu_charge;
   Float_t the_sig_mu_dxy;
   Float_t the_sig_mu_dz;
   Float_t the_sig_mu_ip3d;
   Float_t the_sig_mu_ip3dsig;
   Float_t the_sig_mu_iso03;
   Float_t the_sig_mu_iso03_close;
   Float_t the_sig_mu_iso04;
   Float_t the_sig_mu_iso04_close;
   Float_t the_sig_mu_isloose;
   Float_t the_sig_mu_ismedium;
   Float_t the_sig_mu_istight;
   Float_t the_sig_mu_issoft;
   Float_t the_sig_mu_istriggering;

   Float_t the_sig_pi_pt;
   Float_t the_sig_pi_eta;
   Float_t the_sig_pi_phi;
   Int_t the_sig_pi_charge;
   Float_t the_sig_pi_dcasig;
   Float_t the_sig_pi_dxy;
   Float_t the_sig_pi_dz;
   Float_t the_sig_pi_dxysig;
   Float_t the_sig_pi_dzsig;
   Float_t the_sig_pi_iso03;
   Float_t the_sig_pi_iso03_close;
   Float_t the_sig_pi_iso04;
   Float_t the_sig_pi_iso04_close;

   Float_t the_sig_dimu_mass;
   Float_t the_sig_dimu_pt;
   Float_t the_sig_dimu_lxy;
   Float_t the_sig_dimu_lxyz;
   Float_t the_sig_dimu_vxdiff;
   Float_t the_sig_dimu_vydiff;
   Float_t the_sig_dimu_vzdiff;

   Float_t the_sig_deltar_mu_pi;
   Float_t the_sig_deltar_trgmu_hnl;

   Float_t the_sig_sv_chi2;
   Float_t the_sig_sv_lxy;
   Float_t the_sig_sv_pv_lxy;
   Float_t the_sig_sv_pv_lxyz;
   Float_t the_sig_sv_lxysig;
   Float_t the_sig_sv_prob;
   Float_t the_sig_sv_x;
   Float_t the_sig_sv_y;
   Float_t the_sig_sv_z;

   Float_t the_sig_pi_mu_vzdiff;

   Int_t the_sig_ismatched;

   Float_t the_gen_trgmu_mu_lxy;
   Float_t the_gen_trgmu_mu_lxyz;

   Float_t the_ctrl_b_pt;
   Float_t the_ctrl_b_eta;
   Float_t the_ctrl_b_phi;
   Float_t the_ctrl_b_mass;
   Int_t the_ctrl_b_charge;
   Int_t the_ctrl_b_pdgid;
   Float_t the_ctrl_b_cos2d;
   Float_t the_ctrl_b_iso03;
   Float_t the_ctrl_b_iso03_close;
   Float_t the_ctrl_b_iso04;
   Float_t the_ctrl_b_iso04_close;

   Float_t the_ctrl_k_pt;
   Float_t the_ctrl_k_eta;
   Float_t the_ctrl_k_phi;
   Int_t the_ctrl_k_charge;
   Float_t the_ctrl_k_iso03;
   Float_t the_ctrl_k_iso03_close;
   Float_t the_ctrl_k_iso04;
   Float_t the_ctrl_k_iso04_close;

   Float_t the_ctrl_l1_pt;
   Float_t the_ctrl_l1_eta;
   Float_t the_ctrl_l1_phi;
   Int_t the_ctrl_l1_charge;
   Float_t the_ctrl_l1_iso03;
   Float_t the_ctrl_l1_iso03_close;
   Float_t the_ctrl_l1_iso04;
   Float_t the_ctrl_l1_iso04_close;
   Float_t the_ctrl_l1_istriggering;

   Float_t the_ctrl_l2_pt;
   Float_t the_ctrl_l2_eta;
   Float_t the_ctrl_l2_phi;
   Int_t the_ctrl_l2_charge;
   Float_t the_ctrl_l2_iso03;
   Float_t the_ctrl_l2_iso03_close;
   Float_t the_ctrl_l2_iso04;
   Float_t the_ctrl_l2_iso04_close;
   Float_t the_ctrl_l2_istriggering;

   Float_t the_ctrl_dimu_mass;

   Float_t the_ctrl_sv_x;
   Float_t the_ctrl_sv_y;
   Float_t the_ctrl_sv_z;
   Float_t the_ctrl_sv_lxy;
   Float_t the_ctrl_sv_prob;

   Int_t the_ctrl_ismatched;

   // histograms
   TH1F* ctrlhist_ncand_perevent;
   TH1F* ctrlhist_ncand_matched_perevent;
   TH1F* ctrlhist_selection_efficiency_bpt_allevents;
   TH1F* ctrlhist_selection_efficiency_bpt_eventswithmultcands;
   TH1F* ctrlhist_selection_efficiency_svprob_allevents;
   TH1F* ctrlhist_selection_efficiency_svprob_eventswithmultcands;
   TH1F* ctrlhist_selection_efficiency_cos2d_allevents;
   TH1F* ctrlhist_selection_efficiency_cos2d_eventswithmultcands;
   TH1F* ctrlhist_selection_efficiency_kpt_allevents;
   TH1F* ctrlhist_selection_efficiency_kpt_eventswithmultcands;
   TH1F* ctrlhist_ncand_wtriggeringmuon;

   TH1F* sighist_ncand_perevent;
   TH1F* sighist_ncand_matched_perevent;
   TH1F* sighist_selection_efficiency_hnlpt_allevents;
   TH1F* sighist_selection_efficiency_hnlpt_eventswithmultcands;
   TH1F* sighist_selection_efficiency_bpt_allevents;
   TH1F* sighist_selection_efficiency_bpt_eventswithmultcands;
   TH1F* sighist_selection_efficiency_trgmupt_allevents;
   TH1F* sighist_selection_efficiency_trgmupt_eventswithmultcands;
   TH1F* sighist_selection_efficiency_pipt_allevents;
   TH1F* sighist_selection_efficiency_pipt_eventswithmultcands;
   TH1F* sighist_selection_efficiency_svprob_allevents;
   TH1F* sighist_selection_efficiency_svprob_eventswithmultcands;
   TH1F* sighist_selection_efficiency_svchi2_allevents;
   TH1F* sighist_selection_efficiency_svchi2_eventswithmultcands;
   TH1F* sighist_selection_efficiency_cos2d_allevents;
   TH1F* sighist_selection_efficiency_cos2d_eventswithmultcands;
   TH1F* sighist_selection_efficiency_hnliso4_allevents;
   TH1F* sighist_selection_efficiency_hnliso4_eventswithmultcands;
   TH1F* sighist_selection_efficiency_dr_allevents;
   TH1F* sighist_selection_efficiency_dr_eventswithmultcands;

   ClassDef(NanoDumper,0);
};

#endif

#ifdef NanoDumper_cxx
void NanoDumper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t NanoDumper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef NanoDumper_cxx
