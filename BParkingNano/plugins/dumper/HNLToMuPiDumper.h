//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug  4 13:32:30 2021 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: bparknano.root
//////////////////////////////////////////////////////////

#ifndef HNLToMuPiDumper_h
#define HNLToMuPiDumper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class HNLToMuPiDumper : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> nHNLToMuPi = {fReader, "nHNLToMuPi"};
   TTreeReaderArray<Float_t> HNLToMuPi_eta = {fReader, "HNLToMuPi_eta"};
   TTreeReaderArray<Float_t> HNLToMuPi_fit_mu_eta = {fReader, "HNLToMuPi_fit_mu_eta"};
   TTreeReaderArray<Float_t> HNLToMuPi_fit_mu_mass = {fReader, "HNLToMuPi_fit_mu_mass"};
   TTreeReaderArray<Float_t> HNLToMuPi_fit_mu_phi = {fReader, "HNLToMuPi_fit_mu_phi"};
   TTreeReaderArray<Float_t> HNLToMuPi_fit_mu_pt = {fReader, "HNLToMuPi_fit_mu_pt"};
   TTreeReaderArray<Float_t> HNLToMuPi_fit_pi_eta = {fReader, "HNLToMuPi_fit_pi_eta"};
   TTreeReaderArray<Float_t> HNLToMuPi_fit_pi_mass = {fReader, "HNLToMuPi_fit_pi_mass"};
   TTreeReaderArray<Float_t> HNLToMuPi_fit_pi_phi = {fReader, "HNLToMuPi_fit_pi_phi"};
   TTreeReaderArray<Float_t> HNLToMuPi_fit_pi_pt = {fReader, "HNLToMuPi_fit_pi_pt"};
   TTreeReaderArray<Float_t> HNLToMuPi_gen_lxy = {fReader, "HNLToMuPi_gen_lxy"};
   TTreeReaderArray<Float_t> HNLToMuPi_hnl_cos2D = {fReader, "HNLToMuPi_hnl_cos2D"};
   TTreeReaderArray<Float_t> HNLToMuPi_hnl_eta = {fReader, "HNLToMuPi_hnl_eta"};
   TTreeReaderArray<Float_t> HNLToMuPi_hnl_mass = {fReader, "HNLToMuPi_hnl_mass"};
   TTreeReaderArray<Float_t> HNLToMuPi_hnl_mass_unfitted = {fReader, "HNLToMuPi_hnl_mass_unfitted"};
   TTreeReaderArray<Float_t> HNLToMuPi_diff_mass = {fReader, "HNLToMuPi_diff_mass"};
   TTreeReaderArray<Float_t> HNLToMuPi_hnl_masserr = {fReader, "HNLToMuPi_hnl_masserr"};
   TTreeReaderArray<Float_t> HNLToMuPi_hnl_phi = {fReader, "HNLToMuPi_hnl_phi"};
   TTreeReaderArray<Float_t> HNLToMuPi_hnl_pt = {fReader, "HNLToMuPi_hnl_pt"};
   TTreeReaderArray<Float_t> HNLToMuPi_lxy_reco_gen_reldiff = {fReader, "HNLToMuPi_lxy_reco_gen_reldiff"};
   TTreeReaderArray<Float_t> HNLToMuPi_mass = {fReader, "HNLToMuPi_mass"};
   TTreeReaderArray<Float_t> HNLToMuPi_mupi_mass_reco_gen_reldiff = {fReader, "HNLToMuPi_mupi_mass_reco_gen_reldiff"};
   TTreeReaderArray<Float_t> HNLToMuPi_phi = {fReader, "HNLToMuPi_phi"};
   TTreeReaderArray<Float_t> HNLToMuPi_pt = {fReader, "HNLToMuPi_pt"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_chi2 = {fReader, "HNLToMuPi_sv_chi2"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_lxy = {fReader, "HNLToMuPi_sv_lxy"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_lxy_sig = {fReader, "HNLToMuPi_sv_lxy_sig"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_lxye = {fReader, "HNLToMuPi_sv_lxye"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_prob = {fReader, "HNLToMuPi_sv_prob"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_x = {fReader, "HNLToMuPi_sv_x"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_xe = {fReader, "HNLToMuPi_sv_xe"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_y = {fReader, "HNLToMuPi_sv_y"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_ye = {fReader, "HNLToMuPi_sv_ye"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_z = {fReader, "HNLToMuPi_sv_z"};
   TTreeReaderArray<Float_t> HNLToMuPi_sv_ze = {fReader, "HNLToMuPi_sv_ze"};
   TTreeReaderArray<Int_t> HNLToMuPi_charge = {fReader, "HNLToMuPi_charge"};
   TTreeReaderArray<Int_t> HNLToMuPi_hnl_charge = {fReader, "HNLToMuPi_hnl_charge"};
   TTreeReaderArray<Int_t> HNLToMuPi_isMatched = {fReader, "HNLToMuPi_isMatched"};
   TTreeReaderArray<Int_t> HNLToMuPi_matching_pi_genIdx = {fReader, "HNLToMuPi_matching_pi_genIdx"};
   TTreeReaderArray<Int_t> HNLToMuPi_matching_pi_motherPdgId = {fReader, "HNLToMuPi_matching_pi_motherPdgId"};
   TTreeReaderArray<Int_t> HNLToMuPi_matching_sel_mu_genIdx = {fReader, "HNLToMuPi_matching_sel_mu_genIdx"};
   TTreeReaderArray<Int_t> HNLToMuPi_matching_sel_mu_motherPdgId = {fReader, "HNLToMuPi_matching_sel_mu_motherPdgId"};
   TTreeReaderArray<Int_t> HNLToMuPi_pdgId = {fReader, "HNLToMuPi_pdgId"};
   TTreeReaderArray<Int_t> HNLToMuPi_pi_idx = {fReader, "HNLToMuPi_pi_idx"};
   TTreeReaderArray<Int_t> HNLToMuPi_pi_isMatched = {fReader, "HNLToMuPi_pi_isMatched"};
   TTreeReaderArray<Int_t> HNLToMuPi_sel_mu_idx = {fReader, "HNLToMuPi_sel_mu_idx"};
   TTreeReaderArray<Int_t> HNLToMuPi_sel_mu_isMatched = {fReader, "HNLToMuPi_sel_mu_isMatched"};
   TTreeReaderValue<UInt_t> nGenPart = {fReader, "nGenPart"};
   TTreeReaderArray<Float_t> GenPart_beta = {fReader, "GenPart_beta"};
   TTreeReaderArray<Float_t> GenPart_eta = {fReader, "GenPart_eta"};
   TTreeReaderArray<Float_t> GenPart_gamma = {fReader, "GenPart_gamma"};
   TTreeReaderArray<Float_t> GenPart_mass = {fReader, "GenPart_mass"};
   TTreeReaderArray<Float_t> GenPart_p = {fReader, "GenPart_p"};
   TTreeReaderArray<Float_t> GenPart_phi = {fReader, "GenPart_phi"};
   TTreeReaderArray<Float_t> GenPart_pt = {fReader, "GenPart_pt"};
   TTreeReaderArray<Float_t> GenPart_vx = {fReader, "GenPart_vx"};
   TTreeReaderArray<Float_t> GenPart_vy = {fReader, "GenPart_vy"};
   TTreeReaderArray<Float_t> GenPart_vz = {fReader, "GenPart_vz"};
   TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {fReader, "GenPart_genPartIdxMother"};
   TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "GenPart_pdgId"};
   TTreeReaderArray<Int_t> GenPart_status = {fReader, "GenPart_status"};
   TTreeReaderArray<Int_t> GenPart_statusFlags = {fReader, "GenPart_statusFlags"};
   TTreeReaderValue<Float_t> Generator_binvar = {fReader, "Generator_binvar"};
   TTreeReaderValue<Float_t> Generator_scalePDF = {fReader, "Generator_scalePDF"};
   TTreeReaderValue<Float_t> Generator_weight = {fReader, "Generator_weight"};
   TTreeReaderValue<Float_t> Generator_x1 = {fReader, "Generator_x1"};
   TTreeReaderValue<Float_t> Generator_x2 = {fReader, "Generator_x2"};
   TTreeReaderValue<Float_t> Generator_xpdf1 = {fReader, "Generator_xpdf1"};
   TTreeReaderValue<Float_t> Generator_xpdf2 = {fReader, "Generator_xpdf2"};
   TTreeReaderValue<Int_t> Generator_id1 = {fReader, "Generator_id1"};
   TTreeReaderValue<Int_t> Generator_id2 = {fReader, "Generator_id2"};
   TTreeReaderValue<Float_t> genWeight = {fReader, "genWeight"};
   TTreeReaderValue<UInt_t> nPSWeight = {fReader, "nPSWeight"};
   TTreeReaderArray<Float_t> PSWeight = {fReader, "PSWeight"};
   TTreeReaderValue<UInt_t> nMuon = {fReader, "nMuon"};
   TTreeReaderArray<Float_t> Muon_caloCompatibility = {fReader, "Muon_caloCompatibility"};
   TTreeReaderArray<Float_t> Muon_dxy = {fReader, "Muon_dxy"};
   TTreeReaderArray<Float_t> Muon_dxyS = {fReader, "Muon_dxyS"};
   TTreeReaderArray<Float_t> Muon_dz = {fReader, "Muon_dz"};
   TTreeReaderArray<Float_t> Muon_dzS = {fReader, "Muon_dzS"};
   TTreeReaderArray<Float_t> Muon_eta = {fReader, "Muon_eta"};
   TTreeReaderArray<Float_t> Muon_globalNormalisedChi2 = {fReader, "Muon_globalNormalisedChi2"};
   TTreeReaderArray<Float_t> Muon_kinkFinderChi2 = {fReader, "Muon_kinkFinderChi2"};
   TTreeReaderArray<Float_t> Muon_localPositionChi2 = {fReader, "Muon_localPositionChi2"};
   TTreeReaderArray<Float_t> Muon_mass = {fReader, "Muon_mass"};
   TTreeReaderArray<Float_t> Muon_matched_dpt = {fReader, "Muon_matched_dpt"};
   TTreeReaderArray<Float_t> Muon_matched_dr = {fReader, "Muon_matched_dr"};
   TTreeReaderArray<Float_t> Muon_pfiso03Rel_all = {fReader, "Muon_pfiso03Rel_all"};
   TTreeReaderArray<Float_t> Muon_pfiso03_all = {fReader, "Muon_pfiso03_all"};
   TTreeReaderArray<Float_t> Muon_pfiso04Rel_all = {fReader, "Muon_pfiso04Rel_all"};
   TTreeReaderArray<Float_t> Muon_pfiso04_all = {fReader, "Muon_pfiso04_all"};
   TTreeReaderArray<Float_t> Muon_phi = {fReader, "Muon_phi"};
   TTreeReaderArray<Float_t> Muon_pt = {fReader, "Muon_pt"};
   TTreeReaderArray<Float_t> Muon_segmentCompatibility = {fReader, "Muon_segmentCompatibility"};
   TTreeReaderArray<Float_t> Muon_validHitFraction = {fReader, "Muon_validHitFraction"};
   TTreeReaderArray<Float_t> Muon_vx = {fReader, "Muon_vx"};
   TTreeReaderArray<Float_t> Muon_vy = {fReader, "Muon_vy"};
   TTreeReaderArray<Float_t> Muon_vz = {fReader, "Muon_vz"};
   TTreeReaderArray<Int_t> Muon_charge = {fReader, "Muon_charge"};
   TTreeReaderArray<Int_t> Muon_fired_DST_DoubleMu1_noVtx_CaloScouting_v2 = {fReader, "Muon_fired_DST_DoubleMu1_noVtx_CaloScouting_v2"};
   TTreeReaderArray<Int_t> Muon_fired_DST_DoubleMu3_noVtx_CaloScouting_v6 = {fReader, "Muon_fired_DST_DoubleMu3_noVtx_CaloScouting_v6"};
   TTreeReaderArray<Int_t> Muon_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3 = {fReader, "Muon_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13 = {fReader, "Muon_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1 = {fReader, "Muon_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_L2Mu23NoVtx_2Cha_v1 = {fReader, "Muon_fired_HLT_L2Mu23NoVtx_2Cha_v1"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu10p5_IP3p5 = {fReader, "Muon_fired_HLT_Mu10p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu12_IP6 = {fReader, "Muon_fired_HLT_Mu12_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu7_IP4 = {fReader, "Muon_fired_HLT_Mu7_IP4"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu7p5_Track7_Jpsi_v11 = {fReader, "Muon_fired_HLT_Mu7p5_Track7_Jpsi_v11"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP3 = {fReader, "Muon_fired_HLT_Mu8_IP3"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP5 = {fReader, "Muon_fired_HLT_Mu8_IP5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP6 = {fReader, "Muon_fired_HLT_Mu8_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_v1 = {fReader, "Muon_fired_HLT_Mu8_v1"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_v12 = {fReader, "Muon_fired_HLT_Mu8_v12"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8p5_IP3p5 = {fReader, "Muon_fired_HLT_Mu8p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP4 = {fReader, "Muon_fired_HLT_Mu9_IP4"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP5 = {fReader, "Muon_fired_HLT_Mu9_IP5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP6 = {fReader, "Muon_fired_HLT_Mu9_IP6"};
   TTreeReaderArray<Int_t> Muon_inTimeMuon = {fReader, "Muon_inTimeMuon"};
   TTreeReaderArray<Int_t> Muon_isDSAMuon = {fReader, "Muon_isDSAMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalMuon = {fReader, "Muon_isGlobalMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalNotTrackerMuon = {fReader, "Muon_isGlobalNotTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalOrTrackerMuon = {fReader, "Muon_isGlobalOrTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isPF = {fReader, "Muon_isPF"};
   TTreeReaderArray<Int_t> Muon_isSlimmedMuon = {fReader, "Muon_isSlimmedMuon"};
   TTreeReaderArray<Int_t> Muon_isTrackerMuon = {fReader, "Muon_isTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isTrackerNotGlobalMuon = {fReader, "Muon_isTrackerNotGlobalMuon"};
   TTreeReaderArray<Int_t> Muon_isTriggering = {fReader, "Muon_isTriggering"};
   TTreeReaderArray<Int_t> Muon_isTriggeringBPark = {fReader, "Muon_isTriggeringBPark"};
   TTreeReaderArray<Int_t> Muon_looseId = {fReader, "Muon_looseId"};
   TTreeReaderArray<Int_t> Muon_mediumId = {fReader, "Muon_mediumId"};
   TTreeReaderArray<Int_t> Muon_numberOfPixelLayers = {fReader, "Muon_numberOfPixelLayers"};
   TTreeReaderArray<Int_t> Muon_numberOfStations = {fReader, "Muon_numberOfStations"};
   TTreeReaderArray<Int_t> Muon_numberOfTrackerLayers = {fReader, "Muon_numberOfTrackerLayers"};
   TTreeReaderArray<Int_t> Muon_numberOfValidMuonHits = {fReader, "Muon_numberOfValidMuonHits"};
   TTreeReaderArray<Int_t> Muon_numberOfValidPixelHits = {fReader, "Muon_numberOfValidPixelHits"};
   TTreeReaderArray<Int_t> Muon_pdgId = {fReader, "Muon_pdgId"};
   TTreeReaderArray<Int_t> Muon_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2 = {fReader, "Muon_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2"};
   TTreeReaderArray<Int_t> Muon_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6 = {fReader, "Muon_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6"};
   TTreeReaderArray<Int_t> Muon_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3 = {fReader, "Muon_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13 = {fReader, "Muon_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1 = {fReader, "Muon_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_L2Mu23NoVtx_2Cha_v1 = {fReader, "Muon_prescale_HLT_L2Mu23NoVtx_2Cha_v1"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu10p5_IP3p5 = {fReader, "Muon_prescale_HLT_Mu10p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu12_IP6 = {fReader, "Muon_prescale_HLT_Mu12_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu7_IP4 = {fReader, "Muon_prescale_HLT_Mu7_IP4"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu7p5_Track7_Jpsi_v11 = {fReader, "Muon_prescale_HLT_Mu7p5_Track7_Jpsi_v11"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP3 = {fReader, "Muon_prescale_HLT_Mu8_IP3"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP5 = {fReader, "Muon_prescale_HLT_Mu8_IP5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP6 = {fReader, "Muon_prescale_HLT_Mu8_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_v1 = {fReader, "Muon_prescale_HLT_Mu8_v1"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_v12 = {fReader, "Muon_prescale_HLT_Mu8_v12"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8p5_IP3p5 = {fReader, "Muon_prescale_HLT_Mu8p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP4 = {fReader, "Muon_prescale_HLT_Mu9_IP4"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP5 = {fReader, "Muon_prescale_HLT_Mu9_IP5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP6 = {fReader, "Muon_prescale_HLT_Mu9_IP6"};
   TTreeReaderArray<Int_t> Muon_softId = {fReader, "Muon_softId"};
   TTreeReaderArray<Int_t> Muon_tightId = {fReader, "Muon_tightId"};
   TTreeReaderArray<Int_t> Muon_trackerHighPurityFlag = {fReader, "Muon_trackerHighPurityFlag"};
   TTreeReaderArray<Int_t> Muon_triggerIdLoose = {fReader, "Muon_triggerIdLoose"};
   TTreeReaderArray<UChar_t> Muon_pfIsoId = {fReader, "Muon_pfIsoId"};
   TTreeReaderArray<UChar_t> Muon_tkIsoId = {fReader, "Muon_tkIsoId"};
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
   TTreeReaderValue<Float_t> Pileup_nTrueInt = {fReader, "Pileup_nTrueInt"};
   TTreeReaderValue<Float_t> Pileup_pudensity = {fReader, "Pileup_pudensity"};
   TTreeReaderValue<Float_t> Pileup_gpudensity = {fReader, "Pileup_gpudensity"};
   TTreeReaderValue<Int_t> Pileup_nPU = {fReader, "Pileup_nPU"};
   TTreeReaderValue<Int_t> Pileup_sumEOOT = {fReader, "Pileup_sumEOOT"};
   TTreeReaderValue<Int_t> Pileup_sumLOOT = {fReader, "Pileup_sumLOOT"};
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
   TTreeReaderArray<Float_t> ProbeTracks_dzTrg = {fReader, "ProbeTracks_dzTrg"};
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
   TTreeReaderArray<Int_t> Muon_genPartIdx = {fReader, "Muon_genPartIdx"};
   TTreeReaderArray<Int_t> Muon_genPartFlav = {fReader, "Muon_genPartFlav"};
   TTreeReaderArray<Int_t> TriggerMuon_genPartIdx = {fReader, "TriggerMuon_genPartIdx"};
   TTreeReaderArray<Int_t> TriggerMuon_genPartFlav = {fReader, "TriggerMuon_genPartFlav"};
   TTreeReaderArray<Float_t> SV_chi2 = {fReader, "SV_chi2"};
   TTreeReaderArray<Float_t> SV_eta = {fReader, "SV_eta"};
   TTreeReaderArray<Float_t> SV_mass = {fReader, "SV_mass"};
   TTreeReaderArray<Float_t> SV_ndof = {fReader, "SV_ndof"};
   TTreeReaderArray<Float_t> SV_phi = {fReader, "SV_phi"};
   TTreeReaderArray<Float_t> SV_pt = {fReader, "SV_pt"};
   TTreeReaderArray<Float_t> SV_x = {fReader, "SV_x"};
   TTreeReaderArray<Float_t> SV_y = {fReader, "SV_y"};
   TTreeReaderArray<Float_t> SV_z = {fReader, "SV_z"};
   TTreeReaderArray<Int_t> ProbeTracks_genPartIdx = {fReader, "ProbeTracks_genPartIdx"};
   TTreeReaderArray<Int_t> ProbeTracks_genPartFlav = {fReader, "ProbeTracks_genPartFlav"};


   HNLToMuPiDumper(TTree * /*tree*/ =0) { }
   virtual ~HNLToMuPiDumper() { }
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

   // some branches will be added only if sample is MC
   Bool_t isMC;

   // trees to fill
   TTree* hnl_tree;
   
   // filling variables
   ULong64_t the_event= -99;
   Int_t the_run= -99;
   Int_t the_lumi= -99;

   Int_t the_pv_npvs= -99;

   Float_t the_hnl_pt= -99.;
   Float_t the_hnl_eta= -99.;
   Float_t the_hnl_phi= -99.;
   Float_t the_hnl_mass= -99.;
   Float_t the_hnl_mass_unfitted= -99.;
   Float_t the_diff_mass= -99.;
   Int_t the_hnl_charge= -99;
   Float_t the_hnl_cos2d= -99.;

   Float_t the_mu_pt= -99.;
   Float_t the_mu_eta= -99.;
   Float_t the_mu_phi= -99.;
   Int_t the_mu_charge= -99;
   Float_t the_mu_dxy= -99.;
   Float_t the_mu_dxysig= -99.;
   Float_t the_mu_dz= -99.;
   Float_t the_mu_dzsig= -99.;
   Int_t the_mu_looseid= -99;
   Int_t the_mu_mediumid= -99;
   Int_t the_mu_tightid= -99;
   Int_t the_mu_softid= -99;
   Int_t the_mu_istriggering= -99;
   Int_t the_mu_isslimmed= -99;
   Int_t the_mu_isdsa= -99;

   Float_t the_pi_pt= -99.;
   Float_t the_pi_eta= -99.;
   Float_t the_pi_phi= -99.;
   Int_t the_pi_charge= -99;
   Float_t the_pi_dcasig= -99.;
   Float_t the_pi_dxy= -99.;
   Float_t the_pi_dz= -99.;
   Float_t the_pi_dxysig= -99.;
   Float_t the_pi_dzsig= -99.;

   Float_t the_sv_chi2= -99.;
   Float_t the_sv_lxy= -99.;
   Float_t the_sv_lxysig= -99.;
   Float_t the_sv_prob= -99.;

   Int_t the_ismatched= -99;
   Int_t the_mu_ismatched= -99;
   Int_t the_pi_ismatched= -99;
   Float_t the_mupi_mass_reco_gen_reldiff= -99.;
   Float_t the_lxy_reco_gen_reldiff= -99.;

   // these have the gen information, independently on the matching
   Float_t the_gen_b_pt = -99.;
   Float_t the_gen_b_eta = -99.;
   Float_t the_gen_b_phi = -99.;
   Float_t the_gen_b_mass = -99.;
   Int_t the_gen_b_pdgid = -99;
   Float_t the_gen_hnl_ct = -99.;
   Float_t the_gen_hnl_lxy = -99.;
   Float_t the_gen_hnl_lxyz = -99.;
   Float_t the_gen_hnl_pt = -99.;
   Float_t the_gen_hnl_eta = -99.;
   Float_t the_gen_hnl_phi = -99.;
   Float_t the_gen_hnl_mass = -99.;
   Float_t the_gen_hnl_vx = -99.;
   Float_t the_gen_hnl_vy = -99.;
   Float_t the_gen_hnl_vz = -99.;
   Float_t the_gen_trgmu_pt = -99.;
   Float_t the_gen_trgmu_eta = -99.;
   Float_t the_gen_trgmu_phi = -99.;
   Float_t the_gen_trgmu_vx = -99.;
   Float_t the_gen_trgmu_vy = -99.;
   Float_t the_gen_trgmu_vz = -99.;
   Float_t the_gen_mu_pt = -99.;
   Float_t the_gen_mu_eta = -99.;
   Float_t the_gen_mu_phi = -99.;
   Float_t the_gen_mu_vx = -99.;
   Float_t the_gen_mu_vy = -99.;
   Float_t the_gen_mu_vz = -99.;
   Float_t the_gen_pi_pt = -99.;
   Float_t the_gen_pi_eta = -99.;
   Float_t the_gen_pi_phi = -99.;
   Float_t the_gen_pi_vx = -99.;
   Float_t the_gen_pi_vy = -99.;
   Float_t the_gen_pi_vz = -99.;

   ClassDef(HNLToMuPiDumper,0);

};

#endif

#ifdef HNLToMuPiDumper_cxx
void HNLToMuPiDumper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t HNLToMuPiDumper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef HNLToMuPiDumper_cxx
