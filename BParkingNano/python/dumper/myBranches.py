###
# This file consists in listing the different branches that we want to
# write to the ntupliser
#
# format : ('name_in_flat_file', 'name_in_nano_file')
#
# convention: no capital letters in name_in_flat_file
###

the_common_branches = [
   ('ismatched', 'isMatched'),
]

the_signal_branches = [
   ('b_pt', 'pt'),
   ('b_eta', 'eta'),
   ('b_phi', 'phi'),
   ('b_mass', 'mass'),
   ('b_charge', 'charge'),
   ('b_pdgid', 'pdgId'  ),

   ('hnl_pt', 'hnl_pt'),
   ('hnl_eta', 'hnl_eta'),
   ('hnl_phi', 'hnl_phi'),
   ('hnl_mass', 'hnl_mass'),
   ('hnl_charge', 'hnl_charge'),
   ('hnl_cos2d', 'hnl_cos2D'),
   ('hnl_iso03', 'hnl_iso03'),
   ('hnl_iso03_close', 'hnl_iso03_close'),
   ('hnl_iso04', 'hnl_iso04'),
   ('hnl_iso04_close', 'hnl_iso04_close'),

   ('trgmu_pt', 'trg_mu_pt'),
   ('trgmu_eta', 'trg_mu_eta'),
   ('trgmu_phi', 'trg_mu_phi'),
   ('trgmu_dxy', 'trg_mu_dxy'),
   ('trgmu_dz', 'trg_mu_dz'),
   ('trgmu_ip3d', 'trg_mu_ip3d'),
   ('trgmu_ip3dsig', 'trg_mu_sip3d'),
   ('trgmu_iso03', 'trg_mu_iso03'),
   ('trgmu_iso03_close', 'trg_mu_iso03_close'),
   ('trgmu_iso04', 'trg_mu_iso04'),
   ('trgmu_iso04_close', 'trg_mu_iso04_close'),

   ('mu_pt', 'fit_mu_pt'),
   ('mu_eta', 'fit_mu_eta'),
   ('mu_phi', 'fit_mu_phi'),
   ('mu_dxy', 'sel_mu_dxy'),
   ('mu_dz', 'sel_mu_dz'),
   ('mu_ip3d', 'sel_mu_ip3d'),
   ('mu_ip3dsig', 'sel_mu_sip3d'),
   ('mu_iso03', 'sel_mu_iso03'),
   ('mu_iso03_close', 'sel_mu_iso03_close'),
   ('mu_iso04', 'sel_mu_iso04'),
   ('mu_iso04_close', 'sel_mu_iso04_close'),
   ('mu_isloose', 'sel_mu_isLoose'),
   ('mu_ismedium', 'sel_mu_isMedium'),
   ('mu_istight', 'sel_mu_isTight'),
   ('mu_issoft', 'sel_mu_isSoft'),

   ('pi_pt', 'fit_pi_pt'),
   ('pi_eta', 'fit_pi_eta'),
   ('pi_phi', 'fit_pi_phi'),
   ('pi_dcasig', 'pi_DCASig'),
   ('pi_dxy', 'pi_dxy'),
   ('pi_dz', 'pi_dz'),
   ('pi_dxysig', 'pi_dxyS'),
   ('pi_dzsig', 'pi_dzS'),
   ('pi_iso03', 'pi_iso03'),
   ('pi_iso03_close', 'pi_iso03_close'),
   ('pi_iso04', 'pi_iso04'),
   ('pi_iso04_close', 'pi_iso04_close'),

   ('dimu_mass', 'dilepton_mass'),
   ('dimu_pt', 'dilepton_pt'),
   ('dimu_lxy', 'dimu_Lxy'),
   ('dimu_lxyz', 'dimu_Lxyz'),
   ('dimu_vxdiff', 'dimu_vxdiff'),
   ('dimu_vydiff', 'dimu_vydiff'),
   ('dimu_vzdiff', 'dimu_vzdiff'),

   ('deltar_mu_pi', 'dr_mu_pi'),
   ('deltar_trgmu_hnl', 'dr_trgmu_hnl'),

   ('sv_chi2', 'sv_chi2'),
   ('sv_lxy', 'sv_lxy'),
   ('sv_lxysig', 'sv_lxy_sig'),
   ('sv_prob', 'sv_prob'),
   ('sv_x', 'sv_x'),
   ('sv_y', 'sv_y'),
   ('sv_z', 'sv_z'),

   ('pi_mu_vzdiff', 'pi_mu_vzdiff'),
]

the_control_branches = [
   ('b_pt'    , 'fit_pt'    ),
   ('b_eta'   , 'fit_eta'   ),
   ('b_phi'   , 'fit_phi'   ),
   ('b_mass'  , 'fit_mass'  ),
   ('b_charge'  , 'charge'  ),
   ('b_pdgid'  , 'pdgId'  ),
   ('b_cos2d', 'fit_cos2D'),
   ('b_iso03', 'b_iso03'),
   ('b_iso03_close', 'b_iso03_close'),
   ('b_iso04', 'b_iso04'),
   ('b_iso04_close', 'b_iso04_close'),

   ('k_pt', 'fit_k_pt'),
   ('k_eta', 'fit_k_eta'),
   ('k_phi', 'fit_k_phi'),
   ('k_iso03', 'k_iso03'),
   ('k_iso03_close', 'k_iso03_close'),
   ('k_iso04', 'k_iso04'),
   ('k_iso04_close', 'k_iso04_close'),

   ('l1_pt', 'fit_l1_pt'),
   ('l1_eta', 'fit_l1_eta'),
   ('l1_phi', 'fit_l1_phi'),
   ('l1_iso03', 'l1_iso03'),
   ('l1_iso03_close', 'l1_iso03_close'),
   ('l1_iso04', 'l1_iso04'),
   ('l1_iso04_close', 'l1_iso04_close'),

   ('l2_pt', 'fit_l2_pt'),
   ('l2_eta', 'fit_l2_eta'),
   ('l2_phi', 'fit_l2_phi'),
   ('l2_iso03', 'l2_iso03'),
   ('l2_iso03_close', 'l2_iso03_close'),
   ('l2_iso04', 'l2_iso04'),
   ('l2_iso04_close', 'l2_iso04_close'),

   ('dimu_mass', 'mll_fullfit'),

   ('vtx_x', 'vtx_x'),
   ('vtx_y', 'vtx_y'),
   ('vtx_z', 'vtx_z'),
   ('vtx_lxy', 'l_xy'),
   ('sv_prob', 'svprob'),
]

# branchname will not be prepended
the_extra_branches = [
   ('hlt_mu7_ip4', 'HLT_Mu7_IP4'), # custom column cannot have the same name as branch in original tree
   ('hlt_mu8_ip6', 'HLT_Mu8_IP6'), 
   ('hlt_mu8_ip5', 'HLT_Mu8_IP5'), 
   ('hlt_mu8_ip3', 'HLT_Mu8_IP3'), 
   ('hlt_mu8p5_ip3p5', 'HLT_Mu8p5_IP3p5'), 
   ('hlt_mu9_ip6', 'HLT_Mu9_IP6'), 
   ('hlt_mu9_ip5', 'HLT_Mu9_IP5'), 
   ('hlt_mu9_ip4', 'HLT_Mu9_IP4'), 
   ('hlt_mu10p5_ip3p5', 'HLT_Mu10p5_IP3p5'), 
   ('hlt_mu12_ip6', 'HLT_Mu12_IP6'), 
]

