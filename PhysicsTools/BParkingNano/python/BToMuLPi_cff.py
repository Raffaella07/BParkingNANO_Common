import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import uint, ufloat, Var, CandVars

BToMuMuPi = cms.EDProducer(
    'BToMuMuPiBuilder',
    trgMuons                = cms.InputTag('muonTrgSelector','SelectedMuons'),
    leptons                = cms.InputTag('muonTrgSelector', 'SelectedMuons'), 
    leptonsTransientTracks = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'), 
    pions                   = cms.InputTag('tracksBPark', 'SelectedTracks'),
    pionsTransientTracks    = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks                  = cms.InputTag("packedPFCandidates"), 
    lostTracks              = cms.InputTag("lostTracks"), 
    genParticles            = cms.InputTag("finalGenParticlesBPark"),
    beamSpot                = cms.InputTag('offlineBeamSpot'), 

    isMC = cms.bool(False),
    label = cms.string("muon"),

     pionSelection_loose = cms.string(' && '.join([
        'pt > 0.',
      ])
    ),
    pionSelection = cms.string(' && '.join([
        #'pt > 0.',
        'pt > 0.7',
        'abs(eta)<2.',
        'abs(userFloat("dz")) > 0.005',
        'abs(userFloat("dxy")) > 0.005',
        'abs(userFloat("dzS")) > 1.5',
        'abs(userFloat("dxyS")) > 3.',
        'abs(userFloat("DCASig")) > 5.',
      ])
    ),
    pionSelection_dsa = cms.string(' && '.join([
        #'pt > 0.',
        'pt > 0.6',
        'abs(eta)<2.',
        'abs(userFloat("dz")) > 0.005',
        'abs(userFloat("dxy")) > 0.001',
        'abs(userFloat("dzS")) > 1.5',
        'abs(userFloat("dxyS")) > 0.5',
        'abs(userFloat("DCASig")) > 0.5',
      ])
    ),
    #isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.'),
    isoTracksSelection = cms.string('pt > 0. && abs(eta)<5'),
    trgMuonSelection       = cms.string('pt > 7 && abs(eta) < 1.5 && softId==1'),
    #trgMuonSelection       = cms.string('pt > 7 && abs(eta) < 1.5 '),
    trgMuonSelection_dsa   = cms.string('pt > 7 && abs(eta) < 1.55'),
    trgMuonSelection_loose = cms.string('pt > 0 && abs(eta) < 5'),

    leptonSelection_loose = cms.string(' && '.join([
        'pt > 0.',
      ])
    ),
    leptonSelection = cms.string(' && '.join([
        #'pt > 0.',
        'pt > 1.5',
        'abs(eta) < 2.',
        'abs(userFloat("dz")) > 0.0015',
        'abs(userFloat("dxy")) > 0.001',
        'abs(userFloat("dzS")) > 1.',
        'abs(userFloat("dxyS")) > 1.5',
      ])
    ),
    leptonSelection_dsa = cms.string(' && '.join([
        #'pt > 0.',
        'pt > 2',
        'abs(eta) < 2.',
#       'abs(userFloat("dz")) > 0.0015',
#       'abs(userFloat("dxy")) > 0.1',
      ])
    ),

       preVtxSelection = cms.string(' & '.join([
        'pt > 1',
        'mass > 0.2',
        'mass < 7.0',
        ])
    ),

    # post-fitter preselection
    postVtxSelection_loose = cms.string(' & '.join([
        'userInt("hnl_vtx_OK") == 1',
#       'userFloat("hnl_fitted_cos_theta_2D") > 0.9',
#       'mass < 10',
        ])
    ),
    postVtxSelection = cms.string(' & '.join([
        #'userInt("hnl_vtx_OK") == 1',
        #'userFloat("hnl_fitted_cos_theta_2D") > 0.9',
        #'mass < 10',

        'userInt("hnl_vtx_OK") == 1',
        'userFloat("hnl_vtx_prob") > 0.001',
        'userFloat("hnl_fitted_cos_theta_2D") > 0.99',
        'userFloat("hnl_ls_xy") > 20',
        'mass < 8',
        'userFloat("hnl_fitted_mass")<6.3',
        ])
    ),
    postVtxSelection_dsa = cms.string(' & '.join([
        #'userInt("hnl_vtx_OK") == 1',
        #'userFloat("hnl_fitted_cos_theta_2D") > 0.9',
        #'mass < 10',

        'userInt("hnl_vtx_OK") == 1',
#       'userFloat("hnl_vtx_prob") > 0.001',
#       'userFloat("hnl_fitted_cos_theta_2D") > 0.95',
#       'mass < 8',
#       'userFloat("hnl_fitted_mass")<6.3',
        ])
    ),
)



#MC flag    
BToMuMuPiMC = BToMuMuPi.clone(
    isMC = cms.bool(True),
)

#muchannel table
BToMuMuPiTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToMuMuPi"),
    cut = cms.string(""),
    name = cms.string("BToMuMuPi"),
    doc = cms.string("HNL Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        trg_mu_idx      = uint('trg_mu_idx'),
        sel_mu_idx      = uint('lep_idx'),
        pi_idx          = uint('pi_idx'    ), 
        ## trigger muon
        trg_mu_pt       = ufloat('trg_muon_pt'   ), 
        trg_mu_eta      = ufloat('trg_muon_eta'  ), 
        trg_mu_phi      = ufloat('trg_muon_phi'  ), 
        ## vertex difference between the two muons
        dimu_vxdiff     = ufloat('dilepton_vxdiff' ),
        dimu_vydiff     = ufloat('dilepton_vydiff' ),
        dimu_vzdiff     = ufloat('dilepton_vzdiff' ),
        dimu_Lxy        = ufloat('dilepton_Lxy'    ),
        dimu_Lxyz       = ufloat('dilepton_Lxyz'   ),
        ## vertex difference between the trigger muon and pion
    #    pi_mu_vzdiff    = ufloat('pion_muon_vzdiff'  ),
        # post-fit quantities
        ## vertex information 
        sv_chi2         = ufloat('hnl_vtx_chi2' ),
        sv_prob         = ufloat('hnl_vtx_prob' ),
        sv_lxy          = ufloat('hnl_l_xy'     ),
        sv_lxye         = ufloat('hnl_l_xy_unc' ),
        sv_lxy_sig      = ufloat('hnl_ls_xy'    ),
        sv_x            = ufloat('hnl_vtx_x'    ),
        sv_y            = ufloat('hnl_vtx_y'    ),
        sv_z            = ufloat('hnl_vtx_z'    ),
        sv_xe           = ufloat('hnl_vtx_ex'   ), ## only saving diagonal elements of the cov matrix
        sv_ye           = ufloat('hnl_vtx_ey'   ),
        sv_ze           = ufloat('hnl_vtx_ez'   ),
        ## HNL (only postfit information is saved) 
        hnl_mass        = ufloat('hnl_fitted_mass'   ),
        hnl_masserr     = ufloat('hnl_fitted_massErr'),
        hnl_pt          = ufloat('hnl_fitted_pt'     ),
        hnl_eta         = ufloat('hnl_fitted_eta'    ),
        hnl_phi         = ufloat('hnl_fitted_phi'    ),
        hnl_charge      = Var('daughter("hnl").charge()', int),
        hnl_cos2D       = ufloat('hnl_fitted_cos_theta_2D'   ),
        hnl_cos2D_star  = ufloat('hnl_fitted_cos_theta_2D_star'   ),
        ## daughter muon
        fit_mu_pt       = ufloat('hnl_fitted_lep_pt'  ), 
        fit_mu_eta      = ufloat('hnl_fitted_lep_eta' ),
        fit_mu_phi      = ufloat('hnl_fitted_lep_phi' ),
        fit_mu_mass     = ufloat('hnl_fitted_lep_mass'),
        ## daughter pion
        fit_pi_pt       = ufloat('hnl_fitted_pi_pt'  ),
        fit_pi_eta      = ufloat('hnl_fitted_pi_eta' ),
        fit_pi_phi      = ufloat('hnl_fitted_pi_phi' ),
        fit_pi_mass     = ufloat('hnl_fitted_pi_mass'),
        ## dR quantities
        dr_mu_pi        = ufloat('dr_lep_pi'     ),
        dr_trgmu_hnl    = ufloat('dr_trgmu_hnl' ),
        # Other quantities
        ## ID WP of the selected muon 
  #     sel_mu_isSoft   = ufloat('sel_muon_isSoft'   ),
  #     sel_mu_isTight  = ufloat('sel_muon_isTight'  ),
  #     sel_mu_isMedium = ufloat('sel_muon_isMedium' ),
  #     sel_mu_isLoose  = ufloat('sel_muon_isLoose'  ),
        ## impact paramaters
   #     trg_mu_ip3d     = ufloat('trg_muon_ip3d'  ), 
    #    trg_mu_sip3d    = ufloat('trg_muon_sip3d' ), 
        trg_mu_dxy      = ufloat('trg_muon_dxy'   ), 
        trg_mu_dz       = ufloat('trg_muon_dz'    ), 
  #      sel_mu_ip3d     = ufloat('sel_lep_ip3d'  ), 
   #     sel_mu_sip3d    = ufloat('sel_lep_sip3d' ), 
        sel_mu_dxy      = ufloat('sel_lep_dxy'   ), 
        sel_mu_dz       = ufloat('sel_lep_dz'    ), 
        pi_dz           = ufloat('pion_dz'        ), 
        pi_dxy          = ufloat('pion_dxy'       ), 
        pi_dzS          = ufloat('pion_dzS'       ), 
        pi_dxyS         = ufloat('pion_dxyS'      ), 
        pi_DCASig       = ufloat('pion_DCASig'    ), 
        ## isolation
 #      trg_mu_iso03    = ufloat('trg_mu_iso03'   ), 
 #      trg_mu_iso04    = ufloat('trg_mu_iso04'   ),
 #      sel_mu_iso03    = ufloat('sel_lep_iso03'   ),
 #      sel_mu_iso04    = ufloat('sel_lep_iso04'   ),
 #      pi_iso03        = ufloat('pi_iso03'       ),
 #      pi_iso04        = ufloat('pi_iso04'       ),
 #      hnl_iso03       = ufloat('hnl_iso03'      ),
 #      hnl_iso04       = ufloat('hnl_iso04'      ),
 #  #   trg_mu_iso03_close  = ufloat('trg_mu_iso03_close' ), 
    #   trg_mu_iso04_close  = ufloat('trg_mu_iso04_close' ),
    #   sel_mu_iso03_close  = ufloat('sel_lep_iso03_close' ),
    #   sel_mu_iso04_close  = ufloat('sel_lep_iso04_close' ),
    #   pi_iso03_close      = ufloat('pi_iso03_close'     ),
    #   pi_iso04_close      = ufloat('pi_iso04_close'     ),
    #   hnl_iso03_close     = ufloat('hnl_iso03_close'    ),
    #   hnl_iso04_close     = ufloat('hnl_iso04_close'    ),
    #   ## dilepton mass
        dilepton_mass   = ufloat('dilepton_mass'  ),
        dilepton_pt     = ufloat('dilepton_pt'    ),
        ## gen-matching
        isMatched                   = Var("userInt('isMatched')"                  , int, mcOnly=True),
        matching_sel_mu_genIdx      = Var("userInt('matching_sel_lep_genIdx')"     , int, mcOnly=True),
        matching_trg_mu_genIdx      = Var("userInt('matching_trg_mu_genIdx')"     , int, mcOnly=True),
        matching_pi_genIdx          = Var("userInt('matching_pi_genIdx')"         , int, mcOnly=True),
        matching_sel_mu_motherPdgId = Var("userInt('matching_sel_lep_motherPdgId')", int, mcOnly=True),
        matching_trg_mu_motherPdgId = Var("userInt('matching_trg_mu_motherPdgId')", int, mcOnly=True),
        matching_pi_motherPdgId     = Var("userInt('matching_pi_motherPdgId')"    , int, mcOnly=True),
        ## displacement 
        pimu_mass_reco_gen_reldiff     = Var("userFloat('pilep_mass_reco_gen_reldiff')"    , float, mcOnly=True),
        lxy_reco_gen_reldiff     = Var("userFloat('lxy_reco_gen_reldiff')"    ,float, mcOnly=True),
        #fitter_bs_lxy = ufloat('fitter_bs_lxy'), #same as sv_lxy, with more explicit naming
        ##my_fitter_bs_lxy = ufloat('my_fitter_bs_lxy'), 
        ##disp2DFromBS = ufloat('disp2DFromBS'),
        ##flightdisp = ufloat('flightdisp'),
        #gendisp_trgmu_mu_lxyz = ufloat('gendisp_trgmu_mu_lxyz'),
        #gendisp_trgmu_mu_lxy = ufloat('gendisp_trgmu_mu_lxy'),
    )
)

#count Mu
BToMuEPi = cms.EDProducer(
    'BToMuEPiBuilder',
    trgMuons                = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    leptons                 = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'), 
    leptonsTransientTracks  = cms.InputTag('electronsForAnalysis','SelectedTransientElectrons'), 
    pions                   = cms.InputTag('tracksBPark', 'SelectedTracks'),
    pionsTransientTracks    = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks                  = cms.InputTag("packedPFCandidates"), 
    lostTracks              = cms.InputTag("lostTracks"), 
    genParticles            = cms.InputTag("finalGenParticlesBPark"),
    beamSpot                = cms.InputTag('offlineBeamSpot'), 

    isMC = cms.bool(False),
    label = cms.string("ele"),
 # pre-fitter preselection
   # pionSelection           = cms.string('pt > 0.7 && abs(eta)<2'), # pion preselection to be modified here 
    isoTracksSelection      = cms.string('pt > 0.7 && abs(eta)<2'),
    #trgMuonSelection        = cms.string('pt > 7 && abs(eta) < 1.5 '),
    trgMuonSelection        = cms.string('pt > 7 && abs(eta) < 1.5 && softId==1'),
   # leptonSelection        = cms.string('pt > 1.5 && abs(eta) < 2'),
    pionSelection        = cms.string(' & '.join([
	'pt>0.7',
	'abs(eta)<2.5',
        'abs(userFloat("dz")) > 0.005',
        'abs(userFloat("dxy")) > 0.005',
        'abs(userFloat("dzS")) > 1.5',
        'abs(userFloat("dxyS")) > 3',
	])
	),		
    leptonSelection        = cms.string(' & '.join([
	'abs(eta)<2',
        'abs(userFloat("pfmvaId")) > -3',
        'abs(userFloat("dz")) > 0.008',
        'abs(userFloat("dxy")) > 0.008',
        'abs(userFloat("dz")/abs(edB("PVDZ"))) > 1',
        'abs(userFloat("dxy")/abs(edB("PV2D"))) > 1.5',
	])
	),		
    preVtxSelection = cms.string(' & '.join([
    #    'pt > 2',
        'mass > 0.2',        
        'mass < 7.0',        
        ])
    ), # applied on the HNL cand 

    # post-fitter preselection
    postVtxSelection = cms.string(' & '.join([
        'userInt("hnl_vtx_OK") == 1',
        'userFloat("hnl_vtx_prob") > 0.001',
        'userFloat("hnl_fitted_cos_theta_2D") > 0.99',
        'userFloat("hnl_ls_xy") > 3',
        'mass < 8',
        'userFloat("hnl_fitted_mass")<6.3',
        ])
    ), # applied on the B cand
)
#MC flag    

#MC flag    

BToMuEPiMC = BToMuEPi.clone(
    isMC = cms.bool(True),
)

#muchannel table
BToMuEPiTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToMuEPi"),
    cut = cms.string(""),
    name = cms.string("BToMuEPi"),
    doc = cms.string("HNL Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        trg_mu_idx      = uint('trg_mu_idx'),
        sel_ele_idx      = uint('lep_idx'), 
        pi_idx          = uint('pi_idx'    ), 
        ## trigger muon
        trg_mu_pt       = ufloat('trg_muon_pt'   ), 
        trg_mu_eta      = ufloat('trg_muon_eta'  ), 
        trg_mu_phi      = ufloat('trg_muon_phi'  ), 
        ## vertex difference between the two muons
        dilep_vxdiff     = ufloat('dilepton_vxdiff' ),
        dilep_vydiff     = ufloat('dilepton_vydiff' ),
        dilep_vzdiff     = ufloat('dilepton_vzdiff' ),
        dilep_Lxy        = ufloat('dilepton_Lxy'    ),
        dilep_Lxyz       = ufloat('dilepton_Lxyz'   ),
       ### vertex difference between the trigger muon and pion
        #pi_mu_vzdiff    = ufloat('pion_muon_vzdiff'  ),
        # post-fit quantities
        ## vertex information 
        sv_chi2         = ufloat('hnl_vtx_chi2' ),
        sv_prob         = ufloat('hnl_vtx_prob' ),
        sv_lxy          = ufloat('hnl_l_xy'     ),
        sv_lxye         = ufloat('hnl_l_xy_unc' ),
        sv_lxy_sig      = ufloat('hnl_ls_xy'    ),
        sv_x            = ufloat('hnl_vtx_x'    ),
        sv_y            = ufloat('hnl_vtx_y'    ),
        sv_z            = ufloat('hnl_vtx_z'    ),
        sv_xe           = ufloat('hnl_vtx_ex'   ), ## only saving diagonal elements of the cov matrix
        sv_ye           = ufloat('hnl_vtx_ey'   ),
        sv_ze           = ufloat('hnl_vtx_ez'   ),
        ## HNL (only postfit information is saved) 
        hnl_mass        = ufloat('hnl_fitted_mass'   ),
        hnl_masserr     = ufloat('hnl_fitted_massErr'),
        hnl_pt          = ufloat('hnl_fitted_pt'     ),
        hnl_eta         = ufloat('hnl_fitted_eta'    ),
        hnl_phi         = ufloat('hnl_fitted_phi'    ),
        hnl_charge      = Var('daughter("hnl").charge()', int),
        hnl_cos2D       = ufloat('hnl_fitted_cos_theta_2D'   ),
        hnl_cos2D_star  = ufloat('hnl_fitted_cos_theta_2D_star'),
        ## daughter muon
        fit_e_pt       = ufloat('hnl_fitted_lep_pt'  ), 
        fit_e_eta      = ufloat('hnl_fitted_lep_eta' ),
        fit_e_phi      = ufloat('hnl_fitted_lep_phi' ),
        fit_e_mass     = ufloat('hnl_fitted_lep_mass'),
        ## daughter pion
        fit_pi_pt       = ufloat('hnl_fitted_pi_pt'  ),
        fit_pi_eta      = ufloat('hnl_fitted_pi_eta' ),
        fit_pi_phi      = ufloat('hnl_fitted_pi_phi' ),
        fit_pi_mass     = ufloat('hnl_fitted_pi_mass'),
        ## dR quantities
        dr_ele_pi        = ufloat('dr_lep_pi'     ),
        dr_trgmu_hnl    = ufloat('dr_trgmu_hnl' ),
        # Other quantities
        ## ID WP of the selected ele 
        ## impact paramaters
      #  trg_mu_ip3d     = ufloat('trg_muon_ip3d'  ), 
       # trg_mu_sip3d    = ufloat('trg_muon_sip3d' ), 
        trg_mu_dxy      = ufloat('trg_muon_dxy'   ), 
        trg_mu_dz       = ufloat('trg_muon_dz'    ), 
       # sel_ele_ip3d     = ufloat('sel_lep_ip3d'  ), 
      #  sel_ele_sip3d    = ufloat('sel_lep_sip3d' ), 
        sel_ele_dxy      = ufloat('sel_lep_dxy'   ), 
        sel_ele_dz       = ufloat('sel_lep_dz'    ), 
        pi_dz           = ufloat('pion_dz'        ), 
        pi_dxy          = ufloat('pion_dxy'       ), 
        pi_dzS          = ufloat('pion_dzS'       ), 
        pi_dxyS         = ufloat('pion_dxyS'      ), 
        pi_DCASig       = ufloat('pion_DCASig'    ), 
        ## isolation
     #  trg_mu_iso03    = ufloat('trg_mu_iso03'   ), 
     #  trg_mu_iso04    = ufloat('trg_mu_iso04'   ),
     #  sel_mu_iso03    = ufloat('sel_lep_iso03'   ),
     #  sel_mu_iso04    = ufloat('sel_lep_iso04'   ),
     #  pi_iso03        = ufloat('pi_iso03'       ),
     #  pi_iso04        = ufloat('pi_iso04'       ),
     #  hnl_iso03       = ufloat('hnl_iso03'      ),
     #  hnl_iso04       = ufloat('hnl_iso04'      ),
     #  trg_mu_iso03_close  = ufloat('trg_mu_iso03_close' ), 
     #  trg_mu_iso04_close  = ufloat('trg_mu_iso04_close' ),
     #  sel_ele_iso03_close  = ufloat('sel_lep_iso03_close' ),
     #  sel_ele_iso04_close  = ufloat('sel_lep_iso04_close' ),
     #  pi_iso03_close      = ufloat('pi_iso03_close'     ),
     #  pi_iso04_close      = ufloat('pi_iso04_close'     ),
     #  hnl_iso03_close     = ufloat('hnl_iso03_close'    ),
     #  hnl_iso04_close     = ufloat('hnl_iso04_close'    ),
     #  ## dilepton mass
        dilepton_mass   = ufloat('dilepton_mass'  ),
        dilepton_pt     = ufloat('dilepton_pt'    ),
        ## gen-matching
        isMatched                   = Var("userInt('isMatched')"                  , int, mcOnly=True),
        matching_sel_lep_genIdx      = Var("userInt('matching_sel_lep_genIdx')"     , int, mcOnly=True),
        matching_trg_mu_genIdx      = Var("userInt('matching_trg_mu_genIdx')"     , int, mcOnly=True),
        matching_pi_genIdx          = Var("userInt('matching_pi_genIdx')"         , int, mcOnly=True),
        matching_sel_ele_motherPdgId = Var("userInt('matching_sel_lep_motherPdgId')", int, mcOnly=True),
        matching_trg_mu_motherPdgId = Var("userInt('matching_trg_mu_motherPdgId')", int, mcOnly=True),
        matching_pi_motherPdgId     = Var("userInt('matching_pi_motherPdgId')"    , int, mcOnly=True),
        piele_mass_reco_gen_reldiff     = Var("userFloat('pilep_mass_reco_gen_reldiff')"    , float, mcOnly=True),
        lxy_reco_gen_reldiff     = Var("userFloat('lxy_reco_gen_reldiff')"    , float, mcOnly=True),

        ## displacement 
        #fitter_bs_lxy = ufloat('fitter_bs_lxy'), #same as sv_lxy, with more explicit naming
        ##my_fitter_bs_lxy = ufloat('my_fitter_bs_lxy'), 
        ##disp2DFromBS = ufloat('disp2DFromBS'),
        ##flightdisp = ufloat('flightdisp'),
        #gendisp_trgmu_mu_lxyz = ufloat('gendisp_trgmu_mu_lxyz'),
        #gendisp_trgmu_mu_lxy = ufloat('gendisp_trgmu_mu_lxy'),
    )
)

BToMuEPiMCTable = BToMuEPiTable.clone(
	src = cms.InputTag("BToMuEPiMC"),
        isMC = cms.bool(True),	
)

BToMuMuPiMCTable = BToMuMuPiTable.clone(
	src = cms.InputTag("BToMuMuPiMC"),
        isMC = cms.bool(True),	
)

CountBToMuMuPi = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuMuPi")
)    
CountBToMuMuPiMC = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuMuPiMC")
)    
CountBToMuEPi = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuEPi")
)    
CountBToMuEPiMC = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuEPiMC")
)    
#muon sequences
BToMuMuPiSequence   = cms.Sequence( BToMuMuPi )
BToMuMuPiSequenceMC = cms.Sequence( BToMuMuPiMC )
#count Mu
#muon sequences
BToMuEPiSequence   = cms.Sequence( BToMuEPi )
BToMuEPiSequenceMC = cms.Sequence( BToMuEPiMC )

