import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import uint, ufloat, Var, CandVars

BToMuMuPi = cms.EDProducer(
    'my_BToMuMuPiBuilder',
    trgMuons                = cms.InputTag('muonTrgSelector', 'trgMuons'),
    selMuons                = cms.InputTag('muonTrgSelector', 'SelectedMuons'), 
    selMuonsTransientTracks = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'), 
    pions                   = cms.InputTag('tracksBPark', 'SelectedTracks'),
    pionsTransientTracks    = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks                  = cms.InputTag("packedPFCandidates"), 
    lostTracks              = cms.InputTag("lostTracks"), 
    genParticles            = cms.InputTag("finalGenParticlesBPark"),
    beamSpot                = cms.InputTag('offlineBeamSpot'), 

    isMC = cms.bool(False),

    # pre-fitter preselection
    pionSelection           = cms.string('pt > 0.55 && abs(eta)<2'), # pion preselection to be modified here 
    isoTracksSelection      = cms.string('pt > 0.55 && abs(eta)<2'),
    trgMuonSelection        = cms.string('pt > 5 && abs(eta) < 1.7'),
    selMuonSelection        = cms.string('pt > 1.5 && abs(eta) < 2'),
    preVtxSelection = cms.string(' & '.join([
        'pt > 2',
        'mass > 0.2',        
        'mass < 7.0',        
        ])
    ), # applied on the HNL cand 

    # post-fitter preselection
    postVtxSelection = cms.string(' & '.join([
        'userInt("hnl_vtx_OK") == 1',
        'userFloat("hnl_vtx_prob") > 0.0001',
        'mass < 8',
        'pt > 11',
        'abs(eta) < 1.7',
        'userFloat("hnl_vtx_chi2") < 9',
        'userFloat("trg_muon_sip3d") > 0.8',
        ##'userFloat("sel_muon_ip3d") > 0.0015',
        ##'abs(userFloat("sel_muon_dxy")) > 0.0005',
        ##'abs(userFloat("pion_dz")) > 0.001',
        ##'abs(userFloat("pion_dxy")) > 0.0003',
        'abs(userFloat("pion_dzS")) > 0.5',
        'abs(userFloat("pion_dxyS")) > 0.2',
        'abs(userFloat("pion_DCASig")) > 0.3',
        ##'userFloat("muons_Lxyz") > 0.005',
        ##'userFloat("pion_muon_vzdiff") > 0.001',
        ##'userFloat("dr_mu_pi") < 1.8',
        'userFloat("dr_trgmu_hnl") < 0.8',
        'userFloat("hnl_fitted_mass") > 0.5',
        'userFloat("hnl_fitted_mass") < 6.5',
        'userFloat("hnl_fitted_pt") > 4',
        'abs(userFloat("hnl_fitted_eta")) < 1.8',
        'userFloat("hnl_fitted_cos_theta_2D") > 0.95',
        ])
    ), # applied on the B cand
)
    
BToMuMuPiMC = BToMuMuPi.clone(
    isMC = cms.bool(True),
)


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
        sel_mu_idx      = uint('sel_mu_idx'), 
        pi_idx          = uint('pi_idx'    ), 
        ## trigger muon
        trg_mu_pt       = ufloat('trg_muon_pt'   ), 
        trg_mu_eta      = ufloat('trg_muon_eta'  ), 
        trg_mu_phi      = ufloat('trg_muon_phi'  ), 
        ## vertex difference between the two muons
        dimu_vxdiff     = ufloat('dimuon_vxdiff' ),
        dimu_vydiff     = ufloat('dimuon_vydiff' ),
        dimu_vzdiff     = ufloat('dimuon_vzdiff' ),
        dimu_Lxy        = ufloat('dimuon_Lxy'    ),
        dimu_Lxyz       = ufloat('dimuon_Lxyz'   ),
        ## vertex difference between the trigger muon and pion
        pi_mu_vzdiff    = ufloat('pion_muon_vzdiff'  ),
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
        ## daughter muon
        fit_mu_pt       = ufloat('hnl_fitted_mu_pt'  ), 
        fit_mu_eta      = ufloat('hnl_fitted_mu_eta' ),
        fit_mu_phi      = ufloat('hnl_fitted_mu_phi' ),
        fit_mu_mass     = ufloat('hnl_fitted_mu_mass'),
        ## daughter pion
        fit_pi_pt       = ufloat('hnl_fitted_pi_pt'  ),
        fit_pi_eta      = ufloat('hnl_fitted_pi_eta' ),
        fit_pi_phi      = ufloat('hnl_fitted_pi_phi' ),
        fit_pi_mass     = ufloat('hnl_fitted_pi_mass'),
        ## dR quantities
        dr_mu_pi        = ufloat('dr_mu_pi'     ),
        dr_trgmu_hnl    = ufloat('dr_trgmu_hnl' ),
        # Other quantities
        ## ID WP of the selected muon 
        sel_mu_isSoft   = ufloat('sel_muon_isSoft'   ),
        sel_mu_isTight  = ufloat('sel_muon_isTight'  ),
        sel_mu_isMedium = ufloat('sel_muon_isMedium' ),
        sel_mu_isLoose  = ufloat('sel_muon_isLoose'  ),
        ## impact paramaters
        trg_mu_ip3d     = ufloat('trg_muon_ip3d'  ), 
        trg_mu_sip3d    = ufloat('trg_muon_sip3d' ), 
        trg_mu_dxy      = ufloat('trg_muon_dxy'   ), 
        trg_mu_dz       = ufloat('trg_muon_dz'    ), 
        sel_mu_ip3d     = ufloat('sel_muon_ip3d'  ), 
        sel_mu_sip3d    = ufloat('sel_muon_sip3d' ), 
        sel_mu_dxy      = ufloat('sel_muon_dxy'   ), 
        sel_mu_dz       = ufloat('sel_muon_dz'    ), 
        pi_dz           = ufloat('pion_dz'        ), 
        pi_dxy          = ufloat('pion_dxy'       ), 
        pi_dzS          = ufloat('pion_dzS'       ), 
        pi_dxyS         = ufloat('pion_dxyS'      ), 
        pi_DCASig       = ufloat('pion_DCASig'    ), 
        ## isolation
        trg_mu_iso03    = ufloat('trg_mu_iso03'   ), 
        trg_mu_iso04    = ufloat('trg_mu_iso04'   ),
        sel_mu_iso03    = ufloat('sel_mu_iso03'   ),
        sel_mu_iso04    = ufloat('sel_mu_iso04'   ),
        pi_iso03        = ufloat('pi_iso03'       ),
        pi_iso04        = ufloat('pi_iso04'       ),
        hnl_iso03       = ufloat('hnl_iso03'      ),
        hnl_iso04       = ufloat('hnl_iso04'      ),
        trg_mu_iso03_close  = ufloat('trg_mu_iso03_close' ), 
        trg_mu_iso04_close  = ufloat('trg_mu_iso04_close' ),
        sel_mu_iso03_close  = ufloat('sel_mu_iso03_close' ),
        sel_mu_iso04_close  = ufloat('sel_mu_iso04_close' ),
        pi_iso03_close      = ufloat('pi_iso03_close'     ),
        pi_iso04_close      = ufloat('pi_iso04_close'     ),
        hnl_iso03_close     = ufloat('hnl_iso03_close'    ),
        hnl_iso04_close     = ufloat('hnl_iso04_close'    ),
        ## dilepton mass
        dilepton_mass   = ufloat('dilepton_mass'  ),
        dilepton_pt     = ufloat('dilepton_pt'    ),
        ## gen-matching
        isMatched                   = Var("userInt('isMatched')"                  , int, mcOnly=True),
        matching_sel_mu_genIdx      = Var("userInt('matching_sel_mu_genIdx')"     , int, mcOnly=True),
        matching_trg_mu_genIdx      = Var("userInt('matching_trg_mu_genIdx')"     , int, mcOnly=True),
        matching_pi_genIdx          = Var("userInt('matching_pi_genIdx')"         , int, mcOnly=True),
        matching_sel_mu_motherPdgId = Var("userInt('matching_sel_mu_motherPdgId')", int, mcOnly=True),
        matching_trg_mu_motherPdgId = Var("userInt('matching_trg_mu_motherPdgId')", int, mcOnly=True),
        matching_pi_motherPdgId     = Var("userInt('matching_pi_motherPdgId')"    , int, mcOnly=True),
        ## displacement 
        #fitter_bs_lxy = ufloat('fitter_bs_lxy'), #same as sv_lxy, with more explicit naming
        ##my_fitter_bs_lxy = ufloat('my_fitter_bs_lxy'), 
        ##disp2DFromBS = ufloat('disp2DFromBS'),
        ##flightdisp = ufloat('flightdisp'),
        #gendisp_trgmu_mu_lxyz = ufloat('gendisp_trgmu_mu_lxyz'),
        #gendisp_trgmu_mu_lxy = ufloat('gendisp_trgmu_mu_lxy'),
    )
)


CountBToMuMuPi = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuMuPi")
)    

BToMuMuPiSequence   = cms.Sequence( BToMuMuPi )
BToMuMuPiSequenceMC = cms.Sequence( BToMuMuPiMC )
