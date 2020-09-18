import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import uint, ufloat, Var, CandVars

BToMuMuPi = cms.EDProducer(
    'BToMuMuPiBuilder',
    pionSelection           = cms.string('pt > 0.5 && abs(eta)<2.5'), 
    isoTracksSelection      = cms.string('pt > 0.5 && abs(eta)<2.5'),
    trgMuons                = cms.InputTag('muonTrgSelector', 'trgMuons'),
    selMuons                = cms.InputTag('muonTrgSelector', 'SelectedMuons'), 
    selMuonsTransientTracks = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'), 
    pions                   = cms.InputTag('tracksBPark', 'SelectedTracks'),
    pionsTransientTracks    = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks                  = cms.InputTag("packedPFCandidates"), 
    lostTracks              = cms.InputTag("lostTracks")        , 
    beamSpot                = cms.InputTag("offlineBeamSpot")   , 
#     preVtxSelection         = cms.string(''),
#     postVtxSelection        = cms.string(''),
    preVtxSelection = cms.string(' & '.join([
        'pt > 3.',
        'mass > 0.2',        
        'mass < 7.0',        
        ])
    ),
    postVtxSelection = cms.string(' & '.join([
        'userInt("hnl_vtx_OK") == 1',
        'userFloat("hnl_vtx_prob") > 0.0001',
        'userFloat("hnl_fitted_cos_theta_2D") >= 0.5',
        'userFloat("hnl_fitted_mass") > 0.5',
        'userFloat("hnl_fitted_mass") < 6.5',
        ])
    ),
)
    

BToMuMuPiTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToMuMuPi"),
    cut = cms.string(""),
    name = cms.string("b"),
    doc = cms.string("HNL Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        trg_mu_idx = uint('trg_mu_idx'),
        sel_mu_idx = uint('sel_mu_idx'),
        pi_idx     = uint('pi_idx'    ),
        # mu-pi fit and vtx info
        sv_chi2    = ufloat('hnl_vtx_chi2' ),
        sv_prob    = ufloat('hnl_vtx_prob' ),
        sv_lxy     = ufloat('hnl_l_xy'     ),
        sv_lxye    = ufloat('hnl_l_xy_unc' ),
        sv_lxy_sig = ufloat('hnl_ls_xy'    ),
        sv_x       = ufloat('hnl_vtx_x'    ),
        sv_y       = ufloat('hnl_vtx_y'    ),
        sv_z       = ufloat('hnl_vtx_z'    ),
        sv_xe      = ufloat('hnl_vtx_ex'   ), ## only saving diagonal elements of the cov matrix
        sv_ye      = ufloat('hnl_vtx_ey'   ),
        sv_ze      = ufloat('hnl_vtx_ez'   ),
        # HNL 
        hnl_mass   = Var('daughter("hnl").mass()', float),
        hnl_pt     = Var('daughter("hnl").pt()'  , float),
        hnl_eta    = Var('daughter("hnl").eta()' , float),
        hnl_phi    = Var('daughter("hnl").phi()' , float),
        hnl_charge = Var('daughter("hnl").charge()', int),
        # Cos(theta)
        hnl_cos2D     = ufloat('hnl_cos_theta_2D'       ),
        hnl_fit_cos2D = ufloat('hnl_fitted_cos_theta_2D'),
        # post-fit momentum
        hnl_fit_mass   = ufloat('hnl_fitted_mass'   ),
        hnl_fit_masse = ufloat('hnl_fitted_massErr'),
        hnl_fit_pt    = ufloat('hnl_fitted_pt'     ),
        hnl_fit_eta   = ufloat('hnl_fitted_eta'    ),
        hnl_fit_phi   = ufloat('hnl_fitted_phi'    ),
        fit_mu_pt     = ufloat('hnl_fitted_mu_pt'  ),
        fit_mu_eta    = ufloat('hnl_fitted_mu_eta' ),
        fit_mu_phi    = ufloat('hnl_fitted_mu_phi' ),
        fit_pi_pt     = ufloat('hnl_fitted_pi_pt'  ),
        fit_pi_eta    = ufloat('hnl_fitted_pi_eta' ),
        fit_pi_phi    = ufloat('hnl_fitted_pi_phi' ),
    )
)

CountBToMuMuPi = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuMuPi")
)    

BToMuMuPiSequence = cms.Sequence( BToMuMuPi )
