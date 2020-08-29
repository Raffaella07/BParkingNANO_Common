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
#         'userFloat("hnl_vtx_prob") > 0.00001',
#         'userFloat("hnl_fitted_cos_theta_2D") >= 0',
#         'userFloat("hnl_fitted_mass") > 0.5',
#         'userFloat("hnl_fitted_mass") < 6.5',
        ])
    ),
)
    

BToMuMuPiTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToMuMuPi"),
    cut = cms.string(""),
    name = cms.string("BToMuMuPi"),
    doc = cms.string("BToMuMuPi Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
#         l1Idx    = uint('l1_idx'  ),
#         l2Idx    = uint('l2_idx'  ),
#         kIdx     = uint('k_idx'   ),
#         minDR    = ufloat('min_dr'),
#         maxDR    = ufloat('max_dr'),
        # mu-pi fit and vtx info
        chi2     = ufloat('hnl_vtx_chi2' ),
        svprob   = ufloat('hnl_vtx_prob' ),
        l_xy     = ufloat('hnl_l_xy'     ),
        l_xy_unc = ufloat('hnl_l_xy_unc' ),
        vtx_x    = ufloat('hnl_vtx_x'    ),
        vtx_y    = ufloat('hnl_vtx_y'    ),
        vtx_z    = ufloat('hnl_vtx_z'    ),
        vtx_ex   = ufloat('hnl_vtx_ex'   ), ## only saving diagonal elements of the cov matrix
        vtx_ey   = ufloat('hnl_vtx_ey'   ),
        vtx_ez   = ufloat('hnl_vtx_ez'   ),
        # HNL mass
#         hnl_m_raw     = Var('userCand("hnl").mass()'                     , float),
#         hnl_m_fit     = Var('userCand("hnl").userFloat("fitted_mass")'   , float), # this might not work
#         hnl_m_err_fit = Var('userCand("hnl").userFloat("fitted_massErr")', float), # this might not work
        # Cos(theta)
        cos2D     = ufloat('hnl_cos_theta_2D'       ),
        fit_cos2D = ufloat('hnl_fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass    = ufloat('hnl_fitted_mass'   ),
        fit_massErr = ufloat('hnl_fitted_massErr'),
        fit_pt      = ufloat('hnl_fitted_pt'     ),
        fit_eta     = ufloat('hnl_fitted_eta'    ),
        fit_phi     = ufloat('hnl_fitted_phi'    ),
        fit_mu_pt   = ufloat('hnl_fitted_mu_pt'  ),
        fit_mu_eta  = ufloat('hnl_fitted_mu_eta' ),
        fit_mu_phi  = ufloat('hnl_fitted_mu_phi' ),
        fit_pi_pt   = ufloat('hnl_fitted_pi_pt'  ),
        fit_pi_eta  = ufloat('hnl_fitted_pi_eta' ),
        fit_pi_phi  = ufloat('hnl_fitted_pi_phi' ),
    )
)

CountBToMuMuPi = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuMuPi")
)    

BToMuMuPiSequence = cms.Sequence( BToMuMuPi )
