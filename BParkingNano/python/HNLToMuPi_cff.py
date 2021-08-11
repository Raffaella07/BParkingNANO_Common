import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import uint, ufloat, Var, CandVars

HNLToMuPi = cms.EDProducer(
    'HNLToMuPiBuilder',
    leptons                 = cms.InputTag('muonTrgSelector', 'SelectedMuons'), 
    leptonsTransientTracks  = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'), 
    pions                   = cms.InputTag('tracksBPark', 'SelectedTracks'),
    pionsTransientTracks    = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    genParticles            = cms.InputTag("finalGenParticlesBPark"),
    beamSpot                = cms.InputTag('offlineBeamSpot'), 

    label = cms.string('muon'),
    isMC  = cms.bool(False),

    # pre-fitter preselection
    pionSelection = cms.string(' && '.join([
        'pt > 0.',
        #'pt > 0.7',
        #'abs(eta)<2.',
        #'abs(userFloat("dz")) > 0.005',
        #'abs(userFloat("dxy")) > 0.005',
        #'abs(userFloat("dzS")) > 1.5',  
        #'abs(userFloat("dxyS")) > 3.',
        #'abs(userFloat("DCASig")) > 5.',
      ])
    ),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.'),
    #isoTracksSelection = cms.string('pt > 0. && abs(eta)<5'),
    leptonSelection    = cms.string(' && '.join([
        'pt > 0.',
        #'pt > 1.5',
        #'abs(eta) < 2.',
        #'abs(userFloat("dz")) > 0.0015',
        #'abs(userFloat("dxy")) > 0.001',
        #'abs(userFloat("dzS")) > 1.',
        #'abs(userFloat("dxyS")) > 1.5',
        ##'abs(userFloat("sip3d")) > 7',
      ])
    ),
    preVtxSelection = cms.string(' & '.join([
        #'pt > 1',
        #'mass > 0.2',        
        #'mass < 7.0',        
        ])
    ),  

    # post-fitter preselection
    postVtxSelection = cms.string(' & '.join([
        'userInt("hnl_vtx_OK") == 1',
        'userFloat("hnl_fitted_cos_theta_2D") > 0.9',
        'mass < 10',

        #'userInt("hnl_vtx_OK") == 1',
        ##'abs(userFloat("sel_muon_dz")) > 0.0015',
        ##'abs(userFloat("sel_muon_dxy")) > 0.001',
        ##'userFloat("sel_muon_sip3d") > 7',
        ##'abs(userFloat("pion_dz")) > 0.005',
        ##'abs(userFloat("pion_dzS")) > 1.5',  
        ##'abs(userFloat("pion_dxy")) > 0.005',
        ##'abs(userFloat("pion_dxyS")) > 3',
        ##'abs(userFloat("pion_DCASig")) > 5',
        #'userFloat("hnl_vtx_prob") > 0.001',
        #'userFloat("hnl_fitted_cos_theta_2D") > 0.99',
        #'userFloat("hnl_ls_xy") > 20',
        #'mass < 8',
        #'userFloat("hnl_fitted_mass")<6.3',
        ##'abs(userFloat("deta_pi_fit_pi")) < 0.015',
        ##'abs(userFloat("dphi_pi_fit_pi")) < 0.03',
        ##'userFloat("dr_trgmu_hnl") < 0.5',
        ])
    ), 
)
    
HNLToMuPiMC = HNLToMuPi.clone(
    isMC = cms.bool(True),
)


HNLToMuPiTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("HNLToMuPi"),
    cut = cms.string(""),
    name = cms.string("HNLToMuPi"),
    doc = cms.string("HNL Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        sel_mu_idx      = uint('lep_idx'), 
        pi_idx          = uint('pi_idx'    ), 
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
        hnl_mass_unfitted = ufloat('hnl_mass_unfitted' ),
        hnl_masserr     = ufloat('hnl_fitted_massErr'),
        hnl_pt          = ufloat('hnl_fitted_pt'     ),
        hnl_eta         = ufloat('hnl_fitted_eta'    ),
        hnl_phi         = ufloat('hnl_fitted_phi'    ),
        hnl_charge      = uint('hnl_charge'    ),
        hnl_cos2D       = ufloat('hnl_fitted_cos_theta_2D'   ),
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
        # Other quantities
        ## gen-matching
        isMatched                   = Var("userInt('isMatched')"                  , int, mcOnly=True),
        sel_mu_isMatched            = Var("userInt('sel_mu_isMatched')"           , int, mcOnly=True),
        pi_isMatched                = Var("userInt('pi_isMatched')"               , int, mcOnly=True),
        matching_sel_mu_genIdx      = Var("userInt('matching_sel_mu_genIdx')"     , int, mcOnly=True),
        matching_pi_genIdx          = Var("userInt('matching_pi_genIdx')"         , int, mcOnly=True),
        matching_sel_mu_motherPdgId = Var("userInt('matching_sel_mu_motherPdgId')", int, mcOnly=True),
        matching_pi_motherPdgId     = Var("userInt('matching_pi_motherPdgId')"    , int, mcOnly=True),
        ## reco/gen relative difference
        mupi_mass_reco_gen_reldiff = ufloat('mupi_mass_reco_gen_reldiff'),
        lxy_reco_gen_reldiff = ufloat('lxy_reco_gen_reldiff'),
        ## gen displacement 
        gen_lxy = ufloat('gen_lxy'),
        diff_mass = ufloat('diff_mass'),
    )
)


CountHNLToMuPi = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("HNLToMuPi")
)    

HNLToMuPiSequence   = cms.Sequence( HNLToMuPi )
HNLToMuPiSequenceMC = cms.Sequence( HNLToMuPiMC )
