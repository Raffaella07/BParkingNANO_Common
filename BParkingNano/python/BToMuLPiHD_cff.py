import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import uint, ufloat, Var, CandVars

#comment out BToMuEPi and BToMuEPiMC and their tables when running on BToMuMuPiBuilder
BToMuEPiHD = cms.EDProducer(
     'BToMuTrkPiBuilder',
     trgMuons                = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
     leptons                 = cms.InputTag('tracksBPark', 'SelectedTracks'),
     leptonsTransientTracks  = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
     pions                   = cms.InputTag('tracksBPark', 'SelectedTracks'),
     pionsTransientTracks    = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
     tracks                  = cms.InputTag("packedPFCandidates"), 
     lostTracks              = cms.InputTag("lostTracks"), 
     genParticles            = cms.InputTag("finalGenParticlesBPark"),
     beamSpot                = cms.InputTag('offlineBeamSpot'), 
 
     label = cms.string('ele'),
     isMC  = cms.bool(False),
     # pre-fitter preselection
     #pionSelection           = cms.string('pt > 0.55 && abs(eta)<2'),  
     isoTracksSelection      = cms.string('pt > 0.55 && abs(eta)<2'),
     trgMuonSelection        = cms.string('pt > 7 && abs(eta) < 1.5 '),
     #trgMuonSelection        = cms.string('pt > 7 && abs(eta) < 1.5 && softId==1'),
    # leptonSelection        = cms.string('pt > 1.5 && abs(eta) < 2'),# not fully exploiting lowpt full spectrum, to check

    pionSelection        = cms.string(' & '.join([
#	'pt>0.7',
	'abs(eta)<2',
        'abs(userFloat("dz")) > 0.008',
        'abs(userFloat("dxy")) > 0.008',
        'abs(userFloat("dzS")) > 1',
        'abs(userFloat("dxyS")) > 3',
	])
	),		
    leptonSelection        = cms.string(' & '.join([
       'abs(eta)<2',
       'abs(userFloat("dz")) > 0.008',
       'abs(userFloat("dxy")) > 0.008',
       'abs(userFloat("dzS")) > 1',
       'abs(userFloat("dxyS")) > 3',
	])
	),		
     preVtxSelection = cms.string(' & '.join([
	'max(userFloat("pi_dxy"),userFloat("trk_dxy"))>0.01',
        #'charge ==0',
         'mass > 0.5',        
         'mass < 7.0',        
         ])
     ), # applied on the HNL cand 
 
     # post-fitter preselection
     postVtxSelection = cms.string(' & '.join([
         'userInt("hnl_vtx_OK") == 1',
         'userFloat("hnl_vtx_prob") > 0.01',
        # 'userFloat("hnl_l_xy") > 10 ',
          'mass < 8',
        # 'pt > 11',
        # 'abs(eta) < 1.7',
          'userFloat("hnl_vtx_chi2") < 9', #<----------------ulterior selections have been commented out for a clean slate production 
      #   'userFloat("trg_muon_sip3d") > 0.8',
      #  'userFloat("sel_lep_ip3d") > 0.0015',
      #  'abs(userFloat("sel_lep_dxy")) > 0.0005',
      #  'abs(userFloat("pion_dz")) > 0.001',
      #  'abs(userFloat("pion_dxy")) > 0.0003',
      #  'abs(userFloat("pion_dzS")) > 0.5',
      #  'abs(userFloat("pion_dxyS")) > 0.2',
      #  'abs(userFloat("pion_DCASig")) > 0.3',
         ##'userFloat("muons_Lxyz") > 0.005',
         ##'userFloat("pion_muon_vzdiff") > 0.001',
         ##'userFloat("dr_mu_pi") < 1.8',
       # 'userFloat("dr_trgmu_hnl") < 0.8',
        'userFloat("hnl_fitted_mass") > 0.5',
        'userFloat("hnl_fitted_mass") < 7',
       # 'userFloat("hnl_fitted_pt") > 4',
       # 'abs(userFloat("hnl_fitted_eta")) < 1.8',
         'userFloat("hnl_ls_xy") > 20',
        'userFloat("hnl_fitted_cos_theta_2D") > 0.99',
         ])
     ), # applied on the B cand
)
     
BToMuEPiHDMC = BToMuEPiHD.clone(
     isMC = cms.bool(True),
)

#BToMuMuPiHDTable = cms.EDProducer(
#    'SimpleCompositeCandidateFlatTableProducer',
#    src = cms.InputTag("BToMuMuPiHD"),
#    cut = cms.string(""),
#    name = cms.string("BToMuMuPiHD"),
#    doc = cms.string("HNL Variable"),
#    singleton=cms.bool(False),
#    extension=cms.bool(False),
#    variables=cms.PSet(         #here HNL mu variables are labeled r "sel_mu*"
#        # pre-fit quantities
#        CandVars,
#        trg_mu_idx      = uint('trg_mu_idx'),
#        sel_mu_idx      = uint('trk_idx'), 
#        pi_idx          = uint('pi_idx'    ), 
#        ## trigger muon
#   #    trg_mu_pt       = ufloat('trg_muon_pt'   ), 
#   #    trg_mu_eta      = ufloat('trg_muon_eta'  ), 
#   #    trg_mu_phi      = ufloat('trg_muon_phi'  ), 
#        ## vertex difference between the two muons
#        dimu_vxdiff     = ufloat('dilepton_vxdiff' ),
#        dimu_vydiff     = ufloat('dilepton_vydiff' ),
#:wq
     #   dimu_vzdiff     = ufloat('dilepton_vzdiff' ),
#        dimu_Lxy        = ufloat('dilepton_Lxy'    ),
#        dimu_Lxyz       = ufloat('dilepton_Lxyz'   ),
#        ## vertex difference between the trigger muon and pion
#        pi_mu_vzdiff    = ufloat('pion_trgmuon_vzdiff'  ),
#        # post-fit quantities
#        ## vertex information 
#        sv_chi2         = ufloat('hnl_vtx_chi2' ),
#        sv_prob         = ufloat('hnl_vtx_prob' ),
#        sv_lxy          = ufloat('hnl_l_xy'     ),
#        sv_lxye         = ufloat('hnl_l_xy_unc' ),
#        sv_lxy_sig      = ufloat('hnl_ls_xy'    ),
#        sv_x            = ufloat('hnl_vtx_x'    ),
#        sv_y            = ufloat('hnl_vtx_y'    ),
#        sv_z            = ufloat('hnl_vtx_z'    ),
#        sv_xe           = ufloat('hnl_vtx_ex'   ), ## only saving diagonal elements of the cov matrix
#        sv_ye           = ufloat('hnl_vtx_ey'   ),
#        sv_ze           = ufloat('hnl_vtx_ez'   ),
#        ## HNL (only postfit information is saved) 
#        hnl_mass        = ufloat('hnl_fitted_mass'   ),
#        hnl_masserr     = ufloat('hnl_fitted_massErr'),
#        hnl_pt          = ufloat('hnl_fitted_pt'     ),
#        hnl_eta         = ufloat('hnl_fitted_eta'    ),
#        hnl_phi         = ufloat('hnl_fitted_phi'    ),
#        hnl_charge      = Var('daughter("hnl").charge()', int),
#        hnl_cos2D       = ufloat('hnl_fitted_cos_theta_2D'   ),
#        ## daughter muon
#        fit_mu_pt       = ufloat('hnl_fitted_lep_pt'  ), 
#        fit_mu_eta      = ufloat('hnl_fitted_lep_eta' ),
#        fit_mu_phi      = ufloat('hnl_fitted_lep_phi' ),
#        fit_mu_mass     = ufloat('hnl_fitted_lep_mass'),
#        ## daughter pion
#        fit_pi_pt       = ufloat('hnl_fitted_pi_pt'  ),
#        fit_pi_eta      = ufloat('hnl_fitted_pi_eta' ),
#        fit_pi_phi      = ufloat('hnl_fitted_pi_phi' ),
#        fit_pi_mass     = ufloat('hnl_fitted_pi_mass'),
#        ## dR quantities
#        dr_mu_pi        = ufloat('dr_lep_pi'     ),
#        dr_trgmu_hnl    = ufloat('dr_trgmu_hnl'  ),
#        # Other quantities
#        ## ID WP of the selected muon 
#   #    sel_mu_isSoft   =  ufloat('sel_muon_isSoft'   ),
#   # 	sel_mu_isTight  = ufloat('sel_muon_isTight'  ),
#   #    sel_mu_isMedium = ufloat('sel_muon_isMedium' ),
#   #    sel_mu_isLoose  = ufloat('sel_muon_isLoose'  ),
#   #    ## impact paramaters
#    #   trg_mu_ip3d     = ufloat('trg_muon_ip3d'  ), 
#    #   trg_mu_sip3d    = ufloat('trg_muon_sip3d' ), 
#    #   trg_mu_dxy      = ufloat('trg_muon_dxy'   ), 
#    #   trg_mu_dz       = ufloat('trg_muon_dz'    ), 
#    #   sel_lep_ip3d     = ufloat('sel_lep_ip3d'  ), 
#    #   sel_lep_sip3d    = ufloat('sel_lep_sip3d' ), 
#    #   sel_lep_dxy      = ufloat('sel_lep_dxy'   ), 
#    #   sel_lep_dz       = ufloat('sel_lep_dz'    ), 
#    #   pi_dz           = ufloat('pion_dz'        ), 
#    #   pi_dxy          = ufloat('pion_dxy'       ), 
#    #   pi_dzS          = ufloat('pion_dzS'       ), 
#    #   pi_dxyS         = ufloat('pion_dxyS'      ), 
#    #   pi_DCASig       = ufloat('pion_DCASig'    ), 
#        ## isolation
#        trg_mu_iso03    = ufloat('trg_mu_iso03'   ), 
#        trg_mu_iso04    = ufloat('trg_mu_iso04'   ),
#        sel_mu_iso03    = ufloat('trk_iso03'      ),
#        sel_mu_iso04    = ufloat('trk_iso04'      ),
#        pi_iso03        = ufloat('pi_iso03'       ),
#        pi_iso04        = ufloat('pi_iso04'       ),
#        hnl_iso03       = ufloat('hnl_iso03'      ),
#        hnl_iso04       = ufloat('hnl_iso04'      ),
#        trg_mu_iso03_close  = ufloat('trg_mu_iso03_close' ), 
#        trg_mu_iso04_close  = ufloat('trg_mu_iso04_close' ),
#        sel_mu_iso03_close  = ufloat('sel_trk_iso03_close' ),
#        sel_mu_iso04_close  = ufloat('sel_trk_iso04_close' ),
#        pi_iso03_close      = ufloat('pi_iso03_close'     ),
#        pi_iso04_close      = ufloat('pi_iso04_close'     ),
#        hnl_iso03_close     = ufloat('hnl_iso03_close'    ),
#        hnl_iso04_close     = ufloat('hnl_iso04_close'    ),
#        ## dilepton mass
#        dilepton_mass   = ufloat('dilepton_mass'  ),
#        dilepton_pt     = ufloat('dilepton_pt'    ),
#        ## gen-matching
# #      isMatched                   = Var("userInt('isMatched')"                  , int, mcOnly=True),
# #      matching_sel_mu_genIdx      = Var("userInt('matching_sel_lep_genIdx')"     , int, mcOnly=True),
# #      matching_trg_mu_genIdx      = Var("userInt('matching_trg_mu_genIdx')"     , int, mcOnly=True),
# #      matching_pi_genIdx          = Var("userInt('matching_pi_genIdx')"         , int, mcOnly=True),
# #      matching_sel_mu_motherPdgId = Var("userInt('matching_sel_lep_motherPdgId')", int, mcOnly=True),
# #      matching_trg_mu_motherPdgId = Var("userInt('matching_trg_mu_motherPdgId')", int, mcOnly=True),
# #      matching_pi_motherPdgId     = Var("userInt('matching_pi_motherPdgId')"    , int, mcOnly=True),
# #      ## displacement 
#        #fitter_bs_lxy = ufloat('fitter_bs_lxy'), #same as sv_lxy, with more explicit naming
#        ##my_fitter_bs_lxy = ufloat('my_fitter_bs_lxy'), 
#        ##disp2DFromBS = ufloat('disp2DFromBS'),
#        ##flightdisp = ufloat('flightdisp'),
#        #gendisp_trgmu_mu_lxyz = ufloat('gendisp_trgmu_mu_lxyz'),
#        #gendisp_trgmu_mu_lxy = ufloat('gendisp_trgmu_mu_lxy'),
#    )
#)
#
#BToMuMuPiHDMCTable = BToMuMuPiHDTable.clone(
#	src = cms.InputTag("BToMuMuPiHDMC"),
#        isMC = cms.bool(True),	
#)
#
BToMuEPiHDTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToMuEPiHD"),
    cut = cms.string(""),
    name = cms.string("BToMuEPiHD"),
    doc = cms.string("HNL Variable high displ"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(           #here HNL mu variables are labeled r "sel_e*"
        # pre-fit quantities
        CandVars,
        trg_mu_idx      = uint('trg_mu_idx'),
        sel_e_idx      = uint('lep_idx'), 
        pi_idx          = uint('pi_idx'    ), 
        ## trigger muon
   #    trg_mu_pt       = ufloat('trg_muon_pt'   ), 
   #    trg_mu_eta      = ufloat('trg_muon_eta'  ), 
   #    trg_mu_phi      = ufloat('trg_muon_phi'  ), 
        ## vertex difference between the two muons
        dilep_vxdiff     = ufloat('dilepton_vxdiff' ),
        dilep_vydiff     = ufloat('dilepton_vydiff' ),
        dilep_vzdiff     = ufloat('dilepton_vzdiff' ),
        dilep_Lxy        = ufloat('dilepton_Lxy'    ),
        dilep_Lxyz       = ufloat('dilepton_Lxyz'   ),
        ## vertex difference between the trigger muon and pion
        pi_mu_vzdiff    = ufloat('pion_trgmuon_vzdiff'  ),
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
        fit_l_pt       = ufloat('hnl_fitted_lep_pt'  ), 
        fit_l_eta      = ufloat('hnl_fitted_lep_eta' ),
        fit_l_phi      = ufloat('hnl_fitted_lep_phi' ),
        fit_l_mass     = ufloat('hnl_fitted_lep_mass'),
        ## daughter pion
        fit_pi_pt       = ufloat('hnl_fitted_pi_pt'  ),
        fit_pi_eta      = ufloat('hnl_fitted_pi_eta' ),
        fit_pi_phi      = ufloat('hnl_fitted_pi_phi' ),
        fit_pi_mass     = ufloat('hnl_fitted_pi_mass'),
        ## dR quantities
        dr_mu_pi        = ufloat('dr_lep_pi'     ),
        dr_trgmu_hnl    = ufloat('dr_trgmu_hnl'  ),
        # Other quantities
        ## ID WP of the selected muon 
   #    sel_mu_isSoft   =  ufloat('sel_muon_isSoft'   ),
   # 	sel_mu_isTight  = ufloat('sel_muon_isTight'  ),
   #    sel_mu_isMedium = ufloat('sel_muon_isMedium' ),
   #    sel_mu_isLoose  = ufloat('sel_muon_isLoose'  ),
   #    ## impact paramaters
    #   trg_mu_ip3d     = ufloat('trg_muon_ip3d'  ), 
    #   trg_mu_sip3d    = ufloat('trg_muon_sip3d' ), 
    #   trg_mu_dxy      = ufloat('trg_muon_dxy'   ), 
    #   trg_mu_dz       = ufloat('trg_muon_dz'    ), 
    #   sel_lep_ip3d     = ufloat('sel_lep_ip3d'  ), 
    #   sel_lep_sip3d    = ufloat('sel_lep_sip3d' ), 
    #   sel_lep_dxy      = ufloat('sel_lep_dxy'   ), 
    #   sel_lep_dz       = ufloat('sel_lep_dz'    ), 
    #   pi_dz           = ufloat('pion_dz'        ), 
    #   pi_dxy          = ufloat('pion_dxy'       ), 
    #   pi_dzS          = ufloat('pion_dzS'       ), 
    #   pi_dxyS         = ufloat('pion_dxyS'      ), 
    #   pi_DCASig       = ufloat('pion_DCASig'    ), 
        ## isolation
        trg_mu_iso03    = ufloat('trg_mu_iso03'   ), 
        trg_mu_iso04    = ufloat('trg_mu_iso04'   ),
        sel_e_iso03    = ufloat('trk_iso03'      ),
        sel_e_iso04    = ufloat('trk_iso04'      ),
        pi_iso03        = ufloat('pi_iso03'       ),
        pi_iso04        = ufloat('pi_iso04'       ),
        hnl_iso03       = ufloat('hnl_iso03'      ),
        hnl_iso04       = ufloat('hnl_iso04'      ),
        trg_mu_iso03_close  = ufloat('trg_mu_iso03_close' ), 
        trg_mu_iso04_close  = ufloat('trg_mu_iso04_close' ),
        sel_e_iso03_close  = ufloat('sel_trk_iso03_close' ),
        sel_e_iso04_close  = ufloat('sel_trk_iso04_close' ),
        pi_iso03_close      = ufloat('pi_iso03_close'     ),
        pi_iso04_close      = ufloat('pi_iso04_close'     ),
        hnl_iso03_close     = ufloat('hnl_iso03_close'    ),
        hnl_iso04_close     = ufloat('hnl_iso04_close'    ),
        ## dilepton mass
        dilepton_mass   = ufloat('dilepton_mass'  ),
        dilepton_pt     = ufloat('dilepton_pt'    ),
        ## gen-matching
        isMatched                   = Var("userInt('isMatched')"                  , int, mcOnly=True),
  #     matching_sel_e_genIdx      = Var("userInt('matching_sel_lep_genIdx')"     , int, mcOnly=True),
  #     matching_trg_mu_genIdx      = Var("userInt('matching_trg_mu_genIdx')"     , int, mcOnly=True),
  #     matching_pi_genIdx          = Var("userInt('matching_pi_genIdx')"         , int, mcOnly=True),
  #     matching_sel_e_motherPdgId = Var("userInt('matching_sel_lep_motherPdgId')", int, mcOnly=True),
  #     matching_trg_mu_motherPdgId = Var("userInt('matching_trg_mu_motherPdgId')", int, mcOnly=True),
  #     matching_pi_motherPdgId     = Var("userInt('matching_pi_motherPdgId')"    , int, mcOnly=True),
        ## displacement 
        #fitter_bs_lxy = ufloat('fitter_bs_lxy'), #same as sv_lxy, with more explicit naming
        ##my_fitter_bs_lxy = ufloat('my_fitter_bs_lxy'), 
        ##disp2DFromBS = ufloat('disp2DFromBS'),
        ##flightdisp = ufloat('flightdisp'),
        #gendisp_trgmu_mu_lxyz = ufloat('gendisp_trgmu_mu_lxyz'),
        #gendisp_trgmu_mu_lxy = ufloat('gendisp_trgmu_mu_lxy'),
    )
)


BToMuEPiHDMCTable = BToMuEPiHDTable.clone(
	src = cms.InputTag("BToMuEPiHDMC"),
        isMC = cms.bool(True),	
)


#CountBToMuMuPiHD = cms.EDFilter("PATCandViewCountFilter",
#    minNumber = cms.uint32(1),
#    maxNumber = cms.uint32(999999),
#    src = cms.InputTag("BToMuMuPiHD" )
#)    

CountBToMuEPiHD = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuEPiHD")
)   
 
#ountBToMuMuPiHDMC = cms.EDFilter("PATCandViewCountFilter",
#   minNumber = cms.uint32(1),
#   maxNumber = cms.uint32(999999),
#   src = cms.InputTag("BToMuMuPiHDMC" )
#    

CountBToMuEPiHDMC = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuEPiHDMC")
)   
#BToMuMuPiHDSequence   = cms.Sequence( BToMuMuPiHD ) #maybe double sequence for MC and data can be avoided? Works fine so for now keeping it that way
#BToMuMuPiHDSequenceMC = cms.Sequence( BToMuMuPiHDMC )
BToMuEPiHDSequence   = cms.Sequence( BToMuEPiHD  )
BToMuEPiHDSequenceMC = cms.Sequence( BToMuEPiHDMC )
