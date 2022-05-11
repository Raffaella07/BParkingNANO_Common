import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *


JPsiToMuMu = cms.EDProducer(
    'TagAndProbeJPsiToMuMuBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    genParticles = cms.InputTag("finalGenParticlesBPark"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    isMC = cms.bool(False),
    lep1Selection = cms.string(''), #cms.string('pt > 1.5'),
    lep2Selection = cms.string(''),
    #preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
    #                             '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'),
    #preVtxSelection = cms.string('userFloat("lep_deltaR") > 0.03'),
    preVtxSelection = cms.string(''),
    #postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5 && mass>3.15 && mass<2.95 && userFloat("cos_theta_2D")>0.9'),
    #postVtxSelection = cms.string('userFloat("sv_prob") > 1.e-5 && userFloat("fitted_mass")>2.9 && userFloat("fitted_mass")<3.3 && userFloat("cos_theta_2D")>0.9'),
    postVtxSelection = cms.string('userFloat("sv_prob") > 1.e-5 && userFloat("fitted_mass")>2.9 && userFloat("fitted_mass")<3.3 && userFloat("cos_theta_2D")>0.9'),
    #postVtxSelection = cms.string(''),
)


JPsiToMuMuMC = JPsiToMuMu.clone( 
    isMC = cms.bool(True),
)


JPsiToMuMuTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("JPsiToMuMu"),
    cut = cms.string(""),
    name = cms.string("JPsiToMuMu"),
    doc = cms.string("JPsiToMuMu Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        sv_chi2 = ufloat('sv_chi2'),
        sv_ndof = ufloat('sv_ndof'),
        sv_prob = ufloat('sv_prob'),
        mass = ufloat('fitted_mass'),
        pt = ufloat('fitted_pt'),
        eta = ufloat('fitted_eta'),
        phi = ufloat('fitted_phi'),
        charge = ufloat('charge'),
        cos2D = ufloat('cos_theta_2D'),
        lep1_idx = uint('l1_idx'),
        lep2_idx = uint('l2_idx'),
        lep1_pt = ufloat('fitted_lep1_pt'),
        lep1_eta = ufloat('fitted_lep1_eta'),
        lep1_phi = ufloat('fitted_lep1_phi'),
        lep2_pt = ufloat('fitted_lep2_pt'),
        lep2_eta = ufloat('fitted_lep2_eta'),
        lep2_phi = ufloat('fitted_lep2_phi'),
        lxy = ufloat('l_xy'),
        lxy_sig = ufloat('l_xy_sig'),
        deltaR = ufloat('lep_deltaR'),
        lep_vzdiff = ufloat('lep_vzdiff'),
        isMatched = uint('isMatched'),
    )
)


CountJPsiToMuMu = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("JPsiToMuMu")
) 


JPsiToMuMuSequence = cms.Sequence(JPsiToMuMu)
JPsiToMuMuSequenceMC = cms.Sequence(JPsiToMuMuMC)
JPsiToMuMuTables = cms.Sequence(JPsiToMuMuTable)

