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
        'userInt("hnlVtxOK") == 1',
        'userFloat("hnlVtxProb") > 0.0001',
        'userFloat("hnlFittedCosTheta2D") >= 0.5',
        'userFloat("hnlFittedMass") > 0.5',
        'userFloat("hnlFittedMass") < 6.5',
        ])
    ),
)
    

BToMuMuPiTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToMuMuPi"),
    cut = cms.string(""),
    name = cms.string("Bmeson"),
    doc = cms.string("B candidate from trigger muon, and muon+pion from vertex"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities related to the B
        CandVars,  # pt, eta, phi, mass, charge, pdgId
        # 
        triggerMuonIdx = uint('triggerMuonIdx'), # naming convention of NanoAODs: <collection_name_starting_with_non_capital>Idx
        muonIdx = uint('muonIdx'),
        probeTracksIdx = uint('probeTracksIdx' ), # bad convention in BParkinNanoAOD, it should be ProbeTrack, not ProbeTracks
        # mu-pi fit and vtx info
        hnlVtxChi2    = ufloat('hnlVtxChi2' ),  # naming convention of NanoAODs: use "_" only for variable of a collection
        hnlVtxProb    = ufloat('hnlVtxProb' ),  #                                use "firstSecondThird" where each is a distinct entity
        hnlLxy     = ufloat('hnlLxy'     ),
        hnlLxyUnc    = ufloat('hnlLxyUnc'    ),
        hnlLsxy = ufloat('hnlLsxy'    ),
        hnlVtxX       = ufloat('hnlVtxX'    ),
        hnlVtxY       = ufloat('hnlVtxY'    ),
        hnlVtxZ       = ufloat('hnlVtxZ'    ),
        hnlVtxXE      = ufloat('hnlVtxXE'   ),  # only saving diagonal elements of the cov matrix
        hnlVtxYE      = ufloat('hnlVtxYE'   ),
        hnlVtxZE      = ufloat('hnlVtxZE'   ),
        # HNL: what is the difference between this and the fitted quantities ?
        hnlMass   = Var('daughter("hnl").mass()', float),
        hnlPt     = Var('daughter("hnl").pt()'  , float),
        hnlEta    = Var('daughter("hnl").eta()' , float),
        hnlPhi    = Var('daughter("hnl").phi()' , float),
        hnlCharge = Var('daughter("hnl").charge()', int),
        # Cos(theta)
        hnlCosTheta2D     = ufloat('hnlCosTheta2D'       ),
        hnlFittedCosTheta2D = ufloat('hnlFittedCosTheta2D'),
        # post-fit momentum
        hnlFittedMass  = ufloat('hnlFittedMass'   ), 
        hnlFittedMassErr = ufloat('hnlFittedMassErr'),
        hnlFittedPt    = ufloat('hnlFittedPt'     ),
        hnlFittedEta   = ufloat('hnlFittedEta'    ),
        hnlFittedPhi   = ufloat('hnlFittedPhi'    ),
        hnlFittedMuPt     = ufloat('hnlFittedMuPt'  ),
        hnlFittedMuEta    = ufloat('hnlFittedMuEta' ),
        hnlFittedMuPhi    = ufloat('hnlFittedMuPhi' ),
        hnlFittedPiPt     = ufloat('hnlFittedPiPt'  ),
        hnlFittedPiEta    = ufloat('hnlFittedPiEta' ),
        hnlFittedPiPhi    = ufloat('hnlFittedPiPhi' ),
    )
)

CountBToMuMuPi = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuMuPi")
)    

BToMuMuPiSequence = cms.Sequence( BToMuMuPi )
