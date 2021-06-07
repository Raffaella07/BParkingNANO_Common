import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *


Path=["HLT_Mu7_IP4","HLT_Mu8_IP6","HLT_Mu8_IP5","HLT_Mu8_IP3","HLT_Mu8p5_IP3p5","HLT_Mu9_IP6","HLT_Mu9_IP5","HLT_Mu9_IP4","HLT_Mu10p5_IP3p5","HLT_Mu12_IP6"]

muonTrgSelector = cms.EDProducer("MuonTriggerSelector",
                                 muonCollection = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
                                 displacedStandaloneMuonCollection = cms.InputTag("displacedStandAloneMuons"), #same collection as in NanoAOD                                                           

                                 # trigger muon matching conditions
                                 max_deltaR = cms.double(0.15),
                                 max_deltaPtRel = cms.double(0.15),
                                 
                                 ## for the output selected collection (tag + all compatible in dZ)
                                 # difference of the vz of the trigger muon with selected muon
                                 dzForCleaning_wrtTrgMuon = cms.double(-1), # initial value: 1.8
                                 
                                 # selection for the selected muon
                                 selmu_ptMin = cms.double(1.),
                                 #selmu_ptMin = cms.double(0.3),
                                 selmu_absEtaMax = cms.double(2.3),
                                 #selmu_absEtaMax = cms.double(2.8),
                                 HLTPaths=cms.vstring(Path),
                                 #L1seeds=cms.vstring(Seed),
                             )
#cuts minimun number in B both mu and e, min number of trg, dz muon, dz and dr track, 
countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("muonTrgSelector", "trgMuons")
)


muonBParkTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("muonTrgSelector:SelectedMuons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Muon"),
    doc  = cms.string("slimmedMuons for BPark after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(CandVars,
        vx = Var("vx()", float, doc="x coordinate of vertex position, in cm", precision=6),
        vy = Var("vy()", float, doc="y coordinate of vertex position, in cm", precision=6),
        vz = Var("vz()", float, doc="z coordinate of vertex position, in cm", precision=6),
        ptErr = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
        dz = Var("userFloat('dz')", float, doc="dz (with sign) wrt first PV, in cm", precision=10),
        #dzErr = Var("userFloat('dzErr')", float, doc="dz uncertainty, in cm", precision=6),
        dzS = Var("userFloat('dzS')", float, doc="dz (with sign) significance wrt first PV, in cm", precision=10),
        dxy = Var("userFloat('dxy')", float, doc="dxy (with sign) wrt first PV, in cm", precision=10),
        #dxyErr = Var("userFloat('dxyErr')", float, doc="dxy uncertainty, in cm", precision=6),
        dxyS = Var("userFloat('dxyS')", float, doc="dxy (with sign) significance wrt first PV, in cm", precision=10),
        ip3d = Var("userFloat('ip3d')", float, doc="3D impact parameter wrt first PV, in cm", precision=10),
        sip3d = Var("userFloat('sip3d')", float, doc="3D impact parameter significance wrt first PV", precision=10),
        #pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component"),
        pfiso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))", float, doc="PF isolation dR=0.3, total (deltaBeta corrections)"),
        pfiso03Rel_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt", float, doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
        pfiso03_trk = Var("isolationR03().sumPt", float, doc = "raw isolationR03"),
        pfiso03Rel_trk = Var("isolationR03().sumPt/pt", float, doc = "raw isolationR03 relative"),
        pfiso03_ch = Var("pfIsolationR03().sumChargedHadronPt", float, doc = "raw pf iso charged hadrons"),
        pfiso03Rel_ch = Var("pfIsolationR03().sumChargedHadronPt/pt", float, doc = "raw pf iso charged hadrons"),
        pfiso03_n = Var("pfIsolationR03().sumNeutralHadronEt", float, doc = "raw iso neutral hadron"),
        pfiso03Rel_n = Var("pfIsolationR03().sumNeutralHadronEt/pt", float, doc = "raw iso neutral hadron"),
        pfiso03_pho = Var("pfIsolationR03().sumPhotonEt", float, doc = "raw pf iso photons"),
        pfiso03Rel_pho = Var("pfIsolationR03().sumPhotonEt/pt", float, doc = "raw pf iso photons"),
        pfiso03_pu = Var("pfIsolationR03().sumPUPt", float, doc = "raw pf iso PU"),
        pfiso03Rel_pu = Var("pfIsolationR03().sumPUPt/pt", float, doc = "raw pf iso PU"),
        pfiso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))", float, doc="PF isolation dR=0.4, total (deltaBeta corrections)"),
        pfiso04Rel_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt", float, doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
        #pfiso05_trk = Var("isolationR05().sumPt",float, doc = "raw isolationR03"),
        #pfiso05Rel_trk = Var("isolationR05().sumPt/pt",float, doc = "raw isolationR03 relative"),
        pfiso05_trk = Var("isolationR05().sumPt", float, doc = "raw isolationR05"),
        pfiso05Rel_trk = Var("isolationR05().sumPt/pt", float, doc = "raw relative isolationR05"),
        pfiso04_ch = Var("pfIsolationR04().sumChargedHadronPt", float, doc = "raw pf iso charged hadrons"),
        pfiso04Rel_ch = Var("pfIsolationR04().sumChargedHadronPt/pt", float, doc = "raw relative pf iso charged hadrons"),
        pfiso04_n= Var("pfIsolationR04().sumNeutralHadronEt", float, doc = "raw pf iso neutral hadron"),
        pfiso04Rel_n = Var("pfIsolationR04().sumNeutralHadronEt/pt", float, doc = "raw pf iso neutral hadron"),
        pfiso04_pho = Var("pfIsolationR04().sumPhotonEt", float, doc = "raw pf iso photons"),
        pfiso04Rel_pho = Var("pfIsolationR04().sumPhotonEt/pt", float, doc = "raw pf iso photons"),
        pfiso04_pu = Var("pfIsolationR04().sumPUPt", float, doc = "raw pf iso PU"),
        pfiso04Rel_pu = Var("pfIsolationR04().sumPUPt/pt", float, doc = "raw pf iso PU"),

        #tightCharge = Var("?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0",int,doc="Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)"),
        isPF = Var("isPFMuon", bool, doc="muon is PF candidate"),
        isGlobalMuon = Var("isGlobalMuon", bool, doc="muon is global muon"),
        isTrackerMuon = Var("isTrackerMuon", bool, doc="muon is tracker muon"),
        isGlobalOrTrackerMuon = Var("userInt('isGlobalOrTrackerMuon')", bool, doc="muon is global muon or tracker muon"),
        isGlobalNotTrackerMuon = Var("userInt('isGlobalNotTrackerMuon')", bool, doc="muon is global muon and not tracker muon"),
        isTrackerNotGlobalMuon = Var("userInt('isTrackerNotGlobalMuon')", bool, doc="muon is tracker muon and not global muon"),
        looseId = Var("userInt('looseId')", int, doc="reco muon is Loose"),
        mediumId = Var("passed('CutBasedIdMedium')", bool, doc="cut-based ID, medium WP"),
        #mediumPromptId = Var("passed('CutBasedIdMediumPrompt')",bool,doc="cut-based ID, medium prompt WP"),
        tightId = Var("passed('CutBasedIdTight')", bool, doc="cut-based ID, tight WP"),
        softId = Var("passed('SoftCutBasedId')", bool, doc="soft cut-based ID"),
        pfIsoId = Var("passed('PFIsoVeryLoose')+passed('PFIsoLoose')+passed('PFIsoMedium')+passed('PFIsoTight')+passed('PFIsoVeryTight')+passed('PFIsoVeryVeryTight')", "uint8", doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"),
        tkIsoId = Var("?passed('TkIsoTight')?2:passed('TkIsoLoose')", "uint8", doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"),  
        triggerIdLoose = Var("passed('TriggerIdLoose')", bool, doc="TriggerIdLoose ID"),
        #highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"),
        #softMvaId = Var("passed('SoftMvaId')",bool,doc="soft MVA ID"),
        #mvaId = Var("passed('MvaLoose')+passed('MvaMedium')+passed('MvaTight')","uint8",doc="Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight)"),
        #mvaId = Var("passed('MvaLoose')+passed('MvaMedium')+passed('MvaTight')","uint8",doc="Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight)"),
        #miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"),
        #multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"),

        inTimeMuon = Var("passed('InTimeMuon')", bool, doc="inTimeMuon ID"),
        #segmentCompatibility = Var("segmentCompatibility()", float, doc = "muon segment compatibility: propagating the tracker tracks to the muon system and evaluate the number of matched segments and the closeness of the matching", precision=14), # keep higher precision since people have cuts with 3 digits on this
        segmentCompatibility = Var("userFloat('segmentCompatibility')", float, doc = "muon segment compatibility: propagating the tracker tracks to the muon system and evaluate the number of matched segments and the closeness of the matching", precision=14), # keep higher precision since people have cuts with 3 digits on this
        caloCompatibility = Var("userFloat('caloCompatibility')", float, doc = "calorimetric compatibility"),
        validHitFraction = Var("userFloat('validHitFraction')", float, doc = "fraction of hits a tracker track uses (among inner tracker layers it traverses)"),
        kinkFinderChi2 = Var("userFloat('kinkFinderChi2')", float, doc = "chi2 of kink-finding algorithm: how likely it is that a track is made of more than one single track"),
        globalNormalisedChi2 = Var("userFloat('globalNormalisedChi2')", float, doc = "chi2/ndof of global fit"),
        localPositionChi2 = Var("userFloat('localPositionChi2')", float, doc = "chi2 of the position match between the tracker muon and the standalone muon"),
        trackerHighPurityFlag = Var("userInt('trackerHighPurityFlag')", int, doc = "tracker high-purity flag"), # int? 
        numberOfValidMuonHits = Var("userInt('numberOfValidMuonHits')", int, doc = "number of hits in the muon stations"),
        numberOfValidPixelHits = Var("userInt('numberOfValidPixelHits')", int, doc = "number of pixel hits"),
        numberOfTrackerLayers = Var("userInt('numberOfTrackerLayers')", int, doc = "number of tracker layers with hits"),
        numberOfPixelLayers = Var("userInt('numberOfPixelLayers')", int, doc = "number of pixel layers with hits"),
        numberOfStations = Var("userInt('nStations')", int, doc = "number of matched stations with default arbitration (segment & track)"),

        #nStations = Var("numberOfMatchedStations", int, doc = "number of matched stations with default arbitration (segment & track)"),
        #nTrackerLayers = Var("innerTrack().hitPattern().trackerLayersWithMeasurement()", int, doc = "number of layers in the tracker"),
        isTriggering = Var("userInt('isTriggering')", int, doc="flag the reco muon is also triggering"),
        matched_dr = Var("userFloat('DR')", float, doc="dr with the matched triggering muon" ),
        matched_dpt = Var("userFloat('DPT')", float, doc="dpt/pt with the matched triggering muon" ),   
        fired_HLT_Mu7_IP4 = Var("userInt('HLT_Mu7_IP4')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu8_IP6 = Var("userInt('HLT_Mu8_IP6')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu8_IP5 = Var("userInt('HLT_Mu8_IP5')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu8_IP3 = Var("userInt('HLT_Mu8_IP3')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu8p5_IP3p5 = Var("userInt('HLT_Mu8p5_IP3p5')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu9_IP6 = Var("userInt('HLT_Mu9_IP6')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu9_IP5 = Var("userInt('HLT_Mu9_IP5')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu9_IP4 = Var("userInt('HLT_Mu9_IP4')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu10p5_IP3p5 = Var("userInt('HLT_Mu10p5_IP3p5')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu12_IP6 = Var("userInt('HLT_Mu12_IP6')", int, doc="reco muon fired this trigger"),
    ),
)


muonsBParkMCMatchForTable = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = muonBParkTable.src,                         # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),                             # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.03),                           # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
    minPt = muonTrgSelector.selmu_ptMin,
    maxEta = muonTrgSelector.selmu_absEtaMax,
)

muonBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = muonBParkTable.src,
    mcMap   = cms.InputTag("muonsBParkMCMatchForTable"),
    objName = muonBParkTable.name,
    objType = muonBParkTable.name, 
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

selectedMuonsMCMatchEmbedded = cms.EDProducer(
    'MuonMatchEmbedder',
    src = cms.InputTag('muonTrgSelector:SelectedMuons'),
    matching = cms.InputTag('muonsBParkMCMatchForTable')
)

muonTriggerBParkTable = muonBParkTable.clone(
    src = cms.InputTag("muonTrgSelector:trgMuons"),
    name = cms.string("TriggerMuon"),
    doc  = cms.string("HLT Muons matched with reco muons"), #reco muon matched to triggering muon"),
    variables = cms.PSet(CandVars,
        vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6)####################,
#       trgMuonIndex = Var("userInt('trgMuonIndex')", int,doc="index in trigger muon collection")
   )
)

muonsTriggerBParkMCMatchForTable = cms.EDProducer("MCMatcher",# cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = muonTriggerBParkTable.src,                  # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),                             # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.03),                           # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
    minPt = cms.double(0.),
    maxEta = cms.double(10.),
)

muonTriggerBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = muonTriggerBParkTable.src,
    mcMap   = cms.InputTag("muonsTriggerBParkMCMatchForTable"),
    objName = muonTriggerBParkTable.name,
    objType = muonTriggerBParkTable.name, 
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

triggerMuonsMCMatchEmbedded = cms.EDProducer(
    'TriggerMuonMatchEmbedder',
    src = cms.InputTag('muonTrgSelector', 'trgMuons'),
    matching = cms.InputTag('muonsTriggerBParkMCMatchForTable')
)

muonBParkSequence = cms.Sequence(muonTrgSelector * countTrgMuons)
muonBParkMC = cms.Sequence(muonTrgSelector + muonsBParkMCMatchForTable + selectedMuonsMCMatchEmbedded + muonBParkMCTable + muonsTriggerBParkMCMatchForTable + triggerMuonsMCMatchEmbedded + muonTriggerBParkMCTable* countTrgMuons)
muonBParkTables = cms.Sequence(muonBParkTable)
muonTriggerMatchedTables = cms.Sequence(muonTriggerBParkTable)
