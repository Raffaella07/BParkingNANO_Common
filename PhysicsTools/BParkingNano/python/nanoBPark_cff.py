from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from trgbits_cff import *



##for gen and trigger muon
from PhysicsTools.BParkingNano.genparticlesBPark_cff import finalGenParticlesBPark, genParticleBParkTable, genParticleBParkSequence, genParticleBParkTables
from PhysicsTools.BParkingNano.particlelevelBPark_cff import *
from PhysicsTools.BParkingNano.triggerObjectsBPark_cff import *
from PhysicsTools.BParkingNano.muonsBPark_cff import * 

## filtered input collections
from PhysicsTools.BParkingNano.electronsBPark_cff import * 
from PhysicsTools.BParkingNano.tracksBPark_cff import *

## B collections
from PhysicsTools.BParkingNano.BToMuLPi_cff import * # BToMuMuPi, BToMuMuPiMC, BToMuMuPiTable, BToMuMuPiSequence, BToMuMuPiSequenceMC, CountBToMuMuPi
from PhysicsTools.BParkingNano.BToMuLPiHD_cff import * # BToMuMuPi, BToMuMuPiMC, BToMuMuPiTable, BToMuMuPiSequence, BToMuMuPiSequenceMC, CountBToMuMuPi
from PhysicsTools.BParkingNano.BToKLL_cff import *
from PhysicsTools.BParkingNano.HNLToMuPi_cff import HNLToMuPi, HNLToMuPiMC, HNLToMuPiTable, HNLToMuPiSequence, HNLToMuPiSequenceMC, CountHNLToMuPi
from PhysicsTools.BParkingNano.BToKstarLL_cff import *
from PhysicsTools.BParkingNano.TagAndProbeJPsiToMuMu_cff import * 



nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)

nanoSequence = cms.Sequence(nanoMetadata + 
                            vertexSequence +           
                            globalTables + vertexTables + 
                            triggerObjectBParkTables + l1bits)

nanoSequenceMC = cms.Sequence(particleLevelBParkSequence + genParticleBParkSequence + 
                              globalTablesMC + genWeightsTable + genParticleBParkTables + lheInfoTable) 



def nanoAOD_customizeMuonTriggerBPark(process):
    #process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonBParkTables)#+ muonTriggerMatchedTables)   ###comment in this extra table in case you want to create the TriggerMuon collection again.
    process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonTriggerMatchedTables + muonBParkTables)
    return process

def nanoAOD_customizeTrackFilteredBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + tracksBParkSequence + tracksBParkTables)
    return process

def nanoAOD_customizeElectronFilteredBPark(process):
    process.nanoBKeeSequence     = cms.Sequence( electronsBParkSequence + electronBParkTables)

    process.nanoeSequence     = cms.Sequence( electronsBParkSequence + electronBParkTables) #new sequence defined just for the sake of clarity, to be used in the python config
    process.nanoBKstarEESequence = cms.Sequence( electronsBParkSequence + electronBParkTables)
    return process

def nanoAOD_customizeTriggerBitsBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + trgTables)
    return process
#including here sequences for reconstruction of HNL decaying in electron and muon, separately for MC and data 
def nanoAOD_customizeBToMuLPi(process, isMC=False):
   if isMC == False:
     process.nanoBMuMuPiSequence = cms.Sequence( BToMuMuPiSequence + BToMuMuPiTable + CountBToMuMuPi)
     process.nanoBMuEPiSequence  = cms.Sequence( BToMuEPiSequence  + BToMuEPiTable  + CountBToMuEPi)
#     process.nanoBMuMuPi_HDSequence = cms.Sequence( BToMuMuPi_HDSequence + BToMuMuPi_HDTable + CountBToMuMuPi_HD)
     process.nanoBMuEPiHDSequence  = cms.Sequence( BToMuEPiHDSequence  + BToMuEPiHDTable  + CountBToMuEPiHD)
   else:
     process.nanoBMuMuPiSequence = cms.Sequence( BToMuMuPiSequenceMC + BToMuMuPiMCTable + CountBToMuMuPiMC)
     process.nanoBMuEPiSequence  = cms.Sequence( BToMuEPiSequenceMC  + BToMuEPiMCTable  + CountBToMuEPiMC )
     process.nanoBMuEPiHDSequence  = cms.Sequence( BToMuEPiHDSequenceMC  + BToMuEPiHDMCTable  + CountBToMuEPiHDMC)
   return process

def nanoAOD_customizeHNLToMuPi(process, isMC=False):
    if isMC == False:
      process.nanoHNLToMuPiSequence = cms.Sequence( HNLToMuPiSequence + HNLToMuPiTable )
    else:
      process.nanoHNLToMuPiSequence = cms.Sequence( HNLToMuPiSequenceMC + HNLToMuPiTable )
    return process

def nanoAOD_customizeBToKMuMu(process, isMC=False):
    if isMC == False:
      process.nanoBKMuMuSequence = cms.Sequence( BToKMuMuSequence + BToKmumuTable )
#      process.nanoBMuMuPiSequence = cms.Sequence( BToMuMuPiSequence + BToMuMuPiTable )
    else:
      process.nanoBKMuMuSequence = cms.Sequence( BToKMuMuSequenceMC + BToKmumuTable )
    return process

def nanoAOD_customizeBToKLL(process):
    process.nanoBKeeSequence   = cms.Sequence( process.nanoBKeeSequence + BToKEESequence    + BToKeeTable   )
    process.nanoBKMuMuSequence = cms.Sequence( BToKMuMuSequence + BToKmumuTable )
    return process

def nanoAOD_customizeTagAndProbeJPsiToMuMu(process, isMC=False):
    if isMC == False:
      process.nanoJPsiToMuMuSequence = cms.Sequence( JPsiToMuMuSequence + JPsiToMuMuTable )
    else:
      process.nanoJPsiToMuMuSequence = cms.Sequence( JPsiToMuMuSequenceMC + JPsiToMuMuTable )
    return process

#three possibilities for K*LL
def nanoAOD_customizeBToKstarLL(process):
    process.nanoBKstarLLSequence   = cms.Sequence( KstarToKPiSequence + BToKstarLLSequence + KstarToKPiTable + BToKstarLLTables )
    return process

def nanoAOD_customizeBToKstarEE(process):
    process.nanoBKstarEESequence   = cms.Sequence( process.nanoBKstarEESequence + BToKstarEESequence + BToKstarEETable + KstarToKPiTable )
    return process

def nanoAOD_customizeBToKstarMuMu(process):
    process.nanoBKstarMuMuSequence = cms.Sequence( BToKstarMuMuSequence + BToKstarMuMuTable + KstarToKPiTable )
    return process

from FWCore.ParameterSet.MassReplace import massSearchReplaceAnyInputTag
#def nanoAOD_customizeMC(process, ancestor_particles=[511, 521, 531, 541, 211, 9900015]):  
def nanoAOD_customizeMC(process, ancestor_particles=[511, 521, 531, 541, 9900015]):  
    for name, path in process.paths.iteritems():
        # replace all the non-match embedded inputs with the matched ones
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:SelectedMuons', 'selectedMuonsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:trgMuons', 'triggerMuonsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'electronsForAnalysis:SelectedElectrons', 'selectedElectronsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'tracksBPark:SelectedTracks', 'tracksBParkMCMatchEmbedded')

        # make the BToMuMuPiTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'BToMuMuPi', 'BToMuMuPiMC')
        
        # make the BToKMuMuTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'BToKmumu', 'BToKmumuMC')

        # make the BToMuEPiTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'BToMuEPi', 'BToMuEPiMC') #adding electron channel
        # make the BToMuEPiHDTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'BToMuEPiHD', 'BToMuEPiHDMC') #adding electron channel
        # make the HNLToMuPiTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'HNLToMuPi', 'HNLToMuPiMC')

        # make the JPsiToMuMuTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'JPsiToMuMu', 'JPsiToMuMuMC')

        # save the all descendants of ancestor_particles
        to_save = ' || '.join(['abs(pdgId) == %d'%ipdg for ipdg in ancestor_particles])
        to_save = 'keep++ (%s)'%to_save
        finalGenParticlesBPark.select.append(to_save)

        # modify the path to include mc-specific info
        path.insert(0, nanoSequenceMC)
        path.replace(process.muonBParkSequence, process.muonBParkMC)
        path.replace(process.electronsBParkSequence, process.electronBParkMC)
        path.replace(process.tracksBParkSequence, process.tracksBParkMC)
