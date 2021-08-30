import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.genparticles_cff import finalGenParticles, genParticleTable
from PhysicsTools.BParkingNano.common_cff import Var


# for BPHPark start with merged particles (pruned + packed),
# where pruned contain K* states, but not final states, 
# and packed contain final states (K pi).
# then you save also final states (granddaughters)
finalGenParticlesBPark = finalGenParticles.clone(
  src = cms.InputTag("mergedGenParticles"),
  select = cms.vstring(
	"drop *",
   )
)

genParticleBParkTable = genParticleTable.clone(
  src = cms.InputTag("finalGenParticlesBPark"),
  variables = cms.PSet(
      genParticleTable.variables,
      p = Var("p",  float, precision=8),
      gamma = Var("p4().Gamma()",  float, precision=8),
      beta  = Var("p4().Beta()",  float, precision=8),
      vx = Var("vx()", float, doc="x coordinate of the production vertex position, in cm", precision=10),
      vy = Var("vy()", float, doc="y coordinate of the production vertex position, in cm", precision=10),
      vz = Var("vz()", float, doc="z coordinate of the production vertex position, in cm", precision=10),
  )
)

genParticleBParkSequence = cms.Sequence(finalGenParticlesBPark)
genParticleBParkTables = cms.Sequence(genParticleBParkTable)

