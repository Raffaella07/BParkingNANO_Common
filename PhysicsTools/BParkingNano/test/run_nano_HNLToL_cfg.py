# E.g.: cmsRun run_nano_hnl_cfg.py isMC=true
from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms
from glob import glob

options = VarParsing('python')

options.register('isMC'           ,  True           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run this on real data"                  )
options.register('skipDuplicated' ,  True           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Skip duplicated events. True by default")
options.register('globalTag'      , 'NOTSET'        , VarParsing.multiplicity.singleton, VarParsing.varType.string, "Set global tag"                         )
options.register('wantSummary'    ,  True           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run this on real data"                  )
options.register('wantFullRECO'   ,  False          , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run this on real data"                  )
options.register('reportEvery'    ,  1000           , VarParsing.multiplicity.singleton, VarParsing.varType.int   , "report every N events"                  )
options.register('skip'           ,  0              , VarParsing.multiplicity.singleton, VarParsing.varType.int   , "skip first N events"                    )
options.register('inputFile'      , None            , VarParsing.multiplicity.singleton, VarParsing.varType.string, "inputFile name"                         )
options.register('outFile'        , False           , VarParsing.multiplicity.singleton, VarParsing.varType.string, "outputFile name"                        )


options.setDefault('maxEvents', -1)
options.setDefault('tag', '10215')
options.parseArguments()


globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
if options._beenSet['globalTag']:
	globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['BParkNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['BParkFullEvt', extension[options.isMC], options.tag])+'.root')


if not options.inputFiles:
	options.inputFiles = [
		#'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj999.root'
		#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0B14F6F4-5DD5-B04B-A7EF-7D2B80AD33FA.root'
		'/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/50002/519412AE-918D-C247-B180-2C4E2002A8D0.root'
	] if not options.isMC else \
	[#'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-20to30_EMEnriched_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/0D2BBEDD-FFDA-E843-8620-D9B51558138C.root'
		# '/store/mc/RunIIAutumn18MiniAOD/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStar_BParking_102X_upgrade2018_realistic_v15-v2/100000/9B6F2D16-5F32-1E4D-83BC-8588E08F82AC.root'
#	'/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/100000/005594DA-4AA0-3E48-A8C2-46DECDE2E925.root',

	'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj1.root',
		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj2.root',
		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj3.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj4.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj5.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj6.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj7.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj8.root',
#       #	'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj9.root',
#       #	'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj10.root',
#       #	'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj11.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj12.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj13.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj14.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj15.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj16.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj17.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj18.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj19.root',
#		'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj20.root',
#		        ' /store/mc/RunIIAutumn18MiniAOD/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/100000/0BF89559-F5F7-4D41-A1B6-4037F80E9A4A.root'
		# 'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj278.root',
		# 'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj45.root',
		# 'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BToDMuN_HalfMu_Halfe_mass3p0_ctau184p0/step4_nj898.root'
	]

annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BParkNANO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.BParkingNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
	"PoolSource",
	fileNames = cms.untracked.vstring(options.inputFiles) if not options.inputFile else cms.untracked.vstring('file:{}'.format(options.inputFile)),
	secondaryFileNames = cms.untracked.vstring(),
	skipEvents=cms.untracked.uint32(options.skip),
	duplicateCheckMode = cms.untracked.string('checkEachFile' if options.skipDuplicated else 'noDuplicateCheck'),
)

process.options = cms.untracked.PSet(
	wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
	annotation = cms.untracked.string(annotation),
	name = cms.untracked.string('Applications'),
	version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
	dataset = cms.untracked.PSet(
	dataTier = cms.untracked.string('GEN-SIM-RECO'),
	filterName = cms.untracked.string('')
),
fileName = outputFileFEVT,
outputCommands = (cms.untracked.vstring('keep *',
		'drop *_*_SelectedTransient*_*',
		)),
splitLevel = cms.untracked.int32(0)
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
	compressionAlgorithm = cms.untracked.string('LZMA'),
	compressionLevel = cms.untracked.int32(9),
	dataset = cms.untracked.PSet(
	dataTier = cms.untracked.string('NANOAOD'),
	filterName = cms.untracked.string('')
),
#fileName = outputFileNANO,
	fileName = outputFileNANO if not options.outFile else cms.untracked.string('file:{}'.format(options.outFile)),
		outputCommands = cms.untracked.vstring(
		'drop *',
		'keep nanoaodFlatTable_*Table_*_*',     # event data
		'keep nanoaodUniqueString_nanoMetadata_*_*',   # basic metadata
		'keep nanoaodMergeableCounterTable_*Table_*_*',  # includes gentables
		)

)


# Additional output definition

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')
# this is for the LowPt energy regression
process.GlobalTag.toGet = cms.VPSet(
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
		label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To20_mean"),
		tag = cms.string("lowPtElectron_eb_ecalOnly_05To20_mean_2018V1"),
		connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
		label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To20_mean"),
		tag = cms.string("lowPtElectron_ee_ecalOnly_05To20_mean_2018V1"),
		connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
		label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To20_sigma"),
		tag = cms.string("lowPtElectron_eb_ecalOnly_05To20_sigma_2018V1"),
		connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
		label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To20_sigma"),
		tag = cms.string("lowPtElectron_ee_ecalOnly_05To20_sigma_2018V1"),
		connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
		label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To20_mean"),
		tag = cms.string("lowPtElectron_eb_ecalTrk_05To20_mean_2018V1"),
		connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To20_mean"),
			tag = cms.string("lowPtElectron_ee_ecalTrk_05To20_mean_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To20_sigma"),
			tag = cms.string("lowPtElectron_eb_ecalTrk_05To20_sigma_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To20_sigma"),
			tag = cms.string("lowPtElectron_ee_ecalTrk_05To20_sigma_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_eb_ecalOnly_20To50_mean"),
			tag = cms.string("lowPtElectron_eb_ecalOnly_20To50_mean_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_ee_ecalOnly_20To50_mean"),
			tag = cms.string("lowPtElectron_ee_ecalOnly_20To50_mean_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_eb_ecalOnly_20To50_sigma"),
			tag = cms.string("lowPtElectron_eb_ecalOnly_20To50_sigma_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_ee_ecalOnly_20To50_sigma"),
			tag = cms.string("lowPtElectron_ee_ecalOnly_20To50_sigma_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_eb_ecalTrk_20To50_mean"),
			tag = cms.string("lowPtElectron_eb_ecalTrk_20To50_mean_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_ee_ecalTrk_20To50_mean"),
			tag = cms.string("lowPtElectron_ee_ecalTrk_20To50_mean_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_eb_ecalTrk_20To50_sigma"),
			tag = cms.string("lowPtElectron_eb_ecalTrk_20To50_sigma_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("lowPtElectron_ee_ecalTrk_20To50_sigma"),
			tag = cms.string("lowPtElectron_ee_ecalTrk_20To50_sigma_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_mean"),
			tag = cms.string("gsfElectron_eb_ecalOnly_05To50_mean_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_mean"),
			tag = cms.string("gsfElectron_ee_ecalOnly_05To50_mean_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
			tag = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
			tag = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_mean"),
			tag = cms.string("gsfElectron_eb_ecalTrk_05To50_mean_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_mean"),
			tag = cms.string("gsfElectron_ee_ecalTrk_05To50_mean_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_sigma"),
			tag = cms.string("gsfElectron_eb_ecalTrk_05To50_sigma_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
	cms.PSet(record = cms.string("GBRDWrapperRcd"),
			label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_sigma"),
			tag = cms.string("gsfElectron_ee_ecalTrk_05To50_sigma_2018V1"),
			connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")))






# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.	GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

from PhysicsTools.BParkingNano.nanoBPark_cff import*# nanoAOD_customizeMuonTriggerBPark, nanoAOD_customizeTrackFilteredBPark,nanoAOD_customizeBToKMuMu
process = nanoAOD_customizeMuonTriggerBPark  (process)
process = nanoAOD_customizeTrackFilteredBPark(process)
process = nanoAOD_customizeElectronFilteredBPark(process)
process = nanoAOD_customizeBToMuLPi         (process, isMC=options.isMC)
#process = nanoAOD_customizeBToMuLPi_HD         (process, isMC=options.isMC)
process = nanoAOD_customizeBToKMuMu          (process, isMC=options.isMC) 
process = nanoAOD_customizeTriggerBitsBPark  (process)

# Path and EndPath definitions
process.nanoAOD_MuMuPi_step = cms.Path(process.nanoSequence + process.nanoBMuMuPiSequence  ) #candidate counter has been included in the sequence definition
process.nanoAOD_MuTrkPi_step = cms.Path(process.nanoSequence + process.nanoBMuEPiHDSequence  ) #candidate counter has been included in the sequence definition
process.nanoAOD_MuEPi_step = cms.Path(process.nanoSequence + process.nanoeSequence + process.nanoBMuEPiSequence  ) #candidate counter has been included in the sequence definition
#process.nanoAOD_KMuMu_step  = cms.Path(process.nanoSequence + process.nanoBKMuMuSequence + CountBToKmumu ) 

# customisation of the process.
if options.isMC:
	from PhysicsTools.BParkingNano.nanoBPark_cff import nanoAOD_customizeMC
	nanoAOD_customizeMC(process, ancestor_particles=[511, 521, 531, 541, 9900015]) 

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
	process.nanoAOD_MuMuPi_step,
        process.nanoAOD_MuEPi_step,
 #       process.nanoAOD_MuTrkPi_step,
 #       process.nanoAOD_KMuMu_step, 
	process.endjob_step, 
	process.NANOAODoutput_step
	)

if options.wantFullRECO:
	process.schedule = cms.Schedule(
 		process.nanoAOD_MuMuPi_step,
		process.nanoAOD_MuEPi_step,
 #               process.nanoAOD_MuTrkPi_step,
		#       process.nanoAOD_KMuMu_step, 
		process.endjob_step, 
		process.FEVTDEBUGHLToutput_step, 
		process.NANOAODoutput_step
		)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
	SelectEvents = cms.vstring('nanoAOD_MuMuPi_step','nanoAOD_MuEPi_step') 
#	SelectEvents = cms.vstring('nanoAOD_MuMuPi_step','nanoAOD_MuEPi_step','nanoAOD_MuTrkPi_step') 
#	SelectEvents = cms.vstring('nanoAOD_MuTrkPi_step') 
	)


### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
