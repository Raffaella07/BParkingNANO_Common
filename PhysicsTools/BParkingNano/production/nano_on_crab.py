#Very simple file for CRAB production of NANOAOD thtough BParkingNANO
from CRABClient.UserUtilities import config, ClientException
#from input_crab_data import dataset_files
import yaml
import random
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser



random.seed(2)
production_tag = datetime.date.today().strftime('%Y%b%d')

config = config()
config.section_("General")
config.General.requestName = "privateHNL_nano_part1D_%s" % production_tag
config.General.workArea = 'crab_privateMCnano_part1D'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
#config.JobType.pluginName = 'PrivateMC'
config.JobType.pluginName = 'Analysis'
globaltag = "102X_dataRun2_v11"
isMC = "0"
config.JobType.pyCfgParams = [
    'isMC=%s' % isMC, 'reportEvery=1000',
    'tag=%s' % production_tag,
    'globalTag=%s' % globaltag,
]
config.JobType.maxJobRuntimeMin = 3000
config.JobType.allowUndistributedCMSSW = True
#config.JobType.eventsPerLumi=100
config.JobType.psetName = '/afs/cern.ch/work/r/ratramon/HNL/CommonFrame/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/test/run_nano_HNLToL_cfg.py'
#config.JobType.pyCfgParams = [ 'nThr=1', 'inputFile=BPH-step2_.root', 'outputFile=BPH-step3.root','seedOffset=%i'% random.randint(0,1000)]
config.JobType.disableAutomaticOutputCollection = False
#config.JobType.scriptExe = '/afs/cern.ch/work/r/ratramon/HNL/CMSSW_10_2_15/src/HNLsGen/slurm/V00_v00_n10_njt1/slurm_mass1.0_ctau107.535167851_prod.sh'
#config.JobType.outputFiles = ['MiniAOD.root','GenSimAODSim_step1.log', 'GenSimAODSim_step2.log', 'FrameworkJobReport.xml', 'job.log']
config.JobType.inputFiles = [
#'/afs/cern.ch/work/r/ratramon/HNL/CMSSW_10_2_15/src/HNLsGen/slurm/V00_v00_n10_njt1/BPH_mod_BtoDMuN_NtoE_cfg.py','/afs/cern.ch/work/r/ratramon/HNL/CMSSW_10_2_15/src/HNLsGen/slurm/V00_v00_n10_njt1/step2.py' ,'/afs/cern.ch/work/r/ratramon/HNL/CMSSW_10_2_15/src/HNLsGen/slurm/V00_v00_n10_njt1/step3.py','/afs/cern.ch/work/r/ratramon/HNL/CMSSW_10_2_15/src/HNLsGen/slurm/V00_v00_n10_njt1/step4.py',
'/afs/cern.ch/work/r/ratramon/HNL/CommonFrame/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/test/lowPtEleReg_2018_02062020_nv.db' ]
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.lumiMask = 'LumiMask.txt'
#config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
config.Data.runRange = '320674-321134'
config.Data.splitting = 'LumiBased'
#config.Data.outputPrimaryDataset = 'private_BtoDMuN_NtoL_nano_part1D'
config.Data.inputDataset = '/ParkingBPH1/Run2018D-05May2019promptD-v1/MINIAOD' #implementing input through list of file: works wherever when path file is wired through xrootd  
#config.Data.userInputFiles = open("Part1D_mini_12Files.txt").readlines() #implementing input through list of file: works wherever when path file is wired through xrootd  
#config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 2
#config.Data.unitsPerJob = 100000
#config.Data.totalUnits = 50000000
config.Data.publication = False
config.Data.outputDatasetTag = 'BParkingNano_HNL_%s'% production_tag
#config.Data.inputDBS = 'phys03'
## T3 Beijing
#config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/user/ratramon/HNLGen_ntuples/' #change this to output file: if output is to be sent to a local tier registered on xrootd, the path has to be in the form "/store/.."
## T3 Beijing

config.section_("Data")
#config.Site.ignoreGlobalBlacklist = True
#config.Site.whitelist = ['T2_IT_Rome'] # to be changed accordingly to input site, jobs work way faster if they run on the same site where the input is stored
config.Site.storageSite = 'T2_IT_Rome' #this points to the storage site, for me it was Rome tier 2 

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

print(config)

