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
config.General.requestName = "privateHNL_nano_Unresolved_Mass3_%s" % production_tag
config.General.workArea = 'crab_privateMCnano_Unresolved_Mass3'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
#config.JobType.pluginName = 'PrivateMC'
config.JobType.pluginName = 'Analysis'
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
config.Data.outputPrimaryDataset = 'private_BtoDMuN_NtoL_nano_Unresolved_Mass3'
config.Data.userInputFiles = open("nanoIn_eMu.txt").readlines() #implementing input through list of file: works wherever when path file is wired through xrootd  
#config.Data.splitting = 'EventBased'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.unitsPerJob = 100000
#config.Data.totalUnits = 50000000
config.Data.publication = False
config.Data.outputDatasetTag = 'privateHNLnano_%s'% production_tag
#config.Data.inputDBS = 'phys03'
## T3 Beijing
#config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/user/ratramon/HNLGen_ntuples/' #change this to output file: if output is to be sent to a local tier registered on xrootd, the path has to be in the form "/store/.."
## T3 Beijing

config.section_("Data")
config.Site.whitelist = ['T2_IT_Rome'] # to be changed accordingly to input site, jobs work way faster if they run on the same site where the input is stored
config.Site.storageSite = 'T2_IT_Rome' #this points to the storage site, for me it was Rome tier 2 


print(config)

