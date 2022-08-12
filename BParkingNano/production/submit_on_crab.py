from CRABClient.UserUtilities import config, ClientException, getUsernameFromCRIC
#from input_crab_data i  mport dataset_files
import yaml        
import numpy as np      
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser
from FWCore.PythonUtilities.LumiList import LumiList
import CRABClient
from dbs.apis.dbsClient import DbsApi
dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')


production_tag = datetime.date.today().strftime('%Y%b%d')

config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'BParkingNANO_fast_%s' % production_tag


config.section_('Data')
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/ratramon/HNLGen_ntuples/' #change this to output file: if output is to be sent to a local tier registered on xrootd, the path has to be in the form "/store/.."
config.Data.inputDBS = 'global'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/run_nano_HNLToL_cfg.py'
config.JobType.maxJobRuntimeMin = 3000
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ["../test/lowPtEleReg_2018_02062020_nv.db"]

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Rome'
config.Site.blacklist = ['T2_US_Caltech']
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
if __name__ == '__main__':

  from CRABAPI.RawCommand import crabCommand
  from CRABClient.ClientExceptions import ClientException
  from httplib import HTTPException
  from multiprocessing import Process

  def submit(config):
      try:
          crabCommand('submit', config = config)
      except HTTPException as hte:
          print "Failed submitting task: %s" % (hte.headers)
      except ClientException as cle:
          print "Failed submitting task: %s" % (cle)


  parser = ArgumentParser()
  parser.add_argument('-y', '--yaml', default = 'samples.yml', help = 'File with dataset descriptions')
  parser.add_argument('-f', '--filter', default='*', help = 'filter samples, POSIX regular expressions allowed')
  args = parser.parse_args()

  with open(args.yaml) as f:
    doc = yaml.load(f) # Parse YAML file
    common = doc['common'] if 'common' in doc else {'data' : {}, 'mc' : {}}
    
    # loop over samples
    for sample, info in doc['samples'].iteritems():
        print(info)
      # Given we have repeated datasets check for different parts
        parts = info['parts'] if 'parts' in info else [None]
    	for part in parts:
          name = sample % part if part is not None else sample 
        # filter names according to what we need
          if not fnmatch(name, args.filter):
		 continue
	  
	  isMC = info['isMC']
          common_branch = 'mc' if isMC else 'data'
	  config.Data.unitsPerJob = info.get(
	      'splitting',
	      common[common_branch].get('splitting', None)
	  )
	  globaltag = info.get(
	      'globaltag',
	      common[common_branch].get('globaltag', None)
	  )
	  
	  config.JobType.pyCfgParams = [
	      'isMC=%s' % isMC, 'reportEvery=1000',
	      'tag=%s' % production_tag,
	      'globalTag=%s' % globaltag,
	  ]
          config.Data.splitting = 'FileBased' if isMC else 'Automatic'
	          
	  config.JobType.outputFiles = ['_'.join(['BParkNANO', 'mc' if isMC else 'data', production_tag])+'.root']
	  
	  config.Data.inputDataset = info['dataset'] % part \
	                             if part is not None else \
	                                info['dataset']
	  runs= [run['run_num'] for run in dbs.listRuns(dataset = config.Data.inputDataset)]
	     # print(runs)
          splits = np.array_split(runs,4)
	  print(splits)
	  if not isMC:
	 	 for i in range(1,len(splits)):
	  	
		          config.General.requestName = name+"_section"+str(i) if len(splits)!=0 else sample
		          print 'submitting', config.General.requestName
			  lumiList = LumiList(filename=info.get(
		             'lumimask', 
		              common[common_branch].get('lumimask', None))
		       		)
			      
		          lumiList.selectRuns([x for x in splits[i]])
			  lumiList.writeJSON('lumi_mask_task_'+str(i)+'.json')
		          config.Data.lumiMask = 'lumi_mask_task_'+str(i)+'.json'
		              	
		  
		          
		          print config
	 	          submit(config)
	  else:	
		config.General.requestName = name	
		print 'submitting', config.General.requestName
		common_branch = 'mc' if isMC else 'data'
		config.Data.lumiMask = ''
		print config
	 	submit(config)
  
