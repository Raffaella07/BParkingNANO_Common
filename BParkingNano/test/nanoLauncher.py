import sys 
import os
import os.path
from os import path
import glob
import subprocess
import ROOT

from nanoTools import NanoTools


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to launch the nanoAOD tool on top of privately produced miniAOD files', add_help=True)
  parser.add_argument('--pl'      , type=str, dest='pl'         , help='label of the sample file'                                                            , default=None)
  parser.add_argument('--ds'      , type=str, dest='ds'         , help='[data-mccentral] name of the dataset, e.g /ParkingBPH4/Run2018B-05May2019-v2/MINIAOD', default=None)
  parser.add_argument('--tag'     , type=str, dest='tag'        , help='[optional] tag to be added on the outputfile name'                                   , default=None)
  parser.add_argument('--maxfiles', type=str, dest='maxfiles'   , help='[optional] maximum number of files to process'                                       , default=None)
  parser.add_argument('--user'    , type=str, dest='user'       , help='[optional-mcprivate] specify username where the miniAOD files are stored'            , default=os.environ["USER"])
  parser.add_argument('--mcprivate'         , dest='mcprivate'  , help='run the BParking nano tool on a private MC sample'              , action='store_true', default=False)
  parser.add_argument('--mccentral'         , dest='mccentral'  , help='run the BParking nano tool on a central MC sample'              , action='store_true', default=False)
  parser.add_argument('--data'              , dest='data'       , help='run the BParking nano tool on a data sample'                    , action='store_true', default=False)
  parser.add_argument('--donano'            , dest='donano'     , help='launch the nano tool on top of the minifile'                    , action='store_true', default=False)
  parser.add_argument('--doflat'            , dest='doflat'     , help='launch the ntupliser on top of the nanofile'                    , action='store_true', default=False)
  parser.add_argument('--domergenano'       , dest='domergenano', help='[optional] merge the nanofile steps'                            , action='store_true', default=False)
  parser.add_argument('--doquick'           , dest='doquick'    , help='[optional] run the jobs on the quick partition (t/job<1h)'      , action='store_true', default=False)
  parser.add_argument('--docompile'         , dest='docompile'  , help='[optional] compile the full BParkingNano tool'                  , action='store_true', default=False)
  return parser.parse_args()


def checkParser(opt):
  if opt.pl==None:
    raise RuntimeError('Please indicate the production label: for --mcprivate, it has to correspond to the label of the miniAODl (eg. V01_n9000000_njt300)')

  if opt.data==True and opt.ds==None:
    raise RuntimeError('You are running on data, please indicate the dataset with --ds <dataset>')

  if opt.donano==False and opt.doflat==False:
    raise RuntimeError('Please indicate if you want to run the nano tool (--donano) and/or the ntupliser (--doflat)')

  if opt.mcprivate==False and opt.mccentral==False and opt.data==False:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding either --data or--mcprivate or --mccentral to the command line')

  if opt.mcprivate + opt.mccentral + opt.data > 1:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding only --data or --mcprivate or --mccentral to the command line')


class NanoLauncher(NanoTools):
  def __init__(self, opt):
    self.prodlabel   = vars(opt)['pl']
    self.dataset     = vars(opt)['ds']
    self.tag         = vars(opt)['tag']
    self.maxfiles    = vars(opt)['maxfiles']
    self.mcprivate   = vars(opt)['mcprivate']
    self.mccentral   = vars(opt)['mccentral']
    self.data        = vars(opt)['data']
    self.user        = vars(opt)["user"]
    self.donano      = vars(opt)["donano"]
    self.doflat      = vars(opt)["doflat"]
    self.domergenano = vars(opt)["domergenano"]
    self.doquick     = vars(opt)["doquick"]
    self.docompile   = vars(opt)["docompile"]


  def getLocalFiles(self, point):
    pointdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/{}/'.format(self.user, self.prodlabel, point)
    return [f for f in glob.glob(pointdir+'/step4_nj*.root')]

    
  def checkLocalFile(self, nanofile):
    self.prodlabel = vars(opt)['pl']
    rootFile = ROOT.TNetXNGFile.Open(nanofile, 'r')
    if rootFile and rootFile.GetListOfKeys().Contains('Events'):
      return True
    else: return False


  def getDataLabel(self):
    idx = self.dataset.find('-')
    return self.dataset.replace('/', '_')[1:idx]
  

  def getMCLabel(self):
    idx = self.dataset[1:].find('/')
    return self.dataset[1:idx+1]


  def compile(self):
    import subprocess
    cwd = os.getcwd()
    loc = cwd[0:cwd.find('PhysicsTools')]
    subprocess.call("scram b", cwd=loc, shell=True)


  def getSize(self, file):
    return sum(1 for line in open(file))


  def getJobId(self, job):
    return int(job[job.find('job')+4:])
    #return str(job[job.find('job')+4:])


  def getJobIdsList(self, jobIds):
    listIds = ''
    for jobId in jobIds:
      listIds += '{}:'.format(jobId)
    return listIds[:len(listIds)-1]


  def writeFileList(self, nfiles_perchunk, point=None):
    if not path.exists('./files'):
      os.system('mkdir ./files') 

    if self.mcprivate:
      filename = './files/filelist_{}_{}'.format(self.prodlabel, point)
    else:
      ds_label = self.getDataLabel() if self.data else self.getMCLabel()
      filename = './files/filelist_{dsl}_{pl}'.format(dsl=ds_label, pl=self.prodlabel) 

      if len(glob.glob('{}*.txt'.format(filename))) == 0: # do not create file if already exists
        if self.mcprivate: # fetch miniAOD files
          myfile = open(filename + '.txt', "w+")
          nanofiles = self.getLocalFiles(point)

          for nanofile in nanofiles:
            if self.checkLocalFile(nanofile):
              myfile.write(nanofile + '\n')
            else:
              print '    could not open {} --> skipping'.format(nanofile)
        
          myfile.close()  
        else: # fetch files on DAS
          command = 'dasgoclient --query="file dataset={ds} | grep file.name" > {fn}.txt'.format(ds=self.dataset, fn=filename)
          os.system(command)

        # slurm cannot deal with too large arrays
        # -> submit job arrays of size 750
        if self.getSize(filename + '.txt') > nfiles_perchunk:
          command_split = 'split -l {nfiles} {fn}.txt {fn}_ --additional-suffix=.txt'.format(nfiles=nfiles_perchunk, fn=filename)
          os.system(command_split)
          os.system('rm {fn}.txt'.format(fn=filename))
        
    print '    ---> {}*.txt created'.format(filename)

    return filename 


  def writeMergerSubmitter(self, label, filetype):
    # defining the command
    #type_ = 'mcprivate' if self.mcprivate else 'data'

    prodlabel = ''
    if self.mcprivate: prodlabel = self.prodlabel
    elif self.mccentral: prodlabel = self.getMCLabel() + '_' + self.prodlabel 
    elif self.data: prodlabel = self.getDataLabel() + '_' + self.prodlabel 

    #command = 'python nanoMerger.py --dobatch --pl {pl} --{tp}'.format(
    command = 'python nanoMerger.py --dobatch --pl {pl}'.format(
        pl = prodlabel,
        #tp = type_,
        )
    if self.tag != None:
      command += ' --tag {}'.format(self.tag)
    if filetype == 'nano': command += ' --donano' 
    else: command += ' --doflat'
    
    if self.mcprivate: command += ' --mcprivate'
    if self.mccentral: command += ' --mccentral'
    if self.data: command += ' --data'
    # add for flat

    # defining the workdir
    dirlabel = label
    workdir = '/scratch/{usr}/mergingstep/{lbl}'.format(usr=os.environ["USER"], lbl=dirlabel)

    # content of the submitter
    content = [
      '#!/bin/bash',
      'workdir="{wrkdir}"',
      'mkdir -p $workdir',
      'cp nanoTools.py $workdir',
      'cp nanoMerger.py $workdir',
      'cp haddnano.py $workdir',
      'cd $workdir',
      '{cm}',
      'cd -',
      'rm -r $workdir',
    ]
    content = '\n'.join(content).format(
        wrkdir = workdir,
        cm     = command,
        )

    # create file
    submitter_merger = open('submitter_merger.sh', 'w+')
    submitter_merger.write(content)
    submitter_merger.close()


  def launchNano(self, nfiles, outputdir, logdir, filelist, label):
    if not self.doquick:
      slurm_options = '-p wn --account=t3 -o {ld}/nanostep_nj%a.log -e {ld}/nanostep_nj%a.log --job-name=nanostep_nj%a_{pl} --array {ar} --time=03:00:00'.format(
        ld = logdir,
        pl = label,
        ar = '1-{}'.format(nfiles),
        )
    else:
      slurm_options = '-p quick --account=t3 -o {ld}/nanostep_nj%a.log -e {ld}/nanostep_nj%a.log --job-name=nanostep_nj%a_{pl} --array {ar}'.format(
        ld = logdir,
        pl = label,
        ar = '1-{}'.format(nfiles),
        )

    command = 'sbatch {slurm_opt} submitter.sh {outdir} {usr} {pl} {tag} {isMC} {rmt} {lst} 0'.format(
      slurm_opt = slurm_options,
      pl        = label,
      outdir    = outputdir,
      usr       = os.environ["USER"], 
      tag       = 0 if self.tag == None else self.tag,
      isMC      = 1 if self.mcprivate or self.mccentral else 0,
      rmt       = 0 if self.mcprivate else 1,
      lst       = filelist,
      )

    job = subprocess.check_output(command, shell=True)
    print '\n       ---> {}'.format(job)

    return self.getJobId(job)


  def launchDumper(self, outputdir, logdir, label, jobId):
    if not self.doquick:
      slurm_options = '-p wn --account=t3 -o {ld}/dumperstep.log -e {ld}/dumperstep.log --job-name=dumperstep_{pl} --time=03:00:00 {dp}'.format(
        ld      = logdir,
        pl      = label,
        dp      = '--dependency=afterany:{}'.format(jobId) if jobId != -99 else '',
        )
    else:
      slurm_options = '-p quick --account=t3 -o {ld}/dumperstep.log -e {ld}/dumperstep.log --job-name=dumperstep_{pl} {dp}'.format(
        ld      = logdir,
        pl      = label,
        dp      = '--dependency=afterany:{}'.format(jobId) if jobId != -99 else '',
        )

    command = 'sbatch {slurm_opt} submitter_dumper.sh {outdir} {usr} {pl} {tag} {isMC}'.format(
      slurm_opt = slurm_options,
      pl      = label,
      outdir  = outputdir,
      usr     = os.environ["USER"], 
      tag     = 0 if self.tag == None else self.tag,
      isMC    = 1 if self.mcprivate or self.mccentral else 0,
      )
    
    job_dump = subprocess.check_output(command, shell=True)
    if jobId == -99:
      print '\n       ---> {}'.format(job_dump)
    else:
      print '       ---> (dependency)' 
      print '            {}'.format(job_dump)

    return self.getJobId(job_dump)


  def launchMerger(self, logdir, label, jobIds, filetype):
    self.writeMergerSubmitter(label, filetype)

    if not self.doquick:
      slurm_options = '-p wn --account=t3 -o {ld}/mergerstep.log -e {ld}/mergerstep.log --job-name=mergerstep_{pl} --time=02:00:00 --dependency=afterany:{jobid}'.format(
        ld    = logdir,
        pl    = label,
        jobid = self.getJobIdsList(jobIds),
        )
    else:
      slurm_options = '-p quick --account=t3 -o {ld}/mergerstep.log -e {ld}/mergerstep.log --job-name=mergerstep_{pl} --dependency=afterany:{jobid}'.format(
        ld    = logdir,
        pl    = label,
        jobid = self.getJobIdsList(jobIds),
        )

    command_merge = 'sbatch {slurm_opt} submitter_merger.sh'.format(
        slurm_opt = slurm_options,
        )
    
    n_job_merge = subprocess.check_output(command_merge, shell=True)
    print '       ---> (dependency)' 
    print '            {}'.format(n_job_merge)

    os.system('rm submitter_merger.sh')


  def launchingModule(self, point=None, ds_label=None):
    # declaring some useful quantities
    # ids of the launched jobs, needed for dependency of the merger
    nano_jobIds = []
    flat_jobIds = []

    # counter of processed files
    nfiles_tot = 0
        
    # slurm cannot deal with too large arrays, so does haddnano (keep it hardcoded)
    maxfiles_perchunk = 750
    
    print '\n  --> Fetching the files'
    filelistname = self.writeFileList(maxfiles_perchunk, point)

    # loop on the files (containing at most 750 samples) 
    for iFile, filelist in enumerate(glob.glob('{}*.txt'.format(filelistname))):
      if self.getSize(filelist) == 0:
        print '        WARNING: no files were found with the corresponding production label'
        print '                 Did you set the correct username using --user <username>?'

      # enforcing max files limit
      if self.maxfiles != None and nfiles_tot >= int(self.maxfiles): continue 

      # number of files to process in this chunk
      if self.maxfiles == None or int(self.maxfiles)-nfiles_tot > maxfiles_perchunk:
        nfiles = self.getSize(filelist)
      else:
        nfiles = int(self.maxfiles)-nfiles_tot
      nfiles_tot += nfiles
      
      # merging step (if any) must happen after nano or flat steps are completed
      if self.maxfiles == None:
        merge_cond = (iFile == len(glob.glob('{}*.txt'.format(filelistname)))-1)
      else:
        merge_cond = (nfiles_tot == int(self.maxfiles))
        
      print '\n  --> Creating output directory'
      outputdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/'.format(os.environ["USER"])
      if self.mcprivate:
        outputdir += '{}/{}/nanoFiles/Chunk{}_n{}'.format(self.prodlabel, point, iFile, nfiles)
      else:
        dirname = 'data' if self.data else 'mc_central'
        outputdir += '{}/{}_{}/Chunk{}_n{}/'.format(dirname, ds_label, self.prodlabel, iFile, nfiles)
      if not path.exists(outputdir):
        os.system('mkdir -p {}'.format(outputdir))
      
      print '\n  --> Creating log directory'
      label1 = self.prodlabel if self.mcprivate else ds_label
      label2 = point if self.mcprivate else self.prodlabel
      logdir = './logs/{}/{}/Chunk{}_n{}'.format(label1, label2, iFile, nfiles) if self.tag == None \
               else './logs/{}/{}_{}/Chunk{}_n{}'.format(label1, label2, self.tag, iFile, nfiles)
      if not path.exists(logdir):
        os.system('mkdir -p {}'.format(logdir))

      label = '{}_{}_Chunk{}_n{}'.format(label1, label2 if self.tag==None else label2+'_'+self.tag, iFile, self.getSize(filelist))

      nano_jobId = -99

      if self.donano:
        nano_jobId = self.launchNano(nfiles, outputdir, logdir, filelist, label)

        if self.domergenano:
          nano_jobIds.append(nano_jobId)

          if merge_cond:
            self.launchMerger(logdir, label, nano_jobIds, 'nano')

      if self.doflat:
        flat_outputdir = outputdir + '/flat'
        if not path.exists(flat_outputdir):
          os.system('mkdir -p {}'.format(flat_outputdir))

        flat_jobId = self.launchDumper(outputdir, logdir, label, nano_jobId)

        # merging of flat files happens automatically
        flat_jobIds.append(flat_jobId)

        if merge_cond:
          self.launchMerger(logdir, label, flat_jobIds, 'flat')


  def process(self):
    if self.docompile:
      print '-> Compiling'
      self.compile()

    print '\n------------'
    print ' Processing NanoLauncher on production {} '.format(self.prodlabel if self.mcprivate else self.dataset)
    print ' --> on the batch'
    print '------------'
    
    if self.mcprivate:
      locationSE = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/'.format(self.user, self.prodlabel)

      print '\n-> Getting the different mass points'
      pointsdir = NanoTools.getPointDirs(self, locationSE)
      points    = [point[point.rfind('/')+1:len(point)] for point in pointsdir]
      
      # looping over the signal points
      for point in points:
        print '\n-> Processing mass/ctau point: {}'.format(point)

        if point != 'mass3.0_ctau184.256851021': continue

        self.launchingModule(point=point)

    
    elif self.data or self.mccentral:
      if self.mccentral:
        if 'ext' in self.dataset:
          self.prodlabel += '_ext'

      dataset_label = self.getDataLabel() if self.data else self.getMCLabel()

      self.launchingModule(ds_label=dataset_label)

    print '\n-> Submission completed' 

    

if __name__ == "__main__":

  opt = getOptions()

  checkParser(opt)

  NanoLauncher(opt).process()




