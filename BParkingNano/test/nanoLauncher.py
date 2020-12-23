import sys 
import os
import glob
import ROOT


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to launch the nanoAOD tool on top of privately produced miniAOD files', add_help=True)
  parser.add_argument('--pl'      , type=str, dest='pl'       , help='label of the sample file'                                                       , default=None)
  parser.add_argument('--ds'      , type=str, dest='ds'       , help='[data-mccentral] name of the dataset, e.g /ParkingBPH4/Run2018B-05May2019-v2/MINIAOD', default=None)
  parser.add_argument('--tag'     , type=str, dest='tag'      , help='[optional] tag to be added on the outputfile name'                              , default=None)
  parser.add_argument('--maxfiles', type=str, dest='maxfiles' , help='[optional] maximum number of files to process'                                  , default=None)
  parser.add_argument('--user'    , type=str, dest='user'     , help='[optional-mcprivate] specify username where the miniAOD files are stored'       , default=os.environ["USER"])
  parser.add_argument('--mcprivate'         , dest='mcprivate', help='run the BParking nano tool on a private MC sample'         , action='store_true', default=False)
  parser.add_argument('--mccentral'         , dest='mccentral', help='run the BParking nano tool on a central MC sample'         , action='store_true', default=False)
  parser.add_argument('--data'              , dest='data'     , help='run the BParking nano tool on a data sample'               , action='store_true', default=False)
  parser.add_argument('--doflat'            , dest='doflat'   , help='[optional-mcprivate] launch the ntupliser on top of the nanofile' , action='store_true', default=False)
  return parser.parse_args()


def checkParser(opt):

  if opt.pl==None:
    raise RuntimeError('Please indicate the production label: for --mcprivate, it has to correspond to the label of the miniAODl (eg. V01_n9000000_njt300)')

  if opt.data==True and opt.ds==None:
    raise RuntimeError('You are running on data, please indicate the dataset with --ds <dataset>')

  if opt.mcprivate==False and opt.mccentral==False and opt.data==False:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding either --data or--mcprivate or --mccentral to the command line')

  if opt.mcprivate + opt.mccentral + opt.data > 1:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding only --data or --mcprivate or --mccentral to the command line')


class NanoLauncher(object):
  def __init__(self, opt):
    self.prodlabel = vars(opt)['pl']
    self.dataset   = vars(opt)['ds']
    self.tag       = vars(opt)['tag']
    self.maxfiles  = vars(opt)['maxfiles']
    self.mcprivate = vars(opt)['mcprivate']
    self.mccentral = vars(opt)['mccentral']
    self.data      = vars(opt)['data']
    self.user      = vars(opt)["user"]
    self.doflat    = vars(opt)["doflat"]


  def getPointDirs(self, location):
    return [f for f in glob.glob(location+'/*')]


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


  def writeFileList(self, point=None):

    os.system('mkdir ./files')

    if self.mcprivate:
      os.system('rm ./files/filelist_{}_{}.txt'.format(self.prodlabel, point))

      myfile = open("./files/filelist_{}_{}.txt".format(self.prodlabel, point), "w+")
      
      nanofiles = self.getLocalFiles(point)

      for iNano, nanofile in enumerate(nanofiles):

        if self.checkLocalFile(nanofile):
          iNano += 1
          myfile.write(nanofile + '\n')
        else:
          print '    could not open {} --> skipping'.format(nanofile)
    
      myfile.close()  
      
      if self.maxfiles!=None:
        command_cut = 'head -n {max} ./files/filelist_{pl}_{p}.txt > ./files/tmp.txt && mv ./files/tmp.txt ./files/filelist_{pl}_{p}.txt'.format(max=self.maxfiles, pl=self.prodlabel, p=point)
        os.system(command_cut)
      
      print '    ---> ./files/filelist_{pl}_{p}.txt created'.format(pl=self.prodlabel, p=point)


    elif self.data or self.mccentral:
      ds_label = self.getDataLabel() if self.data else self.getMCLabel()

      os.system('rm ./files/filelist_{ds}_{pl}*.txt'.format(ds=ds_label, pl=self.prodlabel))
      
      command = 'dasgoclient --query="file dataset={ds} | grep file.name" > ./files/filelist_{dsl}_{pl}.txt'.format(ds=self.dataset, dsl=ds_label, pl=self.prodlabel)
      os.system(command)

      if self.maxfiles!=None:
        command_cut = 'head -n {max} ./files/filelist_{ds}_{pl}.txt > ./files/tmp.txt && mv ./files/tmp.txt ./files/filelist_{ds}_{pl}.txt'.format(max=self.maxfiles, ds=ds_label, pl=self.prodlabel)
        os.system(command_cut)

      # slurm cannot deal with too large arrays
      # -> submit job arrays of size 1000
      if self.getSize('./files/filelist_{}_{}.txt'.format(ds_label, self.prodlabel)) > 1000:
        command_split = 'split -l 1000 ./files/filelist_{ds}_{pl}.txt ./files/filelist_{ds}_{pl} --additional-suffix=.txt'.format(ds=ds_label, pl=self.prodlabel)
        os.system(command_split)
        os.system('rm ./files/filelist_{}_{}.txt'.format(ds_label, self.prodlabel))

      print '    ---> ./files/filelist_{ds}_{pl}.txt created'.format(ds=ds_label, pl=self.prodlabel)
      

  def launchNano(self, nNano, outputdir, logdir, filelist, label):
    #command = 'sbatch -p quick --account=t3 --time==00:50:00 -o logs/{pl}/nanostep_nj%a.log -e logs/{pl}/nanostep_nj%a.log --job-name=nanostep_nj%a_{pl} --array {ar} submitter.sh {outdir} {usr} {pl} {tag} {flt}'.format(
    command = 'sbatch -p wn --account=t3 -o {ld}/nanostep_nj%a.log -e {ld}/nanostep_nj%a.log --job-name=nanostep_nj%a_{pl} --array {ar} --time=03:00:00 submitter.sh {outdir} {usr} {pl} {tag} {isMC} {rmt} {flt} {lst}'.format(
      ld      = logdir,
      pl      = label,
      ar      = '1-{}'.format(nNano),
      outdir  = outputdir,
      usr     = os.environ["USER"], 
      tag     = 0 if self.tag == None else self.tag,
      isMC    = 1 if self.mcprivate or self.mccentral else 0,
      rmt     = 0 if self.mcprivate else 1,
      flt     = 1 if self.doflat == True else 0,
      lst     = filelist
    )

    os.system(command)


  def process(self):
  
    print '-> Compiling'
    self.compile()
    

    print '\n------------'
    print ' Processing NanoLauncher on production {} '.format(self.prodlabel if self.mcprivate else self.dataset)
    print ' --> on the batch'
    print '------------'


    if self.mcprivate:
      locationSE = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/'.format(self.user, self.prodlabel)

      print '\n-> Getting the different mass points'
      pointsdir = self.getPointDirs(locationSE)
      points    = [point[point.rfind('/')+1:len(point)] for point in pointsdir]
      
      # looping over all the different points
      for point in points:
     
        print '\n-> Processing mass/ctau point: {}'.format(point)

        print '\n  --> Creating output directory'
        outputdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/{}/nanoFiles/'.format(os.environ["USER"], self.prodlabel, point)
        os.system('mkdir -p {}'.format(outputdir))
        
        print '\n  --> Creating log directory'
        logdir = './logs/{}/{}'.format(self.prodlabel, point) if self.tag == None else './logs/{}/{}_{}'.format(self.prodlabel, point, self.tag)
        os.system('mkdir -p {}'.format(logdir))

        print '\n  --> Fetching the files '
        self.writeFileList(point)
            
        filelist = './files/filelist_{}_{}.txt'.format(self.prodlabel, point)
        label = self.prodlabel + '_' + point if self.tag == None else self.prodlabel + '_' + point + '_' +self.tag

        self.launchNano(self.getSize(filelist), outputdir, logdir, filelist, label) 

    
    elif self.data or self.mccentral:
      if self.mccentral:
        if 'ext' in self.dataset:
          self.prodlabel += '_ext'
      
      print '\n  --> Fetching the files '
      self.writeFileList()

      ds_label = self.getDataLabel() if self.data else self.getMCLabel()

      # loop on the files (containing at most 1k samples) 
      for iFile, filelist in enumerate(glob.glob('./files/filelist_{}_{}*.txt'.format(ds_label, self.prodlabel))):

        print '\n  --> Creating output directory'
        dirname = 'data' if self.data else 'mc_central'
        outputdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/{}_{}/Chunk{}_n{}/'.format(os.environ["USER"], dirname, ds_label, self.prodlabel, iFile, self.getSize(filelist))
        os.system('mkdir -p {}'.format(outputdir))
        
        print '\n-> Creating log directory'
        logdir = './logs/{}/{}/Chunk{}_n{}'.format(ds_label, self.prodlabel, iFile, self.getSize(filelist)) if self.tag == None \
               else './logs/{}/{}_{}/Chunk{}_n{}'.format(ds_label, self.prodlabel, self.tag, iFile, self.getSize(filelist))
        os.system('mkdir -p {}'.format(logdir))
          
        label = '{}_{}_Chunk{}_n{}'.format(ds_label, self.prodlabel, iFile, self.getSize(filelist))

        self.launchNano(self.getSize(filelist), outputdir, logdir, filelist, label)


    print '\n-> Submission completed' 

    

if __name__ == "__main__":

  opt = getOptions()

  checkParser(opt)

  NanoLauncher(opt).process()




