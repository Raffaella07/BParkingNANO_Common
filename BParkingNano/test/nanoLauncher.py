import sys 
import os
import glob
import ROOT


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to launch the nanoAOD tool on top of privately produced miniAOD files', add_help=True)
  parser.add_argument('--pl', type=str, dest='pl', help='label of the sample file', default='V01_n9000000_njt300')
  parser.add_argument('--tag', type=str, dest='tag', help='[optional] tag to be added on the outputfile name', default=None)
  # add isMC, maxEvents
  return parser.parse_args()


class NanoLauncher(object):
  def __init__(self, opt):
    self.prodlabel = vars(opt)['pl']
    self.tag       = vars(opt)['tag']
    self.user      = os.environ["USER"]


  def getPointDirs(self, location):
    return [f for f in glob.glob(location+'/*')]


  def getFiles(self, pointdir):
    return [f for f in glob.glob(pointdir+'/step4_nj*.root')]

    
  def checkFile(self, nanofile):
    rootFile = ROOT.TNetXNGFile.Open(nanofile, 'r')
    if rootFile and rootFile.GetListOfKeys().Contains('Events'):
      return True
    else: return False


  def createOutputDir(self, point='.'):
    outputdir = point + '/nanoFiles/'
    print 'outputdir ',outputdir
    os.system('mkdir {}'.format(outputdir))


  def getStep(self, nanofile):
    return nanofile[nanofile.rfind('nj')+2 :  nanofile.find('.root', 0)]


  def getOutputName(self, nanofile):
    outname = nanofile[0 : nanofile.find('step4', 0)] + 'nanoFiles/bparknano'
    if self.tag != None:
      outname += '_{}'.format(self.tag)
    return outname + '_nj' + self.getStep(nanofile) + '.root'


  def compile(self):
    import subprocess
    cwd = os.getcwd()
    loc = cwd[0:cwd.find('PhysicsTools')]
    subprocess.call("scram b", cwd=loc, shell=True)


  def launchNano(self, nanofile):
    command = 'sbatch -p wn  --account=t3 --time=00:40:00 -o logs/{pl}/nanostep_nj{nj}.log -e logs/{pl}/nanostep_nj{nj}.log --job-name=nanostep_nj{nj}_{pl} submitter.sh {infile} {outfile} {usr} {pl} {step}'.format(
      pl      = self.prodlabel,
      nj      = self.getStep(nanofile),
      infile  = nanofile, 
      outfile = self.getOutputName(nanofile),
      usr     = self.user, 
      step    = self.getStep(nanofile)
    )

    print '\n launching production over {}'.format(nanofile)
    os.system(command)


  def process(self):
  
    print '-> Compiling'
    self.compile()

    print '\n------------'
    print ' Processing NanoLauncher on production {} '.format(self.prodlabel)
    print ' --> on the batch'
    print '------------'

    locationSE = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/'.format(self.user, self.prodlabel)
    
    print '\n-> Creating log directory'
    logdir = './logs/{}'.format(self.prodlabel)
    os.system('mkdir -p {}'.format(logdir))

    print '\n-> Getting the different mass points'
    pointsdir = self.getPointDirs(locationSE)
    
    # looping over all the different points
    for point in pointsdir:
   
      print '\n-> Processing mass/ctau point: {}'.format(point[point.rfind('/')+1:len(point)])

      print '\n  --> Creating output directory'
      self.createOutputDir(point)

      print '\n  --> Fetching the files '.format(point)
      nanofiles = self.getFiles(point)

      for nanofile in nanofiles:

        if self.checkFile(nanofile):
          self.launchNano(nanofile)
        else:
          print '    could not open {} --> skipping'.format(nanofile)

    print '\n-> Submission completed' 

    

if __name__ == "__main__":

  opt = getOptions()

  NanoLauncher(opt).process()




