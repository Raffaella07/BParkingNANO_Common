import sys 
import os
import os.path
from os import path
import glob
import ROOT
from nanoLauncher import NanoLauncher


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to merge the nanoAOD files resulting from a multijob production', add_help=True)
  parser.add_argument('--pl' , type=str, dest='pl'       , help='label of the sample filei, e.g ParkingBPH4_Run2018B_V00'                     , default='V01_n9000000_njt300')
  parser.add_argument('--tag', type=str, dest='tag'      , help='[optional] tag to be added on the outputfile name'                           , default=None)
  parser.add_argument('--mcprivate'    , dest='mcprivate', help='run the BParking nano tool on a private MC sample'      , action='store_true', default=False)
  parser.add_argument('--mccentral'    , dest='mccentral', help='run the BParking nano tool on a central MC sample'      , action='store_true', default=False)
  parser.add_argument('--data'         , dest='data'     , help='run the BParking nano tool on a data sample'            , action='store_true', default=False)
  parser.add_argument('--doflat'       , dest='doflat'   , help='merge as well the files coming out of the truth-mathing', action='store_true', default=False)
  return parser.parse_args()


def checkParser(opt):

  if opt.pl==None:
    raise RuntimeError('Please indicate the production label of the nanoAOD production')

  if opt.mcprivate==False and opt.mccentral==False and opt.data==False:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding either --data or--mcprivate or --mccentral to the command line')

  if opt.mcprivate + opt.mccentral + opt.data > 1:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding only --data or --mcprivate or --mccentral to the command line')


class NanoMerger(NanoLauncher):
  def __init__(self, opt):
    self.prodlabel = vars(opt)['pl']
    self.tag       = vars(opt)['tag']
    self.mcprivate = vars(opt)['mcprivate']
    self.mccentral = vars(opt)['mccentral']
    self.data      = vars(opt)['data']
    self.doflat    = vars(opt)["doflat"]


  def doMerging(self, nanoName, mergedName, locationSE, outputdir, cond):

    print '\n-> Getting the different subdirectories (chunk/signal points)'
    subdirs = [f for f in glob.glob(locationSE+'/*')]

    for subdir in subdirs:
      if 'merged' in subdir: continue
      if '.root' in subdir: continue

      print '\n-> Processing: {}'.format(subdir[subdir.rfind('/')+1:len(subdir)])

      # get files
      nanoFiles = [f for f in glob.glob(subdir+nanoName)]

      # create the outputdir that will contain the merged file
      if not path.exists(subdir+outputdir):
        os.system('mkdir {}'.format(subdir+outputdir))

      command = 'python haddnano.py {}/{}'.format(subdir+outputdir, mergedName)

      print "\n-> Checking the files"
      for iFile, fileName in enumerate(nanoFiles):
        print fileName
        if iFile%100 == 0:
          print '     --> checked {}% of the files'.format(round(float(iFile)/len(nanoFiles)*100, 1))
        rootFile = ROOT.TNetXNGFile.Open(fileName, 'r')
        if not rootFile: continue
        else:
          if cond:
            if not rootFile.GetListOfKeys().Contains('Events'): continue

        command = command + ' {}'.format(fileName)

      print '\n-> Start of the merge'
      os.system(command)

      print '{}/{} created \n'.format(subdir+outputdir, mergedName)


  def doChunkMerging(self, nanoName, mergedName, locationSE):

    print '\n---> Merging the different chunks'

    nanoFiles = [f for f in glob.glob(locationSE+'/Chunk*/'+nanoName)]

    # create the outputdir that will contain the merged file
    if not path.exists('{}/merged'.format(locationSE)):
      os.system('mkdir {}/merged'.format(locationSE))

    filesValid = []
    print "\n-> Checking the files"
    for fileName in nanoFiles:
      rootFile = ROOT.TNetXNGFile.Open(fileName, 'r')
      if rootFile and rootFile.GetListOfKeys().Contains('Events'):
        filesValid.append(fileName)

    print '\n-> Start of the merge'
    command = 'python haddnano.py {}/merged/{}'.format(locationSE, mergedName)
    for iFile, fileName in enumerate(filesValid):
       command = command + ' {}'.format(fileName)

    os.system(command)

    print '{}/merged/{} created \n'.format(locationSE, mergedName)

    # clean
    for f in glob.glob(locationSE+'/Chunk*/merged/'):
      command_clean_file = 'rm -rf {}/bparknano.root'.format(f)
      os.system(command_clean_file)


  def process(self):
    print '---------------------------------'
    print '           Nano Merger           '
    print '---------------------------------'

    user      = os.environ["USER"]

    nanoName   = '/bparknano_nj*.root' if self.tag == None else '/bparknano_{}_nj*.root'.format(self.tag)
    mergedName = 'bparknano.root' if self.tag == None else 'bparknano_{}.root'.format(self.tag)
    outputdir  = '/merged'
        
    nanoName_tot   = 'merged/bparknano.root' if self.tag == None else 'merged/bparknano_{}.root'.format(self.tag)
    mergedName_tot = 'bparknano.root' if self.tag == None else 'bparknano_{}.root'.format(self.tag)

    if self.mcprivate:
      locationSE = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/'.format(user, self.prodlabel)
      
      pointdirs = NanoLauncher.getPointDirs(self, locationSE)

      for pointdir in pointdirs:
        if 'merged' in pointdir: continue

        print '\n --- Mass point: {} --- '.format(pointdir[pointdir.rfind('/')+1:len(pointdir)])

        # per chunk
        self.doMerging(nanoName, mergedName, pointdir+'/nanoFiles/', outputdir, True)

        # inclusive
        self.doChunkMerging(nanoName_tot, mergedName_tot, pointdir+'/nanoFiles/')
  
        if self.doflat:
          nanoName_flat   = '/nanoFiles/flat/flat_bparknano_nj*.root' if self.tag == None else '/nanoFiles/flat/flat_bparknano_{}_nj*.root'.format(self.tag)
          mergedName_flat = 'flat_bparknano.root' if self.tag == None else 'flat_bparknano_{}.root'.format(self.tag)
          outputdir_flat  = '/nanoFiles/flat/merged'

          self.doMerging(nanoName_flat, mergedName_flat, locationSE, outputdir_flat, False)

    elif self.data or self.mccentral:
      locationSE = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/{}'.format(user, 'data' if self.data else 'mc_central', self.prodlabel)
  
      # per chunk
      self.doMerging(nanoName, mergedName, locationSE, outputdir, True)

      # inclusive
      self.doChunkMerging(nanoName_tot, mergedName_tot, locationSE)


if __name__ == "__main__":

  opt = getOptions()

  checkParser(opt)

  NanoMerger(opt).process()

  print 'Done'


