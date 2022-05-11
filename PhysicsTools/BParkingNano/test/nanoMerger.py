import sys 
import os
import os.path
from os import path
import glob
import ROOT
from nanoTools import NanoTools


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to merge the nanoAOD files resulting from a multijob production', add_help=True)
  parser.add_argument('--pl'     , type=str, dest='pl'       , help='label of the nano sample production'                                    , default='V01_n9000000_njt300')
  parser.add_argument('--ds'     , type=str, dest='ds'       , help='[data-mccentral] name of the dataset'                                   , default=None)
  parser.add_argument('--tagnano', type=str, dest='tagnano'  , help='[optional] tag to be added on the outputfile name of the nano sample'   , default=None)
  parser.add_argument('--tagflat', type=str, dest='tagflat'  , help='[optional] tag to be added on the outputfile name of the flat sample'   , default=None)
  parser.add_argument('--mcprivate'        , dest='mcprivate', help='run the BParking nano tool on a private MC sample' , action='store_true', default=False)
  parser.add_argument('--mccentral'        , dest='mccentral', help='run the BParking nano tool on a central MC sample' , action='store_true', default=False)
  parser.add_argument('--data'             , dest='data'     , help='run the BParking nano tool on a data sample'       , action='store_true', default=False)
  parser.add_argument('--donano'           , dest='donano'   , help='merge nano files'                                  , action='store_true', default=False)
  parser.add_argument('--doflat'           , dest='doflat'   , help='merge flat files'                                  , action='store_true', default=False)
  parser.add_argument('--dobatch'          , dest='dobatch'  , help='to be turned on if running the script on the batch', action='store_true', default=False)
  return parser.parse_args()


def checkParser(opt):
  if opt.pl==None:
    raise RuntimeError('Please indicate the production label of the nanoAOD production')

  if opt.donano==False and opt.doflat==False:
    raise RuntimeError('Please indicate if you want to run the nano tool (--donano) and/or the ntupliser (--doflat)')

  if opt.mcprivate==False and opt.mccentral==False and opt.data==False:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding either --data or--mcprivate or --mccentral to the command line')

  if opt.mcprivate + opt.mccentral + opt.data > 1:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding only --data or --mcprivate or --mccentral to the command line')

  if opt.data and opt.ds==None:
    raise RuntimeError('Please indicate the dataset you want to run the tool on using --ds <dataset>')


class NanoMerger(NanoTools):
  def __init__(self, opt):
    self.prodlabel = vars(opt)['pl']
    self.dataset   = vars(opt)['ds']
    self.tagnano   = vars(opt)['tagnano']
    self.tagflat   = vars(opt)['tagflat']
    self.mcprivate = vars(opt)['mcprivate']
    self.mccentral = vars(opt)['mccentral']
    self.data      = vars(opt)['data']
    self.donano    = vars(opt)["donano"]
    self.doflat    = vars(opt)["doflat"]
    self.dobatch   = vars(opt)["dobatch"]


  def doMerging(self, nanoName, mergedName, locationSE, outputdir, cond):
    print '\n-> Getting the different subdirectories (chunk/signal points)'
    subdirs = [f for f in glob.glob(locationSE+'/*')]

    for subdir in subdirs:
      if 'merged' in subdir: continue
      if '.root' in subdir: continue

      print '\n-> Processing: {}'.format(subdir[subdir.rfind('/')+1:len(subdir)])

      # get files
      nanoFilesPlain = [f for f in glob.glob(subdir+nanoName)]

      nanoFiles = map(lambda x: 'root://t3dcachedb.psi.ch:1094/'+x, nanoFilesPlain)

      # create the outputdir that will contain the merged file
      if not path.exists(subdir+outputdir):
        os.system('mkdir {}'.format(subdir+outputdir))

      outputname = '{}/{}'.format(subdir+outputdir, mergedName) if not self.dobatch else 'merge.root'
      command = 'python haddnano.py {}'.format(outputname)

      if len(nanoFiles) == 0: 
        print 'no files of interest in this chunk'

      else:
        print "\n-> Checking the files"
        for iFile, fileName in enumerate(nanoFiles):
          if iFile%100 == 0:              print '     --> checked {}% of the files'.format(round(float(iFile)/len(nanoFiles)*100, 1))
          elif iFile == len(nanoFiles)-1: print '     --> checked 100% of the files'

          if not NanoTools.checkLocalFile(self, fileName, cond, branch_check=True, branchname='nMuon'): continue
          command = command + ' {}'.format(fileName)

        print '\n-> Start of the merge'
        os.system(command)

        if self.dobatch:
          command_xrdcp = 'xrdcp -f merge.root root://t3dcachedb.psi.ch:1094/{}/{}'.format(subdir+outputdir, mergedName)
          os.system(command_xrdcp)

        print '{}/{} created \n'.format(subdir+outputdir, mergedName)


  def doChunkMerging(self, nanoName, mergedName, locationSE, cond=True):
    print '\n---> Merging the different chunks'

    nanoFilesPlain = [f for f in glob.glob(locationSE+'/Chunk*/'+nanoName)]
    nanoFiles = map(lambda x: 'root://t3dcachedb.psi.ch:1094/'+x, nanoFilesPlain)

    # create the outputdir that will contain the merged file
    if not path.exists('{}/merged'.format(locationSE)):
      os.system('mkdir {}/merged'.format(locationSE))

    filesValid = []
    print "\n-> Checking the files"
    for fileName in nanoFiles:
      if cond and not NanoTools.checkLocalFile(self, fileName, cond, branch_check=True, branchname='nMuon'): continue
      elif not cond and not NanoTools.checkLocalFile(self, fileName, cond): continue
      filesValid.append(fileName)

    print '\n-> Start of the merge'
    outputname = '{}/merged/{}'.format(locationSE, mergedName) if not self.dobatch else 'fullmerge.root'
    command = 'python haddnano.py {}'.format(outputname)
    for iFile, fileName in enumerate(filesValid):
       command = command + ' {}'.format(fileName)

    os.system(command)

    if self.dobatch:
      command_xrdcp = 'xrdcp -f fullmerge.root root://t3dcachedb.psi.ch:1094/{}/merged/{}'.format(locationSE, mergedName)
      os.system(command_xrdcp)

    print '{}/merged/{} created \n'.format(locationSE, mergedName)

    # clean
    print 'cleaning'
    for f in glob.glob(locationSE+'/Chunk*/merged/'):
      command_clean_file = 'rm -rf root://t3dcachedb.psi.ch:1094/{}/{}'.format(f, mergedName)
      os.system(command_clean_file)


  def doExtMerging(self, mergedName, location):
    print '\n ---> Merging the extension files'

    no_ext_file = location[:len(location)-4] + '/merged/' + mergedName
    ext_file = location + '/merged/' + mergedName
    extmerged_file = location + '/merged/' + mergedName[:len(mergedName)-5] + '_extmerged.root'

    # copying the files in the workdir
    command_cp_no_ext_file = 'xrdcp {} ./no_ext_file.root'.format(no_ext_file)
    command_cp_ext_file = 'xrdcp {} ./ext_file.root'.format(ext_file)
    os.system(command_cp_no_ext_file)
    os.system(command_cp_ext_file)

    # proceed to the merging
    command = 'hadd -f extmerged_file.root ext_file.root no_ext_file.root'
    os.system(command)

    # copy the merged file
    command_cp_extmerged_file = 'xrdcp -f extmerged_file.root root://t3dcachedb.psi.ch:1094/{}'.format(extmerged_file)
    os.system(command_cp_extmerged_file)

    # erasing the files
    command_rm_no_ext_file = 'rm no_ext_file.root' 
    command_rm_ext_file = 'rm ext_file.root' 
    command_rm_extmerged_file = 'rm extmerged_file.root' 
    os.system(command_rm_no_ext_file)
    os.system(command_rm_ext_file)
    os.system(command_rm_extmerged_file)
    
    print '{} created \n'.format(extmerged_file)


  def runMergingModule(self, location):
    if self.donano:
      nanoName   = '/bparknano_nj*.root' if self.tagnano == None else '/bparknano_{}_nj*.root'.format(self.tagnano)
      mergedName = 'bparknano.root' if self.tagnano == None else 'bparknano_{}.root'.format(self.tagnano)
      outputdir  = '/merged'
          
      nanoName_tot   = 'merged/bparknano.root' if self.tagnano == None else 'merged/bparknano_{}.root'.format(self.tagnano)
      mergedName_tot = 'bparknano.root' if self.tagnano == None else 'bparknano_{}.root'.format(self.tagnano)

      # per chunk
      self.doMerging(nanoName, mergedName, location, outputdir, True)

      # inclusive
      self.doChunkMerging(nanoName_tot, mergedName_tot, location)

    if self.doflat:
      tag = NanoTools.getTag(self, self.tagnano, self.tagflat)

      nanoName_flat   = 'flat/flat_bparknano.root' if self.tagnano == None and self.tagflat == None else 'flat/flat_bparknano_{}.root'.format(tag)
      mergedName_flat = 'flat_bparknano.root' if self.tagnano == None and self.tagflat == None else 'flat_bparknano_{}.root'.format(tag)

      self.doChunkMerging(nanoName_flat, mergedName_flat, location, False)

      if self.mccentral and '_ext' in location:
        self.doExtMerging(mergedName_flat, location)


  def process(self):
    print '---------------------------------'
    print '           Nano Merger           '
    print '---------------------------------'

    user      = os.environ["USER"]

    if self.mcprivate:
      locationSE = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/'.format(user, self.prodlabel)
      
      pointdirs = NanoTools.getPointDirs(self, locationSE)

      for pointdir in pointdirs:
        if 'merged' in pointdir: continue

        point = pointdir[pointdir.rfind('/')+1:len(pointdir)]
        print '\n --- Mass point: {} --- '.format(point)

        self.runMergingModule(pointdir+'/nanoFiles/')

    elif self.data or self.mccentral:
      dataset_label = NanoTools.getDataLabel(self, self.dataset) if self.data else NanoTools.getMCLabel(self, self.dataset)
      locationSE = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/{}/{}'.format(user, 'data' if self.data else 'mc_central', self.prodlabel, dataset_label)

      self.runMergingModule(locationSE)
  


if __name__ == "__main__":

  opt = getOptions()

  checkParser(opt)

  NanoMerger(opt).process()

  print 'Done'


