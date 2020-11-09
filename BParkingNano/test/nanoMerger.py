import sys 
import os
import glob
import ROOT


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to merge the nanoAOD files resulting from a multijob production', add_help=True)
  parser.add_argument('--pl', type=str, dest='pl', help='label of the sample file', default='V01_n9000000_njt300')
  parser.add_argument('--tag', type=str, dest='tag', help='[optional] tag to be added on the outputfile name', default=None)
  parser.add_argument('--doflat', dest='doflat', help='merge as well the files coming out of the truth-mathing', action='store_true', default=False)
  return parser.parse_args()


def doMerging(nanoName, mergedName, outputdir, cond):

  print '\n-> Getting the different mass points'
  pointsdir = [f for f in glob.glob(locationSE+'/*')]

  for point in pointsdir:
    print '\n-> Processing mass/ctau point: {}'.format(point[point.rfind('/')+1:len(point)])

    # get files
    nanoFiles = [f for f in glob.glob(point+nanoName)]

    # create the outputdir that will contain the merged file
    outputdir = point + outputdir
    os.system('mkdir {}'.format(outputdir))

    filesValid = []
    print "\n-> Checking the files"
    for fileName in nanoFiles:
      rootFile = ROOT.TNetXNGFile.Open(fileName, 'r')
      if rootFile: 
        if cond:
          if rootFile.GetListOfKeys().Contains('Events'):
            filesValid.append(fileName)
        else:
          filesValid.append(fileName)

    print '\n-> Start of the merge'
    command = 'python haddnano.py {}/{}'.format(outputdir, mergedName)
    for iFile, fileName in enumerate(filesValid):
       #if iFile > 2: continue
       command = command + ' {}'.format(fileName)

    os.system(command)

  print '{}/{} created'.format(outputdir, mergedName)


if __name__ == "__main__":

  opt = getOptions()

  user      = os.environ["USER"]

  locationSE = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/'.format(user, opt.pl)
 
  nanoName   = '/nanoFiles/bparknano_nj*.root' if opt.tag == None else '/nanoFiles/bparknano_{}_nj*.root'.format(opt.tag)
  mergedName = 'bparknano.root' if opt.tag == None else 'bparknano_{}.root'.format(opt.tag)
  outputdir  = '/nanoFiles/merged'
  
  nanoName_flat   = '/nanoFiles/flat/flat_bparknano_nj*.root' if opt.tag == None else '/nanoFiles/flat/flat_bparknano_{}_nj*.root'.format(opt.tag)
  mergedName_flat = 'flat_bparknano.root' if opt.tag == None else 'flat_bparknano_{}.root'.format(opt.tag)
  outputdir_flat  = '/nanoFiles/flat/merged'

  doMerging(nanoName, mergedName, outputdir, True)
  if opt.doflat:
    doMerging(nanoName_flat, mergedName_flat, outputdir_flat, False)

  print 'Done'




