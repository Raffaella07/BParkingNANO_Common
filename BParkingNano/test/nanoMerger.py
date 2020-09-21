import sys 
import os
import glob
import ROOT


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to merge the nanoAOD files resulting from a multijob production', add_help=True)
  parser.add_argument('--pl', type=str, dest='pl', help='label of the sample file', default='V01_n9000000_njt300')
  parser.add_argument('--tag', type=str, dest='tag', help='[optional] tag to be added on the outputfile name', default=None)
  return parser.parse_args()


if __name__ == "__main__":

  opt = getOptions()

  user      = os.environ["USER"]

  locationSE = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/'.format(user, opt.pl)
  nanoName = '/nanoFiles/bparknano_nj*.root' if opt.tag == None else '/nanoFiles/bparknano_{}_nj*.root'.format(opt.tag)

  print '\n-> Getting the different mass points'
  pointsdir = [f for f in glob.glob(locationSE+'/*')]

  for point in pointsdir:
    print '\n-> Processing mass/ctau point: {}'.format(point[point.rfind('/')+1:len(point)])

    # get files
    nanoFiles = [f for f in glob.glob(point+nanoName)]

    # create the outputdir that will contain the merged file
    outputdir = '{}/nanoFiles/merged'.format(point)
    os.system('mkdir {}'.format(outputdir))

    filesValid = []
    print "\n-> Checking the files"
    for fileName in nanoFiles:
       rootFile = ROOT.TNetXNGFile.Open(fileName, 'r')
       if rootFile and rootFile.GetListOfKeys().Contains('Events'):
          filesValid.append(fileName)

    print '\n-> Start of the merge'
    mergedName = 'bparknano.root' if opt.tag == None else 'bparknano_{}.root'.format(opt.tag)
    command = 'python haddnano.py {}/{}'.format(outputdir, mergedName)
    for fileName in filesValid:
       command = command + ' {}'.format(fileName)

    os.system(command)

  print 'Done'









