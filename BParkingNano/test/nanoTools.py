import glob
import ROOT


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to launch the nanoAOD tool on top of privately produced miniAOD files', add_help=True)
  parser.add_argument('--outdir'      , type=str, dest='outdir'      , help='dir where nanofiles are stored'                      , default=None)
  parser.add_argument('--tag'         , type=str, dest='tag'         , help='[optional] tag of the outputfile name'               , default='0')
  parser.add_argument('--ismc'        , type=str, dest='ismc'        , help='sample is mc'                                        , default='0')
  parser.add_argument('--writestarter'          , dest='writestarter', help='write starter to process dumper', action='store_true', default=False)
  return parser.parse_args()

class NanoTools(object):
  def __init__(self, opt):
    self.ismc    = vars(opt)['ismc']
    self.tag     = vars(opt)['tag']
    self.outdir  = vars(opt)['outdir']


  def getPointDirs(self, location):
    return [f for f in glob.glob(location+'/*')]


  def writeDumperStarter(self):
    nanoname = 'bparknano_nj' if self.tag == '0' else 'bparknano_{}_nj'.format(self.tag) 
    nanofiles = [f for f in glob.glob('{}/{}*'.format(self.outdir, nanoname))]

    event_chain = []
    event_chain.append('TChain* c = new TChain("Events");')
    for nanofile in nanofiles:
      rootFile = ROOT.TNetXNGFile.Open(nanofile, 'r')
      if rootFile and rootFile.GetListOfKeys().Contains('Events'):
        event_chain.append('  c->Add("{}");'.format(nanofile))
    event_chain.append('  c->Process("NanoDumper.C+", outFileName);')
    event_chain = '\n'.join(event_chain)

    run_chain = []
    run_chain.append('TChain* c_run = new TChain("Runs");')
    for nanofile in nanofiles:
      rootFile = ROOT.TNetXNGFile.Open(nanofile, 'r')
      if rootFile and rootFile.GetListOfKeys().Contains('Events'):
        run_chain.append('  c_run->Add("{}");'.format(nanofile))
    run_chain.append('  c_run->Process("NanoRunDumper.C+", outFileName);')
    run_chain = '\n'.join(run_chain)

    content = [
      '#include "TChain.h"',
      '#include <iostream>',
      '#include "TProof.h"\n',
      'void starter(){',
      '  TString outFileName = "flat_bparknano.root";',
      '  {addevt}'.format(addevt = event_chain),
      '  {addrun}'.format(addrun = '' if self.ismc == '0' else run_chain),
      '}',
    ]
    content = '\n'.join(content)
          
    dumper_starter = open('starter.C', 'w+')
    dumper_starter.write(content)
    dumper_starter.close()

    
if __name__ == "__main__":

  opt = getOptions()

  if vars(opt)['writestarter']:
    NanoTools(opt).writeDumperStarter()


