# necessitates python3: source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh

import os
import ROOT
import os.path
from os import path
from myBranches import the_signal_branches, the_control_branches, the_common_branches, the_extra_branches


class NanoDumper(object):
  def __init__(self, nano_file='bparknano.root', flat_file='flat_bparknano.root', ncandtrees=3):
    self.nano_file = nano_file
    self.flat_file = flat_file
    self.ncandtrees  = ncandtrees
  

  def createTree(self, channel, mode, cand_idx=0):
    # number of candidates
    ncand = 'nBToMuMuPi' if channel == 'signal' else 'nBToKMuMu'

    # creating dataframe
    cand_dataframe = ROOT.RDataFrame('Events', self.nano_file).Filter("{ncand} > {idx}".format(ncand=ncand, idx=cand_idx))

    # fetching the branches
    branches_list = []

    the_branches = the_common_branches+the_signal_branches if channel == 'signal' else the_common_branches+the_control_branches

    branch_name = 'BToMuMuPi' if channel == 'signal' else 'BToKMuMu'
    for flat_name, nano_name in the_branches:
      full_nano_name = '{branch}_{nano}'.format(branch=branch_name, nano=nano_name) 
      if mode == 'percand':
        full_nano_name += '[{idx}]'.format(idx=cand_idx)
      branches_list.append((flat_name, full_nano_name))
     
    ## add extra branches
    if mode == 'allcands': 
      for flat_name, nano_name in the_extra_branches:
        full_nano_name = nano_name 
        branches_list.append((flat_name, full_nano_name))

      ### add number of candidates information
      ncand = 'n'+branch_name
      branches_list.append(
            ('ncand', ncand)
      )

    # add branches to dataframe
    branches = ROOT.vector('string')()
    for branchName, branchDef in branches_list:
      cand_dataframe = cand_dataframe.Define(branchName, branchDef) 
      branches.push_back(branchName)

    # snapshot  
    ## disable overwriting
    snapopts = ROOT.RDF.RSnapshotOptions()
    snapopts.fMode = 'UPDATE'

    tree_name = 'allcands_tree' if mode == 'allcands' else 'cand{idx}_tree'.format(idx=cand_idx)
    print('\n --> Creating tree {ch}_channel/{tree}'.format(ch=channel, tree=tree_name))
    cand_dataframe.Snapshot('{ch}_channel/{tree}'.format(ch=channel, tree=tree_name), self.flat_file, branches, snapopts);


  def process(self):
     print('-----------------------')
     print('     Nano Dumper       ')
     print('-----------------------')

     if path.exists(self.flat_file):
        command_rm = 'rm {}'.format(self.flat_file)
        os.system(command_rm)

     self.createTree('signal', 'allcands')
     for idx in range(0, self.ncandtrees):
        self.createTree('signal', 'percand', idx)

     self.createTree('control', 'allcands')
     for idx in range(0, self.ncandtrees):
        self.createTree('control', 'percand', idx)


if __name__ == '__main__':

  nano_file = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/ParkingBPH1_Run2018B_F1/Chunk0_n750/bparknano_nj95.root'
  NanoDumper(nano_file=nano_file, ncandtrees=5).process()

