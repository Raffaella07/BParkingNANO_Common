import uproot
import pandas as pd
import numpy as np
import awkward
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('f_old', help='file path')
parser.add_argument('f_new', help='file path')
parser.add_argument('--legacy', action='store_true', help='compare against legacy version')
args = parser.parse_args()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

if not os.path.isdir('validation'):
  os.makedirs('validation')

legacy_mapping = { #mapping between George's Ntuples and Nano
  'nmuon' : 'nMuon',
  'muon_eta' : 'Muon_eta',
  # TODO: Add variables as they are validated and produced
}

def load_nano(infile, legacy = False):
  uf = uproot.open(infile)
  tt = uf['demo/mytree'] if legacy else uf['Events']
  arrs = tt.arrays(tt.keys())
  if legacy:
    return {j : arrs[i] for i, j in legacy_mapping.iteritems()}
  return arrs

def byval_validation(v1, v2):
  try:
    return ((np.abs(v1 - v2) / abs(v1)) < 0.001).all()
  except ValueError:
    return False

def stat_validation(v1, v2, name = '', nbins = 20):
  if v1.shape[0] == 0 and v2.shape[0] == 0:
    return True
  elif v1.shape[0] == 0 or v2.shape[0] == 0:
    return False
  M = max(v1.max(), v2.max())*1.2
  m = min(v1.min(), v2.min())*0.9
  plt.clf()
  h1, _, _ = plt.hist(v1, range = (m,M), bins = nbins, label = 'old', histtype = 'step')
  h2, _, _ = plt.hist(v2, range = (m,M), bins = nbins, label = 'new', histtype = 'step')
  plt.legend(loc='best')
  plt.savefig('validation/%s.png' % name)
  return (h1 == h2).all()

old = load_nano(args.f_old, args.legacy)
new = load_nano(args.f_new)

old_k = set(old.keys())
new_k = set(new.keys())
intersection = old_k.intersection(new_k)

for branch in intersection:
  v_old = old[branch]
  v_new = new[branch]
  if hasattr(v_old, 'flatten'):
    v_old = v_old.flatten()
    v_new = v_new.flatten()
  stat_valid = stat_validation(v_old, v_new, branch)
  val_valid  = byval_validation(v_old, v_new)
  if val_valid and stat_valid:
    print '\033[1;32m', branch, '--> OK!\033[0m'
  elif stat_valid:
    print '\033[1;35m', branch, '--> FAILS BY VALUE CHECK ONLY!\033[0m'
  else:
    print '\033[1;31m', branch, '--> FAILS ALL CHECKS!\033[0m'
