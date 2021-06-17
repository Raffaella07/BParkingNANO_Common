from __future__ import print_function
import ROOT
import argparse
import numpy as np
from time import time
from datetime import datetime, timedelta
from array import array
from glob import glob
from collections import OrderedDict
from scipy.constants import c as speed_of_light
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi, bestMatch, matchObjectCollection3
# https://pypi.org/project/particle/
#from particle import Particle


# load custom library to ROOT. This contains the kinematic vertex fitter class
ROOT.gSystem.Load('libPhysicsToolsBParkingNano')
from ROOT import HNLKalmanVertexFitter as VertexFitter
vf = VertexFitter()
tofit = ROOT.std.vector('reco::Track')()
tofit.clear()


class MuTk(object):
  def __init__(self,mu,tk):
    self.mu = mu
    self.tk = tk
    self.mu_p4 = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<float>')(mu.pt(),mu.eta(),mu.phi(),0.105658) # assumption that this is a muon
    self.tk_p4 = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<float>')(tk.pt(),tk.eta(),tk.phi(),0.139571) # assumption that this is a pion... 
    self.p4 =   self.mu_p4 + self.tk_p4
    self.mass = self.p4.mass()
    self.isSS = int(self.mu.charge()==self.tk.charge())


    
handles = OrderedDict()
handles['genp'    ] = ('genParticles'   , Handle('std::vector<reco::GenParticle>'))
handles['genInfo' ] = ('generator'      , Handle('GenEventInfoProduct'           ))
handles['beamspot'] = ('offlineBeamSpot', Handle('reco::BeamSpot '               ))
handles['muons'        ] = ('muons'                    , Handle('std::vector<reco::Muon>'       ))
handles['muons_dsa'    ] = ('displacedStandAloneMuons' , Handle('std::vector<reco::Track>'       ))
handles['muons_dg'     ] = ('displacedGlobalMuons'     , Handle('std::vector<reco::Track>'       ))
handles['tracks'       ] = ('generalTracks'            , Handle('std::vector<reco::Track>'      ))   
handles['tracks_disp'  ] = ('displacedTracks'          , Handle('std::vector<reco::Track>'      ))   

events = Events( ['root://t3dcachedb.psi.ch:1094/' + i for i in glob('/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V20_test/mass3.0_ctau184.256851021/step3_nj*.root')])

doAddDisplaced = True
doDispGMuons   = True
doDispSTMuons = True
doDisplTracks  = True

suffix = ''
if doDispGMuons:
  suffix = '_dg'
  if doDispSTMuons:
    suffix = '_dsadg'
    if doDisplTracks:
      suffix = '_displAll'

fout = ROOT.TFile('bhnl_recoeff{}.root'.format(suffix), 'recreate')

branches = [
    'run',
    'lumi',
    'event',
    'hnl_mass'    ,
    'hnl_pt'      ,
    'hnl_eta'     ,
    #'hnl_y'       ,
    'hnl_phi'     ,
    'hnl_q'       ,
    'hnl_vx'      ,
    'hnl_vy'      ,
    'hnl_vz'      ,

    'mu_pt'      ,
    'mu_eta'     ,
    'mu_y'       ,
    'mu_phi'     ,
    'mu_q'       ,
    'mu_vx'      ,
    'mu_vy'      ,
    'mu_vz'      ,

    'tk_pt'      ,
    'tk_eta'     ,
    'tk_y'       ,
    'tk_phi'     ,
    'tk_q'       ,
    'tk_vx'      ,
    'tk_vy'      ,
    'tk_vz'      ,

    'hnl_gen_mass',
    'hnl_gen_pt'  ,
    'hnl_gen_eta' ,
    'hnl_gen_y'   ,
    'hnl_gen_phi' ,
    'hnl_gen_q'   ,
    'hnl_gen_vx'  ,
    'hnl_gen_vy'  ,
    'hnl_gen_vz'  ,

    'mu_gen_pt'  ,
    'mu_gen_eta' ,
    'mu_gen_y'   ,
    'mu_gen_phi' ,
    'mu_gen_q'   ,
    'mu_gen_vx'  ,
    'mu_gen_vy'  ,
    'mu_gen_vz'  ,

    'pi_gen_pt'  ,
    'pi_gen_eta' ,
    'pi_gen_y'   ,
    'pi_gen_phi' ,
    'pi_gen_q'   ,
    'pi_gen_vx'  ,
    'pi_gen_vy'  ,
    'pi_gen_vz'  ,

    'lxy_gen'     ,
    'lxy'         ,
]

ntuple = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))


#ntuple = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
#tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))
bins = [0.,1.,5.,10.,15.,30.,50.,100.]
hnum = ROOT.TH1F('hnum', 'hnum', len(bins)-1, array('d',bins))
hden = ROOT.TH1F('hden', 'hden', len(bins)-1, array('d',bins))

bins_pT = [1.,5.,10.,20.,100.]
hnum_pT = ROOT.TH1F('hnum_pT', 'hnum_pT', len(bins_pT)-1, array('d',bins_pT))
hden_pT = ROOT.TH1F('hden_pT', 'hden_pT', len(bins_pT)-1, array('d',bins_pT))

hnum_2D = ROOT.TH2F('hnum_2D', 'hnum_2D', len(bins)-1, array('d',bins),len(bins_pT)-1, array('d',bins_pT))
hden_2D = ROOT.TH2F('hden_2D', 'hden_2D', len(bins)-1, array('d',bins),len(bins_pT)-1, array('d',bins_pT))

start = time()
maxevents = -1
maxevents = maxevents if maxevents>=0 else events.size() # total number of events in the files

for i, event in enumerate(events):

    if (i+1)>maxevents:
        break
        
    if i%100==0:
    #if i%1==0:
        percentage = float(i)/maxevents*100.
        speed = float(i)/(time()-start)
        eta = datetime.now() + timedelta(seconds=(maxevents-i) / max(0.1, speed))
        print('\t===> processing %d / %d event \t completed %.1f%s \t %.1f ev/s \t ETA %s s' %(i, maxevents, percentage, '%', speed, eta.strftime('%Y-%m-%d %H:%M:%S')))

    # initialise ntuple
    for k, v in tofill.items(): 
      tofill[k] = np.nan

    # access the handles
    for k, v in handles.iteritems():
        event.getByLabel(v[0], v[1])
        setattr(event, k, v[1].product())

    
    print('PAT MUONS:')
    stamp(event.muons)
    print('DSA MUONS:')
    stamp(event.muons_dsa)
    print('DSA MUONS:')
    stamp(event.muons_dg)

    #event.muons = event.muons_reco
    all_reco_muons = [imu for imu in event.muons]
    all_reco_tracks = [itk for itk in event.tracks]

    if doDispGMuons:
      all_reco_muons +=   [imu for imu in event.muons_dg]
      if doDispSTMuons:
        all_reco_muons += [imu for imu in event.muons_dsa] 
    if doDisplTracks:
      all_reco_tracks +=  [itk for itk in event.tracks_disp] 

           
    hnls_gen = [ip for ip in event.genp if abs(ip.pdgId())==9900015]
    hnl_gen = hnls_gen[0]

    hnl_gen_daus = sorted([hnl_gen.daughter(0), hnl_gen.daughter(1)], key = lambda x : x.pt(), reverse=True)
    mu_gen  = hnl_gen_daus[0] if abs(hnl_gen_daus[0].pdgId()) == 13  else hnl_gen_daus[1]
    pi_gen =  hnl_gen_daus[0] if abs(hnl_gen_daus[0].pdgId()) == 211 else hnl_gen_daus[1]

    # loose preselection on the gen
    if abs(pi_gen.pt()) < 1 or abs(mu_gen.pt()) < 1: continue
    if abs(pi_gen.eta())> 2.5 or abs(mu_gen.eta()) > 2.5: continue

    # fill denominator
    gen_lxy = np.sqrt( (pi_gen.vx() - hnl_gen.vx() )**2 + (pi_gen.vy() - hnl_gen.vy() )**2 ) 
    hden.Fill(gen_lxy)
    hden_pT.Fill(hnl_gen.pt())
    hden_2D.Fill(gen_lxy,hnl_gen.pt())

    # gen matching 
    # first a loose preselection, only consider objects with Dr in the vicinity
    reco_muons = [imu for imu in all_reco_muons  if deltaR(imu,mu_gen) < 0.3 ]
    reco_tks =   [itk for itk in all_reco_tracks if deltaR(itk,pi_gen) < 0.3 ]
    #reco_mutks = [(imu,itk) for imu in reco_muons for itk in tks]
    reco_mutks = [MuTk(imu,itk) for imu in reco_muons for itk in reco_tks if imu.pt() != itk.pt()]
#    print(len(reco_mutks))

    if len(reco_mutks)>0:

      # then sort the couples by mass 
      reco_mutks = sorted(reco_mutks, key = lambda x: (x.isSS,abs(x.mass-3.)), reverse=False)
      to_print = ['{}_{}'.format(x.isSS,abs(x.mass-3.)) for x in reco_mutks]
      #best_mutk = reco_mutks[0]
     
      for imutk in reco_mutks:
      
        if abs(imutk.mass-3) < 0.5:  # and this is loose..
  
          tofit.clear()
          mu = imutk.mu
          tk = imutk.tk
          if hasattr(mu,'muonBestTrack'):
            tofit.push_back(mu.muonBestTrack().get())
          else:
            tofit.push_back(mu)
          tofit.push_back(tk)
          svkalman = vf.Fit(tofit)
          sv = svkalman
          if not sv.isValid():
              continue
  
          reco_lxy = np.sqrt( (event.beamspot.x0() - sv.position().x())**2 + (event.beamspot.y0() - sv.position().y())**2 )
          #print(reco_lxy,gen_lxy)
  
          if(abs(reco_lxy/gen_lxy-1)<0.5 or abs(gen_lxy/reco_lxy-1)<0.5):
           
            # fill histogram
            hnum.Fill(gen_lxy)
            hnum_pT.Fill(hnl_gen.pt())
            hnum_2D.Fill(gen_lxy,hnl_gen.pt())
            print(reco_lxy,gen_lxy,imutk.mass)
            # Fill ntuple
            tofill['run'     ] = event.eventAuxiliary().run()
            tofill['lumi'    ] = event.eventAuxiliary().luminosityBlock()
            tofill['event'   ] = event.eventAuxiliary().event()
     
            tofill['hnl_mass' ] = imutk.mass
            tofill['hnl_pt'   ] = imutk.p4.pt()
            tofill['hnl_eta'  ] = imutk.p4.eta()
            tofill['hnl_phi'  ] = imutk.p4.phi()

            tofill['hnl_q'   ] = mu.charge() + tk.charge()
            tofill['hnl_vx'  ] = sv.position().x()
            tofill['hnl_vy'  ] = sv.position().y()
            tofill['hnl_vz'  ] = sv.position().z()
        
            tofill['mu_pt'  ] = mu.pt()
            tofill['mu_eta' ] = mu.eta()
            #tofill['mu_y'   ] = mu1.y()
            tofill['mu_phi' ] = mu.phi()
            tofill['mu_q'   ] = mu.charge()
            tofill['mu_vx'  ] = mu.vx()
            tofill['mu_vy'  ] = mu.vy()
            tofill['mu_vz'  ] = mu.vz()
        
            tofill['tk_pt'  ] = tk.pt()
            tofill['tk_eta' ] = tk.eta()
            #tofill['tk_y'   ] = pi2.y()
            tofill['tk_phi' ] = tk.phi()
            tofill['tk_q'   ] = tk.charge()
            tofill['tk_vx'  ] = tk.vx()
            tofill['tk_vy'  ] = tk.vy()
            tofill['tk_vz'  ] = tk.vz()
        
            tofill['hnl_gen_mass'] = hnl_gen.mass()
            tofill['hnl_gen_pt'  ] = hnl_gen.pt()
            tofill['hnl_gen_eta' ] = hnl_gen.eta()
            tofill['hnl_gen_y'   ] = hnl_gen.y()
            tofill['hnl_gen_phi' ] = hnl_gen.phi()
            tofill['hnl_gen_q'   ] = mu_gen.charge() + pi_gen.charge()
            tofill['hnl_gen_vx'  ] = hnl_gen.vx()
            tofill['hnl_gen_vy'  ] = hnl_gen.vy()
            tofill['hnl_gen_vz'  ] = hnl_gen.vz()
        
            tofill['mu_gen_pt'  ] = mu_gen.pt()
            tofill['mu_gen_eta' ] = mu_gen.eta()
            tofill['mu_gen_y'   ] = mu_gen.y()
            tofill['mu_gen_phi' ] = mu_gen.phi()
            tofill['mu_gen_q'   ] = mu_gen.charge()
            tofill['mu_gen_vx'  ] = mu_gen.vx()
            tofill['mu_gen_vy'  ] = mu_gen.vy()
            tofill['mu_gen_vz'  ] = mu_gen.vz()
        
            tofill['pi_gen_pt'  ] = pi_gen.pt()
            tofill['pi_gen_eta' ] = pi_gen.eta()
            tofill['pi_gen_y'   ] = pi_gen.y()
            tofill['pi_gen_phi' ] = pi_gen.phi()
            tofill['pi_gen_q'   ] = pi_gen.charge()
            tofill['pi_gen_vx'  ] = pi_gen.vx()
            tofill['pi_gen_vy'  ] = pi_gen.vy()
            tofill['pi_gen_vz'  ] = pi_gen.vz()
        
            tofill['lxy_gen'     ] = np.sqrt( (pi_gen.vx()        - hnl_gen.vx()     )**2 + (pi_gen.vy()        - hnl_gen.vy()     )**2 ) 
            tofill['lxy'         ] = np.sqrt( (event.beamspot.x0() - sv.position().x())**2 + (event.beamspot.y0() - sv.position().y())**2 )
      
            ntuple.Fill(array('f', tofill.values()))

            break # as soon as you find a decent one, go to next event




fout.cd()
ntuple.Write()
fout.Write()
fout.Close()

