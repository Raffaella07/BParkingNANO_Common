from root_numpy import root2array, tree2array, fill_hist,array2tree
from root_numpy import testdata
from rootpy.tree import Tree, TreeModel
from rootpy.tree import IntCol, FloatCol, FloatArrayCol, CharCol, CharArrayCol
from rootpy.io import root_open
from ctypes import *
import ROOT
import numpy as np
import glob
import math
from ROOT import TH1D
from operator import itemgetter
import multiprocessing
import sys
sys.path.append("/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/python/2.7.14-omkpbe4/lib/python2.7/site-packages/")
from sklearn.externals.joblib import Parallel, delayed
from tqdm import tqdm


class Event(TreeModel):
    evtIdx = IntCol()
    nCand= IntCol()
    Type = IntCol() # idx for channel category : pi-mu :0 , pi-PF :1, pi-lowpT: 2, track track 3 
    LxyBin = IntCol()	# idx for displacement bin : lxy< 3 : 0, lxy <10 : 1 , lxy < 20:2, lxy >20:3
    B_mass= FloatCol()
    B_pt= FloatCol()
    B_eta= FloatCol()
    B_phi= FloatCol()
    TrgMu_pt= FloatCol()
    TrgMu_eta= FloatCol()
    TrgMu_phi= FloatCol()
    hnl_mass= FloatCol()
    hnl_pt= FloatCol()
    hnl_lxy= FloatCol()
    hnl_cos2D= FloatCol()
    hnl_lxy_sig= FloatCol()
    hnl_drLepPi = FloatCol()
    dr_trgmu_hnl = FloatCol()
    hnl_vtxProb = FloatCol()
    hnl_vtxChi2 =  FloatCol()
    hnl_vx =  FloatCol()
    hnl_ex =  FloatCol()
    hnl_vy =  FloatCol()
    hnl_ey =  FloatCol()
    hnl_vz =  FloatCol()
    hnl_ez =  FloatCol()
    hnl_l_pt =  FloatCol()
    hnl_l_eta =  FloatCol()
    hnl_l_phi =  FloatCol()
    hnl_l_dz =  FloatCol()
    hnl_l_dzS =  FloatCol()
    hnl_l_dxy =  FloatCol()
    hnl_l_dxyS =  FloatCol()
    hnl_l_DCAS =  FloatCol()
    hnl_l_mvaId =  FloatCol()
    hnl_pi_pt =  FloatCol()
    hnl_pi_eta =  FloatCol()
    hnl_pi_phi =  FloatCol()
    hnl_pi_dz =  FloatCol()
    hnl_pi_dzS =  FloatCol()
    hnl_pi_dxy =  FloatCol()
    hnl_pi_dxyS =  FloatCol()
    hnl_pi_DCAS =  FloatCol()

#def createHistos(full, chString):
#  listHist = []
#  #print len(l[0][0])
#  for j in range(0,len(lxy_lables)):
#	temp = []
#   	for i in range(0,len(in_branches)):
#	
##	temp.resize(4)
##	print i
#		temp.append(ROOT.TH1D(chString+str(i)+str(j),chString+in_branches[i]+lxy_lables[j], Nbin[i], binMin[i], binMax[i]))
##	print len(temp)
#	listHist.append(temp)
#  #print len(l[0])	
##  print len(listHist)
#  return listHist   

def sortBranches(l,full, chString):
	J = 1.0*j
	if len(l[0]!=0):
 	   
         test = []
#  	 print "BEFORE SORT"
 #    	 print l
 	 test = l
 #	if j % len(rta_charge[0]) == 1000:
 	#	print('{:.3f}'.format(J/len(rta_charge[0])))
 #	print k
 #	print test[1]
	 if not full:
        	 test = zip(*test) 
      	         test = sorted(test,key= itemgetter(1) ,reverse=True)
                 test = zip(*test)
         l = test



n_branches = 8
binMin = [0,0,0,0,0,0,0] 
binMax = [7,20,100,1,10,50,50] 
Nbin = [100,100,100,100,100,100,100]
sig = ["BToMuMuPi","BToMuPFPi","BToMuLowPtPi","BToMutt"]
lxyMu_bins = ["BToMuMuPi_sv_lxy<3","(BToMuMuPi_sv_lxy>3 && BToMuMuPi_sv_lxy<10)","(BToMuMuPi_sv_lxy>10 && BToMuMuPi_sv_lxy<20)","BToMuMuPi_sv_lxy>20"]
lxyEle_bins = ["BToMuEPi_sv_lxy<3","BToMuEPi_sv_lxy>3 && BToMuEPi_sv_lxy<10","BToMuEPi_sv_lxy>10 && BToMuEPi_sv_lxy<20","BToMuEPi_sv_lxy>20"]
lxytt_bins = ["BToMuEPiHD_sv_lxy<3","BToMuEPiHD_sv_lxy>3 && BToMuEPiHD_sv_lxy<10","BToMuEPiHD_sv_lxy>10 && BToMuEPiHD_sv_lxy<20","BToMuEPiHD_sv_lxy>20","ProbeTracks_dz[BToMuLPiHD_pi_idx]"]
lxy_lables = ["lxyUnder3","lxyUnder10","lxyUnder20","lxyOver20"]
leaf_lables = ["Hnl_mass","Hnl_pt","lxy","cos2D","Bmass","Bpt","lxy_sig"]
mu_branches = ["BToMuMuPi_mass","BToMuMuPi_pt","BToMuMuPi_eta","BToMuMuPi_phi",
               "TriggerMuon_pt[0]","TriggerMuon_eta[0]","TriggerMuon_phi[0]",
               "BToMuMuPi_hnl_mass","BToMuMuPi_hnl_pt","BToMuMuPi_sv_lxy","BToMuMuPi_hnl_cos2D","BToMuMuPi_sv_lxy_sig","BToMuMuPi_dr_mu_pi",
               "BToMuMuPi_dr_trgmu_hnl","BToMuMuPi_sv_prob","BToMuMuPi_sv_chi2","BToMuMuPi_sv_x","BToMuMuPi_sv_xe","BToMuMuPi_sv_y","BToMuMuPi_sv_ye","BToMuMuPi_sv_z","BToMuMuPi_sv_ze",
               "Muon_pt[BToMuMuPi_sel_mu_idx]","Muon_eta[BToMuMuPi_sel_mu_idx]","Muon_phi[BToMuMuPi_sel_mu_idx]","Muon_dz[BToMuMuPi_sel_mu_idx]","Muon_dzErr[BToMuMuPi_sel_mu_idx]",
               "Muon_dxy[BToMuMuPi_sel_mu_idx]","Muon_dxyErr[BToMuMuPi_sel_mu_idx]","Muon_matched_dpt[BToMuMuPi_sel_mu_idx]","Muon_matched_dr[BToMuMuPi_sel_mu_idx]",
               "ProbeTracks_pt[BToMuMuPi_pi_idx]","ProbeTracks_eta[BToMuMuPi_pi_idx]","ProbeTracks_phi[BToMuMuPi_pi_idx]","ProbeTracks_dz[BToMuMuPi_pi_idx]","ProbeTracks_dzS[BToMuMuPi_pi_idx]",
               "ProbeTracks_dxy[BToMuMuPi_pi_idx]","ProbeTracks_dxyS[BToMuMuPi_pi_idx]","ProbeTracks_DCASig[BToMuMuPi_pi_idx]"]

ele_branches = ["BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_eta","BToMuEPi_phi",
               "TriggerMuon_pt[0]","TriggerMuon_eta[0]","TriggerMuon_phi[0]",
               "BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_sv_lxy_sig","BToMuEPi_dr_mu_pi",
               "BToMuEPi_dr_trgmu_hnl","BToMuEPi_sv_prob","BToMuEPi_sv_chi2","BToMuEPi_sv_x","BToMuEPi_sv_xe","BToMuEPi_sv_y","BToMuEPi_sv_ye","BToMuEPi_sv_z","BToMuEPi_sv_ze",
               "Electron_pt[BToMuEPi_sel_e_idx]","Electron_eta[BToMuEPi_sel_e_idx]","Electron_phi[BToMuEPi_sel_e_idx]","Electron_dz[BToMuEPi_sel_e_idx]","Electron_dzErr[BToMuEPi_sel_e_idx]",
               "Electron_dxy[BToMuEPi_sel_e_idx]","Electron_dxyErr[BToMuEPi_sel_e_idx]","Electron_pfmvaId[BToMuEPi_sel_e_idx]","Electron_mvaId[BToMuEPi_sel_e_idx]",
               "ProbeTracks_pt[BToMuEPi_pi_idx]","ProbeTracks_eta[BToMuEPi_pi_idx]","ProbeTracks_phi[BToMuEPi_pi_idx]","ProbeTracks_dz[BToMuEPi_pi_idx]","ProbeTracks_dzS[BToMuEPi_pi_idx]",
               "ProbeTracks_dxy[BToMuEPi_pi_idx]","ProbeTracks_dxyS[BToMuEPi_pi_idx]","ProbeTracks_DCASig[BToMuEPi_pi_idx]"
		]

tt_branches = ["BToMuEPiHD_mass","BToMuEPiHD_pt","BToMuEPiHD_eta","BToMuEPiHD_phi",
               "TriggerMuon_pt[0]","TriggerMuon_eta[0]","TriggerMuon_phi[0]",
               "BToMuEPiHD_hnl_mass","BToMuEPiHD_hnl_pt","BToMuEPiHD_sv_lxy","BToMuEPiHD_hnl_cos2D","BToMuEPiHD_sv_lxy_sig","BToMuEPiHD_dr_mu_pi",
               "BToMuEPiHD_dr_trgmu_hnl","BToMuEPiHD_sv_prob","BToMuEPiHD_sv_chi2","BToMuEPiHD_sv_x","BToMuEPiHD_sv_xe","BToMuEPiHD_sv_y","BToMuEPiHD_sv_ye","BToMuEPiHD_sv_z","BToMuEPiHD_sv_ze",
               "ProbeTracks_pt[BToMuEPiHD_sel_e_idx]","ProbeTracks_eta[BToMuEPiHD_sel_e_idx]","ProbeTracks_phi[BToMuEPiHD_sel_e_idx]","ProbeTracks_dz[BToMuEPiHD_sel_e_idx]","ProbeTracks_dzS[BToMuEPiHD_sel_e_idx]",
               "ProbeTracks_dxy[BToMuEPiHD_sel_e_idx]","ProbeTracks_dxyS[BToMuEPiHD_sel_e_idx]","ProbeTracks_DCASig[BToMuEPiHD_sel_e_idx]","ProbeTracks_nValidHits[BToMuEPiHD_sel_e_idx]",
               "ProbeTracks_pt[BToMuEPiHD_pi_idx]","ProbeTracks_eta[BToMuEPiHD_pi_idx]","ProbeTracks_phi[BToMuEPiHD_pi_idx]","ProbeTracks_dz[BToMuEPiHD_pi_idx]","ProbeTracks_dzS[BToMuEPiHD_pi_idx]",
               "ProbeTracks_dxy[BToMuEPiHD_pi_idx]","ProbeTracks_dxyS[BToMuEPiHD_pi_idx]","ProbeTracks_DCASig[BToMuEPiHD_pi_idx]"
		]
#tt_branches = ["BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_sv_lxy_sig"]
in_branches = ["BToMuMuPi_hnl_mass","BToMuMuPi_hnl_pt","BToMuMuPi_sv_lxy","BToMuMuPi_hnl_cos2D","BToMuMuPi_mass","BToMuMuPi_pt","BToMuMuPi_sv_lxy_sig"]
print str(sys.argv)
datatype = 0
#print datatype
print "initialize list of input"
if datatype == 1:
	filenames = glob.glob('../data/Mass_2tracks/*111.root')
	mu_sig= "BToMuMuPi_isMatched ==1"# &&  BToMuMuPi_hnl_charge ==0 "
	pf_sig= "BToMuEPi_isMatched ==1 &&  BToMuEPi_hnl_charge ==0 && Electron_isPF[BToMuEPi_sel_e_idx]"
	lowpT_sig= "BToMuEPi_isMatched &&  BToMuEPi_hnl_charge ==0 && Electron_isLowPt[BToMuEPi_sel_e_idx] && Electron_isPFoverlap[BToMuEPi_sel_e_idx]==0"
	tt_sig= "BToMuEPiHD_isMatched &&  BToMuEPiHD_hnl_charge ==0 "
elif datatype == 0:

	mu_sig= " (BToMuMuPi_mass<4.8 || BToMuMuPi_mass>5.7) && BToMuMuPi_hnl_charge ==0 "
	pf_sig="(BToMuEPi_mass<4.8 || BToMuEPi_mass>5.7) &&  BToMuEPi_hnl_charge ==0 && Electron_isPF[BToMuEPi_sel_e_idx]" 
	lowpT_sig="(BToMuEPi_mass<4.8 || BToMuEPi_mass>5.7) &&  BToMuEPi_hnl_charge ==0 && Electron_isLowPt[BToMuEPi_sel_e_idx] && Electron_isPFoverlap[BToMuEPi_sel_e_idx]==0 " 
	tt_sig= "(BToMuEPiHD_mass<4.8 || BToMuEPiHD_mass>5.7)&&  BToMuEPiHD_hnl_charge ==0 "
	filenames = glob.glob(str(sys.argv[2]))
#	sel= "(BToMuEPi_mass<2.2 || BToMuEPi_mass>7) && BToMuEPi_hnl_charge ==0"
#	filenames = glob.glob('/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/ParkingBPH1/crab_data_Run2018A_part1/210729_073310/0000*/*121.root')#/bparknano_56.root'
print "list of input completed"
#intree = filenames.Get('Events')

# and convert the TTree into an array
rta_mu =[] 
rta_muSel = []

rta_pf =[] 
rta_pfSel = []
rta_lowPt =[] 
rta_lowPtSel = []
rta_tt =[] 
rta_ttSel = []


mu_sel = " && TriggerMuon_pt[0]>0.3 &&  fabs(TriggerMuon_eta[0])<2.8 &&  Muon_pt[BToMuMuPi_sel_mu_idx]>0.3 &&  fabs(Muon_eta[BToMuMuPi_sel_mu_idx])<2 &&  ProbeTracks_pt[BToMuMuPi_pi_idx]>0.7  &&  fabs(ProbeTracks_eta[BToMuMuPi_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuMuPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuMuPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuMuPi_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuMuPi_pi_idx])>1.5  &&  BToMuMuPi_sv_prob>0.001 &&  BToMuMuPi_hnl_pt>1 &&  fabs(BToMuMuPi_hnl_cos2D)>0.9 &&  BToMuMuPi_sv_lxy/BToMuMuPi_sv_lxye>3  &&  fabs(Muon_dxy[BToMuMuPi_sel_mu_idx])>0.001 &&  fabs(Muon_dz[BToMuMuPi_sel_mu_idx])>0.0015 &&  fabs(Muon_dxy[BToMuMuPi_sel_mu_idx]/Muon_dxyErr[BToMuMuPi_sel_mu_idx])>1.5 &&  fabs(Muon_dz[BToMuMuPi_sel_mu_idx]/Muon_dzErr[BToMuMuPi_sel_mu_idx])>1  &&  fabs(ProbeTracks_eta[BToMuMuPi_pi_idx]-BToMuMuPi_fit_pi_eta )<0.015";
pf_sel = "  &&  TriggerMuon_pt[0]>7 &&  fabs(TriggerMuon_eta[0])<1.5 && Electron_pfmvaId[BToMuEPi_sel_e_idx]> -3  && ProbeTracks_pt[BToMuEPi_pi_idx]>0.7  &&  fabs(ProbeTracks_eta[BToMuEPi_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuEPi_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuEPi_pi_idx])>1.5  &&  fabs(ProbeTracks_DCASig[BToMuEPi_pi_idx])>5  &&  BToMuEPi_sv_prob>0.001 &&  fabs(BToMuEPi_hnl_cos2D)>0.99  &&  BToMuEPi_sv_lxy/BToMuEPi_sv_lxye>3  &&  fabs(Electron_dxy[BToMuEPi_sel_e_idx])>0.008 &&  fabs(Electron_dz[BToMuEPi_sel_e_idx])>0.008 &&  fabs(Electron_dxy[BToMuEPi_sel_e_idx]/Electron_dxyErr[BToMuEPi_sel_e_idx])>1.5 &&  fabs(Electron_dz[BToMuEPi_sel_e_idx]/Electron_dzErr[BToMuEPi_sel_e_idx])>1"# &&  deltaR_trgMuPi<2;
#lowpt missing cut on lowptmvaId
lowPt_sel = "  &&  TriggerMuon_pt[0]>7 &&  fabs(TriggerMuon_eta[0])<1.5 && Electron_mvaId[BToMuEPi_sel_e_idx]> -3  && ProbeTracks_pt[BToMuEPi_pi_idx]>0.7  &&  fabs(ProbeTracks_eta[BToMuEPi_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuEPi_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuEPi_pi_idx])>1.5  &&  fabs(ProbeTracks_DCASig[BToMuEPi_pi_idx])>5  &&  BToMuEPi_sv_prob>0.001 &&  fabs(BToMuEPi_hnl_cos2D)>0.99  &&  BToMuEPi_sv_lxy/BToMuEPi_sv_lxye>3  &&  fabs(Electron_dxy[BToMuEPi_sel_e_idx])>0.008 &&  fabs(Electron_dz[BToMuEPi_sel_e_idx])>0.008 &&  fabs(Electron_dxy[BToMuEPi_sel_e_idx]/Electron_dxyErr[BToMuEPi_sel_e_idx])>1.5 &&  fabs(Electron_dz[BToMuEPi_sel_e_idx]/Electron_dzErr[BToMuEPi_sel_e_idx])>1"# &&  deltaR_trgMuPi<2;

tt_sel = "  &&  TriggerMuon_pt[0]>7 &&  fabs(TriggerMuon_eta[0])<1.5 &&  fabs(ProbeTracks_eta[BToMuEPiHD_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuEPiHD_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuEPiHD_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuEPiHD_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuEPiHD_pi_idx])>1.5  &&  fabs(ProbeTracks_DCASig[BToMuEPiHD_pi_idx])>5  &&  BToMuEPiHD_sv_prob>0.001 &&  fabs(BToMuEPiHD_hnl_cos2D)>0.99  &&  BToMuEPiHD_sv_lxy/BToMuEPiHD_sv_lxye>3  &&  fabs(ProbeTracks_dxy[BToMuEPiHD_sel_e_idx])>0.008 &&  fabs(ProbeTracks_dz[BToMuEPiHD_sel_e_idx])>0.008 &&  fabs(ProbeTracks_dxyS[BToMuEPiHD_sel_e_idx])>1.5 &&  fabs(ProbeTracks_dzS[BToMuEPiHD_sel_e_idx])>1"# &&  deltaR_trgMuPi<2;
inputFiles = tqdm(filenames)
#ROOT.EnableImplicitMT(); # enable multi-threading
for j in range(0,len(lxyMu_bins)):
	print "root to array"
#	rta_mu.append(root2array(inputFiles,"Events",branches= mu_branches,object_selection={mu_sig+" && "+lxyMu_bins[j]:['BToMuMuPi_hnl_mass','BToMuMuPi_hnl_pt','BToMuMuPi_sv_lxy','BToMuMuPi_hnl_cos2D','BToMuMuPi_mass','BToMuMuPi_pt','BToMuMuPi_sv_lxy_sig']}))

 	rta_muSel.append(root2array(inputFiles,"Events",branches= mu_branches,object_selection={mu_sig+mu_sel+" && "+lxyMu_bins[j]:["BToMuMuPi_mass","BToMuMuPi_pt","BToMuMuPi_eta","BToMuMuPi_phi",
           
               "BToMuMuPi_hnl_mass","BToMuMuPi_hnl_pt","BToMuMuPi_sv_lxy","BToMuMuPi_hnl_cos2D","BToMuMuPi_sv_lxy_sig","BToMuMuPi_dr_mu_pi",
               "BToMuMuPi_dr_trgmu_hnl","BToMuMuPi_sv_prob","BToMuMuPi_sv_chi2","BToMuMuPi_sv_x","BToMuMuPi_sv_xe","BToMuMuPi_sv_y","BToMuMuPi_sv_ye","BToMuMuPi_sv_z","BToMuMuPi_sv_ze",
               "Muon_pt[BToMuMuPi_sel_mu_idx]","Muon_eta[BToMuMuPi_sel_mu_idx]","Muon_phi[BToMuMuPi_sel_mu_idx]","Muon_dz[BToMuMuPi_sel_mu_idx]","Muon_dzErr[BToMuMuPi_sel_mu_idx]",
               "Muon_dxy[BToMuMuPi_sel_mu_idx]","Muon_dxyErr[BToMuMuPi_sel_mu_idx]","Muon_matched_dpt[BToMuMuPi_sel_mu_idx]","Muon_matched_dr[BToMuMuPi_sel_mu_idx]",
               "ProbeTracks_pt[BToMuMuPi_pi_idx]","ProbeTracks_eta[BToMuMuPi_pi_idx]","ProbeTracks_phi[BToMuMuPi_pi_idx]","ProbeTracks_dz[BToMuMuPi_pi_idx]","ProbeTracks_dzS[BToMuMuPi_pi_idx]",
               "ProbeTracks_dxy[BToMuMuPi_pi_idx]","ProbeTracks_dxyS[BToMuMuPi_pi_idx]","ProbeTracks_DCASig[BToMuMuPi_pi_idx]"]
}))
 #	rta_pf.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+" && "+lxyEle_bins[j]:["BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_sv_lxy_sig"]}))
 	rta_pfSel.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={pf_sig+pf_sel+" && "+lxyEle_bins[j]:["BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_eta","BToMuEPi_phi",
               "BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_sv_lxy_sig","BToMuEPi_dr_mu_pi",
               "BToMuEPi_dr_trgmu_hnl","BToMuEPi_sv_prob","BToMuEPi_sv_chi2","BToMuEPi_sv_x","BToMuEPi_sv_xe","BToMuEPi_sv_y","BToMuEPi_sv_ye","BToMuEPi_sv_z","BToMuEPi_sv_ze",
               "Electron_pt[BToMuEPi_sel_e_idx]","Electron_eta[BToMuEPi_sel_e_idx]","Electron_phi[BToMuEPi_sel_e_idx]","Electron_dz[BToMuEPi_sel_e_idx]","Electron_dzErr[BToMuEPi_sel_e_idx]",
               "Electron_dxy[BToMuEPi_sel_e_idx]","Electron_dxyErr[BToMuEPi_sel_e_idx]","Electron_pfmvaId[BToMuEPi_sel_e_idx]","Electron_mvaId[BToMuEPi_sel_e_idx]",
               "ProbeTracks_pt[BToMuEPi_pi_idx]","ProbeTracks_eta[BToMuEPi_pi_idx]","ProbeTracks_phi[BToMuEPi_pi_idx]","ProbeTracks_dz[BToMuEPi_pi_idx]","ProbeTracks_dzS[BToMuEPi_pi_idx]",
               "ProbeTracks_dxy[BToMuEPi_pi_idx]","ProbeTracks_dxyS[BToMuEPi_pi_idx]","ProbeTracks_DCASig[BToMuEPi_pi_idx]"
		]
}))
 #	rta_`:wq
#lowPt.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={lowpT_sig+" && "+lxyEle_bins[j]:["BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_sv_lxy_sig"]}))
 	rta_lowPtSel.append(root2array(inputFiles,"Events",branches= ele_branches,object_selection={lowpT_sig+lowPt_sel+" && "+lxyEle_bins[j]:["BToMuEPi_mass","BToMuEPi_pt","BToMuEPi_eta","BToMuEPi_phi",
               "BToMuEPi_hnl_mass","BToMuEPi_hnl_pt","BToMuEPi_sv_lxy","BToMuEPi_hnl_cos2D","BToMuEPi_sv_lxy_sig","BToMuEPi_dr_mu_pi",
               "BToMuEPi_dr_trgmu_hnl","BToMuEPi_sv_prob","BToMuEPi_sv_chi2","BToMuEPi_sv_x","BToMuEPi_sv_xe","BToMuEPi_sv_y","BToMuEPi_sv_ye","BToMuEPi_sv_z","BToMuEPi_sv_ze",
               "Electron_pt[BToMuEPi_sel_e_idx]","Electron_eta[BToMuEPi_sel_e_idx]","Electron_phi[BToMuEPi_sel_e_idx]","Electron_dz[BToMuEPi_sel_e_idx]","Electron_dzErr[BToMuEPi_sel_e_idx]",
               "Electron_dxy[BToMuEPi_sel_e_idx]","Electron_dxyErr[BToMuEPi_sel_e_idx]","Electron_pfmvaId[BToMuEPi_sel_e_idx]","Electron_mvaId[BToMuEPi_sel_e_idx]",
               "ProbeTracks_pt[BToMuEPi_pi_idx]","ProbeTracks_eta[BToMuEPi_pi_idx]","ProbeTracks_phi[BToMuEPi_pi_idx]","ProbeTracks_dz[BToMuEPi_pi_idx]","ProbeTracks_dzS[BToMuEPi_pi_idx]",
               "ProbeTracks_dxy[BToMuEPi_pi_idx]","ProbeTracks_dxyS[BToMuEPi_pi_idx]","ProbeTracks_DCASig[BToMuEPi_pi_idx]"
		]
}))
 #	rta_tt.append(root2array(inputFiles,"Events",branches= tt_branches,object_selection={tt_sig+" && "+lxytt_bins[j]:["BToMuEPiHD_hnl_mass","BToMuEPiHD_hnl_pt","BToMuEPiHD_sv_lxy","BToMuEPiHD_hnl_cos2D","BToMuEPiHD_mass","BToMuEPiHD_pt","BToMuEPiHD_sv_lxy_sig"]}))
 	rta_ttSel.append(root2array(inputFiles,"Events",branches= tt_branches,object_selection={tt_sig+tt_sel+" && "+lxytt_bins[j]:["BToMuEPiHD_mass","BToMuEPiHD_pt","BToMuEPiHD_eta","BToMuEPiHD_phi",
               "BToMuEPiHD_hnl_mass","BToMuEPiHD_hnl_pt","BToMuEPiHD_sv_lxy","BToMuEPiHD_hnl_cos2D","BToMuEPiHD_sv_lxy_sig","BToMuEPiHD_dr_mu_pi",
               "BToMuEPiHD_dr_trgmu_hnl","BToMuEPiHD_sv_prob","BToMuEPiHD_sv_chi2","BToMuEPiHD_sv_x","BToMuEPiHD_sv_xe","BToMuEPiHD_sv_y","BToMuEPiHD_sv_ye","BToMuEPiHD_sv_z","BToMuEPiHD_sv_ze",
               "ProbeTracks_pt[BToMuEPiHD_sel_e_idx]","ProbeTracks_eta[BToMuEPiHD_sel_e_idx]","ProbeTracks_phi[BToMuEPiHD_sel_e_idx]","ProbeTracks_dz[BToMuEPiHD_sel_e_idx]","ProbeTracks_dzS[BToMuEPiHD_sel_e_idx]",
               "ProbeTracks_dxy[BToMuEPiHD_sel_e_idx]","ProbeTracks_dxyS[BToMuEPiHD_sel_e_idx]","ProbeTracks_DCASig[BToMuEPiHD_sel_e_idx]","ProbeTracks_nValidHits[BToMuEPiHD_sel_e_idx]",
               "ProbeTracks_pt[BToMuEPiHD_pi_idx]","ProbeTracks_eta[BToMuEPiHD_pi_idx]","ProbeTracks_phi[BToMuEPiHD_pi_idx]","ProbeTracks_dz[BToMuEPiHD_pi_idx]","ProbeTracks_dzS[BToMuEPiHD_pi_idx]",
               "ProbeTracks_dxy[BToMuEPiHD_pi_idx]","ProbeTracks_dxyS[BToMuEPiHD_pi_idx]","ProbeTracks_DCASig[BToMuEPiHD_pi_idx]"

		]
}))
 #	rta_charge.append(root2array(inputFiles,"Events",branches= in_branches,selection= "BToMuEPiHD_hnl_charge ==2 && (BToMuEPiHD_hnl_cos2D>0.99 ||BToMuEPiHD_hnl_cos2D<-0.99) "))

rta = []
rta.append(rta_muSel)
rta.append(rta_pfSel)
rta.append(rta_lowPtSel)
rta.append(rta_ttSel)
#num_cores = multiprocessing.cpu_count()

HistPath = "../hists/"
#iprint "printing array..."
#print rta_mu[0]
#print len(rta_muSel[0])
#print len(rta_pfSel[0])
#print len(rta_lowPtSel[0])
#print len(rta_ttSel[0])
#print len(rta_mu[0][0][0])
#print len(rta_mu[0][0][0][0])
#histMu = []
histMuSel = []
#histPF = []
histPFSel = []
#histLowpT = []
histLowpTSel = []
#histTT = []
histTTSel = []
#
##inputs = tqdm(rta_mu)
#
##processed_list = Parallel(n_jobs=num_cores)(delayed(fillHistos)(i) for i in rta_mu)
#
#
#
#histMu = createHistos(1, "BToMuMuPi")
#histMuSel= createHistos(0,"BToMuMuPi_sel")
#histPF = createHistos( 1, "BToMuPFPi")
#histPFSel = createHistos( 0,"BToMuPFPi_sel")
#histLowpT=createHistos(1,"BToMuLowPtPi")
#histLowpTSel= createHistos( 0,"BToMuLowPtPi_sel")
#histTT=createHistos(1,"BToMuTPi")
#histTTSel= createHistos( 0,"BToMuTPi_sel")
#print len(histMu)
#tree = ROOT.TTree("Events", "Events")

#lenght = []
#cat = [lenght, lenght, lenght, lenght, lenght, lenght,lenght]
#cat =  ROOT.std.vector(ROOT.std.vector('float'))()
#cat.resize(len(in_branches)) #branches definition for output tree
#for i in range (0,len(in_branches)):
#	temp = ROOT.std.vector('float')()
#	cat.push_back(temp)
#	tree._v = cat.at(i)
#	tree.Branch(in_branches[i],cat.at(i))
rfileSel = root_open('BToMuLPiSmall_Chunk'+str(sys.argv[1])+'.root', 'w')
treeSel = Tree("Events", model=Event)

for k in range(0,len(rta[0])):
  for j in range(0,len(rta[0][0])):
    for SigIdx in range (0,4):
  #   if (len(rta_mu[k][j][0])!=0 ):
  #    # print "loop mu"
  #    # fillHistos(rta_mu[k][j], histMu[k],1, "BToMuMuPi")
# #     test = zip(*rta_mu[k][j])
# #     test = zip(*test)
  #     tree.evtIdx =j 
  #     tree.nCand=len(rta_mu[k][j][0]) 
  #     tree.Type=0
  #     tree.LxyBin=k
# #      for idx in range(0,len(in_branches)):
  #     for l in range(0,len(rta_mu[k][j][0])):	
  #       	tree.hnl_mass[l] =  rta_mu[k][j][0][l]
  #      	tree.hnl_pt[l] =  rta_mu[k][j][1][l]
# #      	print rta_mu[k][j][1][l]

  #      	tree.hnl_lxy[l] =  rta_mu[k][j][2][l]
  #      	tree.hnl_cos2D[l] =  rta_mu[k][j][3][l]
  #      	tree.B_mass[l] =  rta_mu[k][j][4][l]
  #      	tree.B_pt[l] =  rta_mu[k][j][5][l]
  #      	tree.hnl_lxy_sig[l] =  rta_mu[k][j][6][l]
  #   else:
 
  #    if (len(rta_pf[k][j][0])!=0 ):
  #     tree.evtIdx =j 
  #     tree.nCand=len(rta_pf[k][j][0]) 
  #     tree.Type=1
  #     tree.LxyBin=k
  #     for l in range(0,len(rta_pf[k][j][0])):	
  #       	tree.hnl_mass[l] =  rta_pf[k][j][0][l]
  #      	tree.hnl_pt[l] =  rta_pf[k][j][1][l]
  #      	tree.hnl_lxy[l] =  rta_pf[k][j][2][l]
  #      	tree.hnl_cos2D[l] =  rta_pf[k][j][3][l]
  #      	tree.B_mass[l] =  rta_pf[k][j][4][l]
  #      	tree.B_pt[l] =  rta_pf[k][j][5][l]
  #      	tree.hnl_lxy_sig[l] =  rta_pf[k][j][6][l]
  #    # print "loop pf"
  #    # fillHistos(rta_pf[k][j], histPF[k],1, "BToMuPFPi")
  #    else:
 
  #     if (len(rta_lowPt[k][j][0])!=0 ):
  #      tree.evtIdx =j 
  #      tree.nCand=len(rta_lowPt[k][j][0]) 
  #      tree.Type=2
  #      tree.LxyBin=k
  #      
  #      for l in range(0,len(rta_lowPt[k][j][0])):	
  #       	tree.hnl_mass[l] =  rta_lowPt[k][j][0][l]
  #      	tree.hnl_pt[l] =  rta_lowPt[k][j][1][l]
  #      	tree.hnl_lxy[l] =  rta_lowPt[k][j][2][l]
  #      	tree.hnl_cos2D[l] =  rta_lowPt[k][j][3][l]
  #      	tree.B_mass[l] =  rta_lowPt[k][j][4][l]
  #      	tree.B_pt[l] =  rta_lowPt[k][j][5][l]
  #      	tree.hnl_lxy_sig[l] =  rta_lowPt[k][j][6][l]
  #   #  print "loop pf"
  #    #  fillHistos(rta_lowPt[k][j], histLowpT[k],1,"BToMuLowPtPi")
  #     else:
  #      	 
  #      if (len(rta_tt[k][j][0])!=0 ):
  #       tree.evtIdx =j 
  #       tree.nCand=len(rta_tt[k][j][0]) 
  #       tree.Type=3
  #       tree.LxyBin=k
  #       for l in range(0,len(rta_tt[k][j][0])):	
  #        	tree.hnl_mass[l] =  rta_tt[k][j][0][l]
  #       	tree.hnl_pt[l] =  rta_tt[k][j][1][l]
  #       	tree.hnl_lxy[l] =  rta_tt[k][j][2][l]
  #       	tree.hnl_cos2D[l] =  rta_tt[k][j][3][l]
  #       	tree.B_mass[l] =  rta_tt[k][j][4][l]
  #       	tree.B_pt[l] =  rta_tt[k][j][5][l]
  #       	tree.hnl_lxy_sig[l] =  rta_tt[k][j][6][l]
  #    #  print "loop pf"
  #    #   fillHistos(rta_tt[k][j], histTT[k],1,"BToMuTPi")

#trees after selection
      if ( len(rta[SigIdx][k][j][0])!=0):
#	print("for event %d and displacement %d, signal category %d",k,j,SigIdx)
        sortBranches(rta[SigIdx][k][j],0, sig[SigIdx])
        treeSel.evtIdx =j 
        treeSel.nCand=1
        treeSel.Type=SigIdx
        treeSel.LxyBin=k
        treeSel.B_mass= rta[SigIdx][k][j][0][0]                                                      
        treeSel.B_pt= rta[SigIdx][k][j][1][0]                                                        
        treeSel.B_eta= rta[SigIdx][k][j][2][0]                                                       
        treeSel.B_phi= rta[SigIdx][k][j][3][0]                                                       
        treeSel.TrgMu_pt= rta[SigIdx][k][j][4][0]                                                    
        treeSel.TrgMu_eta= rta[SigIdx][k][j][5][0]                                                   
        treeSel.TrgMu_phi= rta[SigIdx][k][j][6][0]                                                   
        treeSel.hnl_mass= rta[SigIdx][k][j][7][0]                                                    
        treeSel.hnl_pt= rta[SigIdx][k][j][8][0]                                                      
        treeSel.hnl_lxy= rta[SigIdx][k][j][9][0]                                                     
        treeSel.hnl_cos2D= rta[SigIdx][k][j][10][0]                                                  
        treeSel.hnl_lxy_sig = rta[SigIdx][k][j][11][0]                                               
        treeSel.hnl_drLepPi = rta[SigIdx][k][j][12][0]                                               
        treeSel.dr_trgmu_hnl = rta[SigIdx][k][j][13][0]                                              
        treeSel.hnl_vtxProb= rta[SigIdx][k][j][14][0]                                                
        treeSel.hnl_vtxChi2= rta[SigIdx][k][j][15][0]                                                
        treeSel.hnl_vx =  rta[SigIdx][k][j][16][0]                                                   
        treeSel.hnl_ex =  rta[SigIdx][k][j][17][0]                                                   
        treeSel.hnl_vy =  rta[SigIdx][k][j][18][0]                                                   
        treeSel.hnl_ey =  rta[SigIdx][k][j][19][0]                                                   
        treeSel.hnl_vz = rta[SigIdx][k][j][20][0]                                                    
        treeSel.hnl_ez = rta[SigIdx][k][j][21][0]                                                    
        treeSel.hnl_l_pt = rta[SigIdx][k][j][22][0]                                                  
        treeSel.hnl_l_eta =rta[SigIdx][k][j][23][0]                                                  
        treeSel.hnl_l_phi =rta[SigIdx][k][j][24][0]                                                  
        treeSel.hnl_l_dz = rta[SigIdx][k][j][25][0]   
	if (SigIdx <3):
						                                              
          treeSel.hnl_l_dzS =rta[SigIdx][k][j][25][0]/rta[SigIdx][k][j][26][0]
	else:                                             
          treeSel.hnl_l_dzS =rta[SigIdx][k][j][26][0]
        treeSel.hnl_l_dxy =rta[SigIdx][k][j][27][0]                                                  
	if (SigIdx <3):
						                                               
          treeSel.hnl_l_dxyS =rta[SigIdx][k][j][27][0]/rta[SigIdx][k][j][28][0]
	else:                                            
          treeSel.hnl_l_dxyS =rta[SigIdx][k][j][28][0]
	
	if (SigIdx <3):
          treeSel.hnl_l_DCAS = -99
	else:  
          treeSel.hnl_l_DCAS =rta[SigIdx][k][j][29][0] 

	if (SigIdx ==1):
         treeSel.hnl_l_mvaId= rta[SigIdx][k][j][29][0]    
  	elif (SigIdx ==2):
         treeSel.hnl_l_mvaId= rta[SigIdx][k][j][30][0]    
  	else:
         treeSel.hnl_l_mvaId= -99
        treeSel.hnl_pi_pt =rta[SigIdx][k][j][31][0]                                                  
        treeSel.hnl_pi_eta= rta[SigIdx][k][j][32][0]                                                 
        treeSel.hnl_pi_phi = rta[SigIdx][k][j][33][0]                                                
        treeSel.hnl_pi_dz =rta[SigIdx][k][j][34][0]                                                  
        treeSel.hnl_pi_dzS = rta[SigIdx][k][j][35][0]                                                
        treeSel.hnl_pi_dxy= rta[SigIdx][k][j][36][0]                                                 
        treeSel.hnl_pi_dxyS = rta[SigIdx][k][j][37][0]                                      
        treeSel.hnl_pi_DCAS = rta[SigIdx][k][j][38][0]                                      
        treeSel.fill()
      else:
	continue
 
##saveHistos(rta_mu, histMu,1, "BToMuMuPi")
#saveHistos(rta_muSel, histMuSel,0,"BToMuMuPi_sel")
#saveHistos(rta_pf, histPF,1, "BToMuPFPi")
#saveHistos(rta_pfSel, histPFSel,0,"BToMuPFPi_sel")
#saveHistos(rta_lowPt, histLowpT,1,"BToMuLowPtPi")
#saveHistos(rta_lowPtSel, histLowpTSel,0,"BToMuLowPtPi_sel")
#saveHistos(rta_tt, histTT,1,"BToMuTPi")
#saveHistos(rta_ttSel, histTTSel,1,"BToMuTPi")
#ree.Write()
treeSel.write()
rfileSel.close()

#array2root(rta_mu[0],"test.root","recreate")
