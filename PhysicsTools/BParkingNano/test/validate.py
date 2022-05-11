from ROOT import gROOT, TH1F, TH1, kBlack, kRed, kBlue, TCanvas, gStyle, TLegend, TLatex, TFile, TChain, TPad, kWhite, TMath
import ROOT
import os
import os.path

quantities = [
{'hname':'BToMuMuPi_hnl_pt', 'title':'hnl pT [GeV]', 'binning':(90, 0, 35)},
{'hname':'BToMuMuPi_hnl_eta', 'title':'hnl eta', 'binning':(50, -3.5, 3.5)},
{'hname':'BToMuMuPi_hnl_phi', 'title':'hnl phi', 'binning':(50, -3.2, 3.2)},
{'hname':'BToMuMuPi_hnl_mass', 'title':'hnl mass [GeV]', 'binning':(90, 0, 6.5)},
{'hname':'BToMuMuPi_hnl_cos2D', 'title':'hnl cos2D [GeV]', 'binning':(90, 0.9, 1)},
{'hname':'BToMuMuPi_hnl_iso03', 'title':'hnl iso03', 'binning':(50, 0, 50)},
{'hname':'BToMuMuPi_hnl_iso04', 'title':'hnl iso04', 'binning':(50, 0, 50)},

#{'hname':'hnl_pt', 'title':'hnl pT [GeV]', 'binning':(90, 0, 35)},
#{'hname':'hnl_eta', 'title':'hnl eta', 'binning':(50, -3.5, 3.5)},
#{'hname':'hnl_phi', 'title':'hnl phi', 'binning':(50, -3.2, 3.2)},
#{'hname':'hnl_mass', 'title':'hnl mass [GeV]', 'binning':(90, 0, 6.5)},
#{'hname':'hnl_cos2D', 'title':'hnl cos2D [GeV]', 'binning':(90, 0.9, 1)},

{'hname':'BToMuMuPi_pt', 'title':'b pT [GeV]', 'binning':(90, 0, 65)},
{'hname':'BToMuMuPi_eta', 'title':'b eta', 'binning':(50, -3.5, 3.5)},
{'hname':'BToMuMuPi_phi', 'title':'b phi', 'binning':(50, -3.2, 3.2)},
{'hname':'BToMuMuPi_mass', 'title':'b mass [GeV]', 'binning':(90, 0, 8)},

{'hname':'BToMuMuPi_trg_mu_pt', 'title':'trigger muon pT [GeV]', 'binning':(90, 0, 30)},
{'hname':'BToMuMuPi_trg_mu_eta', 'title':'trigger muon eta', 'binning':(50, -3.5, 3.5)},
{'hname':'BToMuMuPi_trg_mu_phi', 'title':'trigger muon phi', 'binning':(50, -3.2, 3.2)},
{'hname':'BToMuMuPi_trg_mu_ip3d', 'title':'trigger muon ip3d [cm]', 'binning':(50, 0, 1.5)},
{'hname':'BToMuMuPi_trg_mu_sip3d', 'title':'trigger muon sip3d', 'binning':(50, 0, 300)},
{'hname':'BToMuMuPi_trg_mu_dxy', 'title':'trigger muon dxy [cm]', 'binning':(50, -1.5, 1.5)},
{'hname':'BToMuMuPi_trg_mu_dz', 'title':'trigger muon dz [cm]', 'binning':(50, -1.5, 1.5)},
{'hname':'BToMuMuPi_trg_mu_iso03', 'title':'trigger muon iso03', 'binning':(50, 0, 50)},
{'hname':'BToMuMuPi_trg_mu_iso04', 'title':'trigger muon iso04', 'binning':(50, 0, 50)},

{'hname':'BToMuMuPi_fit_mu_pt', 'title':'muon pT [GeV]', 'binning':(90, 0, 35)},
{'hname':'BToMuMuPi_fit_mu_eta', 'title':'muon eta', 'binning':(50, -3.5, 3.5)},
{'hname':'BToMuMuPi_fit_mu_phi', 'title':'muon phi', 'binning':(50, -3.2, 3.2)},
{'hname':'BToMuMuPi_sel_mu_ip3d', 'title':'muon ip3d [cm]', 'binning':(50, 0, 1.5)},
{'hname':'BToMuMuPi_sel_mu_sip3d', 'title':'muon sip3d', 'binning':(50, 0, 300)},
{'hname':'BToMuMuPi_sel_mu_dxy', 'title':'muon dxy [cm]', 'binning':(50, -1.5, 1.5)},
{'hname':'BToMuMuPi_sel_mu_dz', 'title':'muon dz [cm]', 'binning':(50, -1.5, 1.5)},
{'hname':'BToMuMuPi_sel_mu_iso03', 'title':'muon iso03', 'binning':(50, 0, 50)},
{'hname':'BToMuMuPi_sel_mu_iso04', 'title':'muon iso04', 'binning':(50, 0, 50)},

{'hname':'BToMuMuPi_fit_pi_pt', 'title':'pion pT [GeV]', 'binning':(90, 0, 30)},
{'hname':'BToMuMuPi_fit_pi_eta', 'title':'pion eta', 'binning':(50, -3.5, 3.5)},
{'hname':'BToMuMuPi_fit_pi_phi', 'title':'pion phi', 'binning':(50, -3.2, 3.2)},
{'hname':'BToMuMuPi_pi_dxy', 'title':'pion dxy [cm]', 'binning':(50, -3, 3)},
{'hname':'BToMuMuPi_pi_dz', 'title':'pion dz [cm]', 'binning':(50, -3, 3)},
{'hname':'BToMuMuPi_pi_dxyS', 'title':'pion dxyS [cm]', 'binning':(50, -3, 3)},
{'hname':'BToMuMuPi_pi_dzS', 'title':'pion dzS [cm]', 'binning':(50, -3, 3)},
{'hname':'BToMuMuPi_pi_DCASig', 'title':'pion DCASig [cm]', 'binning':(50, 0, 30)},
{'hname':'BToMuMuPi_pi_iso03', 'title':'pion iso03', 'binning':(50, 0, 50)},
{'hname':'BToMuMuPi_pi_iso04', 'title':'pion iso04', 'binning':(50, 0, 50)},

{'hname':'BToKMuMu_fit_pt', 'title':'BToKMuMu b pT [GeV]', 'binning':(90, 0, 65)},
{'hname':'BToKMuMu_fit_eta', 'title':'BToKMuMu b eta', 'binning':(50, -3.5, 3.5)},
{'hname':'BToKMuMu_fit_phi', 'title':'BToKMuMu b phi', 'binning':(50, -3.2, 3.2)},
{'hname':'BToKMuMu_fit_mass', 'title':'BToKMuMu b mass [GeV]', 'binning':(90, 0, 8)},
{'hname':'BToKMuMu_fit_cos2D', 'title':'BToKMuMu cos2D', 'binning':(90, 0.9, 1)},
{'hname':'BToKMuMu_b_iso03', 'title':'BToKMuMu b iso03', 'binning':(50, 0, 50)},
{'hname':'BToKMuMu_b_iso04', 'title':'BToKMuMu b iso04', 'binning':(50, 0, 50)},

{'hname':'BToKMuMu_fit_k_pt', 'title':'BToKMuMu fit k pT [GeV]', 'binning':(90, 0, 65)},
{'hname':'BToKMuMu_fit_k_eta', 'title':'BToKMuMu fit k eta', 'binning':(50, -3.5, 3.5)},
{'hname':'BToKMuMu_fit_k_phi', 'title':'BToKMuMu fit k phi', 'binning':(50, -3.2, 3.2)},
{'hname':'BToKMuMu_k_iso03', 'title':'BToKMuMu k iso03', 'binning':(50, 0, 50)},
{'hname':'BToKMuMu_k_iso04', 'title':'BToKMuMu k iso04', 'binning':(50, 0, 50)},

{'hname':'BToKMuMu_fit_l1_pt', 'title':'BToKMuMu fit l1 pT [GeV]', 'binning':(90, 0, 65)},
{'hname':'BToKMuMu_fit_l1_eta', 'title':'BToKMuMu fit l1 eta', 'binning':(50, -3.5, 3.5)},
{'hname':'BToKMuMu_fit_l1_phi', 'title':'BToKMuMu fit l1 phi', 'binning':(50, -3.2, 3.2)},
{'hname':'BToKMuMu_l1_iso03', 'title':'BToKMuMu l1 iso03', 'binning':(50, 0, 50)},
{'hname':'BToKMuMu_l1_iso04', 'title':'BToKMuMu l1 iso04', 'binning':(50, 0, 50)},

{'hname':'BToKMuMu_fit_l2_pt', 'title':'BToKMuMu fit l2 pT [GeV]', 'binning':(90, 0, 65)},
{'hname':'BToKMuMu_fit_l2_eta', 'title':'BToKMuMu fit l2 eta', 'binning':(50, -3.5, 3.5)},
{'hname':'BToKMuMu_fit_l2_phi', 'title':'BToKMuMu fit l2 phi', 'binning':(50, -3.2, 3.2)},
{'hname':'BToKMuMu_l2_iso03', 'title':'BToKMuMu l2 iso03', 'binning':(50, 0, 50)},
{'hname':'BToKMuMu_l2_iso04', 'title':'BToKMuMu l2 iso04', 'binning':(50, 0, 50)},

{'hname':'BToKMuMu_mll_fullfit', 'title':'BToKMuMu_mll_fullfit', 'binning':(50, 0, 10)},
{'hname':'BToKMuMu_l_xy', 'title':'BToKMuMu l_xy', 'binning':(90, 0, 1)},

{'hname':'Pileup_nPU', 'title':'Pileup_nPU', 'binning':(90, 0, 90)},
{'hname':'PV_npvs', 'title':'PV_npvs', 'binning':(90, 0, 90)},
{'hname':'PV_npvsGood', 'title':'PV_npvsGood', 'binning':(90, 0, 90)},
#{'hname':'fixedGridRhoFastjetAll', 'title':'rho', 'binning':(90, 0, 50)},

{'hname':'TriggerMuon_pt', 'title':'TriggerMuon pT [GeV]', 'binning':(90, 0, 30)},
{'hname':'TriggerMuon_eta', 'title':'TriggerMuon eta', 'binning':(50, -3.5, 3.5)},
{'hname':'TriggerMuon_phi', 'title':'TriggerMuon phi', 'binning':(50, -3.2, 3.2)},
{'hname':'TriggerMuon_vx', 'title':'TriggerMuon vx', 'binning':(45, 0, 15)},
{'hname':'TriggerMuon_vy', 'title':'TriggerMuon vy', 'binning':(45, 0, 15)},
{'hname':'TriggerMuon_vz', 'title':'TriggerMuon vz', 'binning':(45, 0, 15)},

{'hname':'Muon_pt', 'title':'Muon pT [GeV]', 'binning':(90, 0, 30)},
{'hname':'Muon_eta', 'title':'Muon eta', 'binning':(50, -3.5, 3.5)},
{'hname':'Muon_phi', 'title':'Muon phi', 'binning':(50, -3.2, 3.2)},
{'hname':'Muon_vx', 'title':'Muon vx', 'binning':(45, 0, 15)},
{'hname':'Muon_vy', 'title':'Muon vy', 'binning':(45, 0, 15)},
{'hname':'Muon_vz', 'title':'Muon vz', 'binning':(45, 0, 15)},
{'hname':'Muon_dz', 'title':'Muon dz', 'binning':(100, 0, 0.1)},
{'hname':'Muon_dzS', 'title':'Muon dzS', 'binning':(100, 0, 0.1)},
{'hname':'Muon_dxy', 'title':'Muon dxy', 'binning':(100, 0, 0.1)},
{'hname':'Muon_dxyS', 'title':'Muon dxyS', 'binning':(100, 0, 0.1)},
{'hname':'Muon_isPF', 'title':'Muon isPF', 'binning':(3, -1, 1)},
{'hname':'Muon_isGlobalMuon', 'title':'Muon isGlobalMuon', 'binning':(3, -1, 1)},
{'hname':'Muon_isTrackerMuon', 'title':'Muon isTrackerMuon', 'binning':(3, -1, 1)},
{'hname':'Muon_looseId', 'title':'Muon looseId', 'binning':(3, -1, 1)},
{'hname':'Muon_softId', 'title':'Muon softId', 'binning':(3, -1, 1)},
{'hname':'Muon_mediumId', 'title':'Muon mediumId', 'binning':(3, -1, 1)},
{'hname':'Muon_tightId', 'title':'Muon tightId', 'binning':(3, -1, 1)},

{'hname':'ProbeTracks_pt', 'title':'ProbeTracks pT [GeV]', 'binning':(90, 0, 30)},
{'hname':'ProbeTracks_eta', 'title':'ProbeTracks eta', 'binning':(50, -3.5, 3.5)},
{'hname':'ProbeTracks_phi', 'title':'ProbeTracks phi', 'binning':(50, -3.2, 3.2)},
{'hname':'ProbeTracks_vx', 'title':'ProbeTracks vx', 'binning':(45, 0, 15)},
{'hname':'ProbeTracks_vy', 'title':'ProbeTracks vy', 'binning':(45, 0, 15)},
{'hname':'ProbeTracks_vz', 'title':'ProbeTracks vz', 'binning':(45, 0, 15)},

{'hname':'SV_pt', 'title':'SV pT [GeV]', 'binning':(90, 0, 30)},
{'hname':'SV_chi2', 'title':'SV chi2', 'binning':(15, 0, 7)},
{'hname':'SV_x', 'title':'SV x', 'binning':(45, 0, 15)},
{'hname':'SV_y', 'title':'SV y', 'binning':(45, 0, 15)},
{'hname':'SV_z', 'title':'SV z', 'binning':(45, 0, 15)},
]


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Validation of pfcluster analysis', add_help=True)

  parser.add_argument('--f1', type=str, dest='fname1', help='', default = '../bparknano.root')
  parser.add_argument('--f2', type=str, dest='fname2', help='', default = '../bparknano.root')
  parser.add_argument('--t1', '--tree1', type=str, dest='treename1', help='', default = 'Events')
  parser.add_argument('--t2', '--tree2', type=str, dest='treename2', help='', default = 'Events')
  parser.add_argument('--l1', type=str, dest='label1', help='', default = 'A')
  parser.add_argument('--l2', type=str, dest='label2', help='', default = 'B')
  parser.add_argument('--o',  type=str, dest='outdirname', help='output dir', default = 'test')
  parser.add_argument('--doShape', dest='doShape', help='Do shape comparison', action='store_true', default=False)
  parser.add_argument('--doLog', dest='doLog', help='Use Log scale on y axis', action='store_true', default=False)

  return parser.parse_args()


def getOverflowedHisto(h):
  htemp = h.Clone()
  nbins = htemp.GetNbinsX()
  last_plus_overflow = htemp.GetBinContent(nbins) + htemp.GetBinContent(nbins+1)
  last_plus_overflow_error = TMath.Sqrt( htemp.GetBinError(nbins)*htemp.GetBinError(nbins)  + htemp.GetBinError(nbins+1)*htemp.GetBinError(nbins+1))
  htemp.SetBinContent(nbins,last_plus_overflow )
  htemp.SetBinError(nbins,last_plus_overflow_error)
  return htemp


def makeHistoFromNtuple(infilename, intreename, outhistoname, outhistobinning, outhistoquantity, outhistoselection="(1)", outhistosmooth=False, addOverflow=True):
  TH1.SetDefaultSumw2()
  histo = TH1F(outhistoname, outhistoname, *(outhistobinning))
  chain = TChain(intreename)
  chain.Add(infilename)
  ret = chain.Project(histo.GetName(), outhistoquantity, '('+outhistoselection + ')') #+ ')*(' + outhistoweight + ')')

  if outhistosmooth:
    histo.Scale(1, 'width')

  if addOverflow:
    hret = getOverflowedHisto(histo)
  else:
    hret = histo

  return ret,hret


def makeHistSettings(h):
  h.GetYaxis().SetNdivisions(510);
  h.GetYaxis().SetLabelSize(0.06);
  h.GetXaxis().SetLabelSize(0.06);
  h.GetYaxis().SetTitleSize(0.06);
  h.GetXaxis().SetTitleSize(0.08);
  h.GetYaxis().SetTitleOffset(0.96);
  h.GetXaxis().SetTitleOffset(1.0);
  h.GetYaxis().SetLabelOffset(0.01);
  h.GetXaxis().SetLabelOffset(0.01);
  h.GetXaxis().SetTickLength(0.06);
  #h.SetMarkerStyle(1)


def makeRatioSettings(ratioMC):
  #ratioMC.GetYaxis().SetRangeUser(0.75,1.25)
  ratioMC.GetYaxis().SetRangeUser(0.0,2.0)
  ratioMC.GetYaxis().SetNdivisions(504);
  #ratioMC.GetYaxis().SetNdivisions(508);
  ratioMC.GetYaxis().SetLabelSize(0.15);
  ratioMC.GetXaxis().SetLabelSize(0.15);
  ratioMC.GetYaxis().SetTitleSize(0.15);
  ratioMC.GetXaxis().SetTitleSize(0.15);
  ratioMC.GetYaxis().SetTitleOffset(0.35);
  ratioMC.GetXaxis().SetTitleOffset(1.2);
  ratioMC.GetYaxis().SetLabelOffset(0.01);
  ratioMC.GetXaxis().SetLabelOffset(0.03);
  ratioMC.GetXaxis().SetTickLength(0.06);
  ratioMC.SetFillStyle(3004)
  ratioMC.SetLineWidth(2)


def makeRatioPlot(hNum, hDen, hDen2="", nameNum="num", nameDen="den", nameDen2="", xtitle="pt",ytitle="Entries", ratiotitle="Ratio", shape=False, log=True, plotName="ratio", outDir='out'):
  TH1.SetDefaultSumw2()

  # prepare settings of histos
  hNum.SetLineColor(kBlack)
  hNum.SetLineWidth(2)
  hNum.SetMarkerStyle(20)
  #hNum.SetMarkerStyle(1)
  hNum.SetMarkerColor(kBlack)
  hNum.GetYaxis().SetTitle(ytitle)

  hDen.SetLineColor(kRed+1)
  hDen.SetMarkerColor(kRed+1)
  hDen.SetLineWidth(2)
  hDen.SetMarkerStyle(21)
  #hDen.SetMarkerStyle(1)

  if nameDen2 != "":
    hDen2.SetLineColor(kBlue)
    hDen2.SetMarkerColor(kBlue)
    hDen2.SetLineWidth(2)
    hDen2.SetMarkerStyle(22)
    #hDen2.SetMarkerStyle(1)
    makeHistSettings(hDen2)#

  makeHistSettings(hNum)
  makeHistSettings(hDen)

  # prepare canva
  canvas=TCanvas(plotName, plotName, 600, 600)
  ROOT.SetOwnership(canvas, False) # Crucial to avoid crashes due to the way python deletes the objects
  canvas.cd()
  yMinP1=0.305;
  bottomMarginP1=0.005;
  pad1 = TPad('pad1','pad1',0,yMinP1,0.99,1)
  if log: pad1.SetLogy()
  pad1.SetBottomMargin(bottomMarginP1)
  pad1.SetFillColor(kWhite)
  pad1.SetTickx()
  pad1.SetTicky()
  pad2 = TPad('pad2','pad2',0.,0.01,.99,0.300)
  pad2.SetTopMargin(0.010)
  pad2.SetBottomMargin(0.40)
  pad2.SetFillColor(kWhite)
  pad2.SetGridy()
  pad1.SetNumber(1)
  pad2.SetNumber(2)
  pad1.Draw()
  pad2.Draw()

  # prepare legend
  leg = TLegend(0.65,0.68,0.82,0.88,'')
  leg.SetBorderSize(0)
  leg.SetTextSize(0.05)

  # Draw
  pad1.cd()
  histo_saver = []

  #hNumNorm = hNum.Clone()
  if shape:
    hNumNorm = hNum.DrawNormalized('LP')
    histo_saver.append(hNumNorm)
  else:
    hNum.Draw('PE')

  if log:
    hNum.SetMaximum(hNum.GetMaximum()*4)
    if shape:
      hNumNorm.SetMaximum(hNumNorm.GetMaximum()*4)
  else:
    hNum.SetMaximum(hNum.GetMaximum()*2)
    if shape:
      hNumNorm.SetMaximum(hNumNorm.GetMaximum()*2)
  #leg.Draw('same')

  #hDenNorm = hDen.Clone()
  if shape:
    hDenNorm = hDen.DrawNormalized('samePE')
    histo_saver.append(hDenNorm)
  else:
    hDen.Draw('samePE')

  #hDenNorm2 = hDen2.Clone()
  if nameDen2 != "":
    if shape:
      hDenNorm2 = hDen2.DrawNormalized('samePE')
      histo_saver.append(hDenNorm2)
    else:
      hDen2.Draw('samePE')
    leg.AddEntry(hDen2, nameDen2, 'PE')

  leg.AddEntry(hDen, nameDen, 'PE')
  leg.AddEntry(hNum, nameNum, 'PE')
  leg.Draw('same')

  #print hNumNorm.Integral(), hDenNorm.Integral(), hDenNorm2.Integral()
  #print hNum.Integral(), hDen.Integral(), hDen2.Integral()

  #######################################################
  ### RATIO PAD
  #######################################################
  pad2.cd()
  hRatio = hNum.Clone()
  if nameDen2 != "": hRatio2 = hNum.Clone()

  if shape:
    hRatio = hNumNorm.Clone()
    hRatio2 = hNumNorm.Clone()

    hRatio.Divide(hDenNorm)
    if nameDen2 != "": hRatio2.Divide(hDenNorm2)

  else:
    hRatio.Divide(hDen)
    if nameDen2 != "": hRatio2.Divide(hDen2)

  #hRatio.SetLineColor(kRed+2)
  if nameDen2 != "": hRatio2.SetLineColor(kBlue)
  makeRatioSettings(hRatio)
  if nameDen2 != "": makeRatioSettings(hRatio2)

  hRatio.GetXaxis().SetTitle(xtitle)
  hRatio.GetYaxis().SetTitle(ratiotitle)
  hRatio.SetTitle('')

  hRatio.Draw('PE')
  if nameDen2 != "": hRatio2.Draw('PEsame')
  #hRatio.GetXaxis().SetRangeUser(200.,2000.)

  canvas.SaveAs('{d}/{name}.pdf'.format(d=outDir, name = plotName))
  canvas.SaveAs('{d}/{name}.png'.format(d=outDir, name = plotName))


if __name__ == "__main__":

  options = getOptions()

  #ROOT.gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  ROOT.gROOT.ProcessLine('setTDRStyle()')
  ROOT.gROOT.SetBatch(True)
  
  import os
  os.system('mkdir {d}'.format(d=options.outdirname))

  if os.path.isfile(options.fname1)==False or os.path.isfile(options.fname2)==False: raise RuntimeError('One of the two files is not available, \n{} \n{}'.format(options.fname1, options.fname2))

  for q in quantities:
    print '\n\n=> Investigating quantity', q['title']
    #q['name'] = q['hname'].split('/')[1]
    q['name'] = q['hname']
    ret1,histo1 = makeHistoFromNtuple(infilename=options.fname1, intreename=options.treename1, outhistoname=q['name']+'_1', outhistobinning=q['binning'], outhistoquantity=q['name'])
    ret2,histo2 = makeHistoFromNtuple(infilename=options.fname2, intreename=options.treename2, outhistoname=q['name']+'_2', outhistobinning=q['binning'], outhistoquantity=q['name'])
   
    # now you should have the histograms
    if ret1 != -1 and ret2 != -1:
      # do the ratio plot
      ytitle = 'Entries' #if not options.doNorm else 'Entries (norm)'
      makeRatioPlot(hNum=histo1, hDen=histo2, nameNum=options.label1, nameDen=options.label2, xtitle=q['title'],ytitle=ytitle, ratiotitle="Ratio", shape=options.doShape, log=options.doLog, outDir=options.outdirname, plotName=q['name'])

      # entries
      print '     Entries histo1', histo1.GetEntries(), ' histo2', histo2.GetEntries()

    else:
      print '     Skipping ', q['name']




