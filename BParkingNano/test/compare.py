import ROOT
import os
import sys
#sys.path.append('/work/mratti/plotting/myplotting')
#from spares import *

if __name__ == "__main__":

  #### ROOT Options
  ROOT.gROOT.SetBatch(True)
  ROOT.TH1.SetDefaultSumw2()
  ROOT.TH1.StatOverflows(ROOT.kTRUE) # consider overflows for mean and rms calculation
  ROOT.gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  ROOT.gROOT.ProcessLine('setTDRStyle()')
  #gStyle.SetTitleXOffset(1.1);
  #gStyle.SetTitleYOffset(1.45);

  path_1 = 'LxyCompTwoDir_imrovedMatchEff_v2/bhnl_recoeff.root'
  path_2 = 'LxyCompTwoDir_imrovedMatchEff_v2/bhnl_recoeff_dg.root'
  path_3 = 'LxyCompTwoDir_imrovedMatchEff_v2/bhnl_recoeff_dsadg.root'
  path_4 = 'LxyCompTwoDir_imrovedMatchEff_v2/bhnl_recoeff_displAll.root'

  f1 = ROOT.TFile.Open(path_1)
  f2 = ROOT.TFile.Open(path_2)
  f3 = ROOT.TFile.Open(path_3)
  f4 = ROOT.TFile.Open(path_4)
 
  hnum1 = f1.Get('hnum')
  hden1 = f1.Get('hden')
  hnum2 = f2.Get('hnum')
  hden2 = f2.Get('hden')
  hnum3 = f3.Get('hnum')
  hden3 = f3.Get('hden')
  hnum4 = f4.Get('hnum')
  hden4 = f4.Get('hden')

  if ROOT.TEfficiency.CheckConsistency(hnum1,hden1):
    eff1 = ROOT.TEfficiency(hnum1,hden1)
  if ROOT.TEfficiency.CheckConsistency(hnum2,hden2):
    eff2 = ROOT.TEfficiency(hnum2,hden2)
  if ROOT.TEfficiency.CheckConsistency(hnum3,hden3):
    eff3 = ROOT.TEfficiency(hnum3,hden3)
  if ROOT.TEfficiency.CheckConsistency(hnum4,hden4):
    eff4 = ROOT.TEfficiency(hnum4,hden4)

  eff1.SetLineColor(ROOT.kBlack)
  eff2.SetLineColor(ROOT.kBlue)
  eff3.SetLineColor(ROOT.kRed)
  eff4.SetLineColor(ROOT.kOrange)
  eff1.SetMarkerColor(ROOT.kBlack)
  eff2.SetMarkerColor(ROOT.kBlue)
  eff3.SetMarkerColor(ROOT.kRed)
  eff4.SetMarkerColor(ROOT.kOrange)


  eff1.SetTitle('(generalTracks) x (pat);gen L_{xy} (cm);Vtx reconstruction efficiency')
  eff2.SetTitle('(generalTracks) x (pat+dg);gen L_{xy} (cm);Vtx reconstruction efficiency')
  eff3.SetTitle('(generalTracks) x (pat+dg+dsa);gen L_{xy} (cm);Vtx reconstruction efficiency')
  eff4.SetTitle('(generalTracks+dispTracks) + (muons+dg+dsa);gen L_{xy} (cm);Vtx reconstruction efficiency')

  c = ROOT.TCanvas()
  eff1.Draw('AP')
  eff2.Draw('same')
  eff3.Draw('same')
  eff4.Draw('same')

  ROOT.gPad.Update(); 
  graph = eff1.GetPaintedGraph(); 
  graph.SetMinimum(0);
  graph.SetMaximum(0.5); 
  ROOT.gPad.Update(); 
  #eff2.GetXaxis().SetTitle('L_{xy} (cm)')
  #eff2.GetYaxis().SetTitle('Vtx reconstruction efficiency')

  c.BuildLegend(0.4,0.7,0.92,0.9)
  c.SaveAs('recoeff_comparison_Lxy.pdf')


  #### same exact thing for the pT
  hnum1 = f1.Get('hnum_pT')
  hden1 = f1.Get('hden_pT')
  hnum2 = f2.Get('hnum_pT')
  hden2 = f2.Get('hden_pT')
  hnum3 = f3.Get('hnum_pT')
  hden3 = f3.Get('hden_pT')
  hnum4 = f4.Get('hnum_pT')
  hden4 = f4.Get('hden_pT')

  if ROOT.TEfficiency.CheckConsistency(hnum1,hden1):
    eff1 = ROOT.TEfficiency(hnum1,hden1)
  if ROOT.TEfficiency.CheckConsistency(hnum2,hden2):
    eff2 = ROOT.TEfficiency(hnum2,hden2)
  if ROOT.TEfficiency.CheckConsistency(hnum3,hden3):
    eff3 = ROOT.TEfficiency(hnum3,hden3)
  if ROOT.TEfficiency.CheckConsistency(hnum4,hden4):
    eff4 = ROOT.TEfficiency(hnum4,hden4)

  eff1.SetLineColor(ROOT.kBlack)
  eff2.SetLineColor(ROOT.kBlue)
  eff3.SetLineColor(ROOT.kRed)
  eff4.SetLineColor(ROOT.kOrange)
  eff1.SetMarkerColor(ROOT.kBlack)
  eff2.SetMarkerColor(ROOT.kBlue)
  eff3.SetMarkerColor(ROOT.kRed)
  eff4.SetMarkerColor(ROOT.kOrange)


  eff1.SetTitle('(generalTracks) x (pat);gen HNL p_{T} (GeV);Vtx reconstruction efficiency')
  eff2.SetTitle('(generalTracks) x (pat+dg);gen HNL p_{T} (GeV);Vtx reconstruction efficiency')
  eff3.SetTitle('(generalTracks) x (pat+dg+dsa);gen HNL p_{T} (GeV);Vtx reconstruction efficiency')
  eff4.SetTitle('(generalTracks+dispTracks) + (muons+dg+dsa);gen HNL p_{T} (GeV);Vtx reconstruction efficiency')

  c = ROOT.TCanvas()
  eff1.Draw('AP')
  eff2.Draw('same')
  eff3.Draw('same')
  eff4.Draw('same')

  ROOT.gPad.Update(); 
  graph = eff1.GetPaintedGraph(); 
  graph.SetMinimum(0);
  graph.SetMaximum(0.5); 
  ROOT.gPad.Update(); 
  #eff2.GetXaxis().SetTitle('L_{xy} (cm)')
  #eff2.GetYaxis().SetTitle('Vtx reconstruction efficiency')

  c.BuildLegend(0.4,0.25,0.92,0.45)
  c.SaveAs('recoeff_comparison_pT.pdf')


  ### 2D plot

  ROOT.gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle2D.C')
  ROOT.gROOT.ProcessLine('setTDRStyle2D()')


  hnum1 = f1.Get('hnum_2D')
  hden1 = f1.Get('hden_2D')
  hnum2 = f2.Get('hnum_2D')
  hden2 = f2.Get('hden_2D')
  hnum3 = f3.Get('hnum_2D')
  hden3 = f3.Get('hden_2D')
  hnum4 = f4.Get('hnum_2D')
  hden4 = f4.Get('hden_2D')

  if ROOT.TEfficiency.CheckConsistency(hnum1,hden1):
    eff1 = ROOT.TEfficiency(hnum1,hden1)
  if ROOT.TEfficiency.CheckConsistency(hnum2,hden2):
    eff2 = ROOT.TEfficiency(hnum2,hden2)
  if ROOT.TEfficiency.CheckConsistency(hnum3,hden3):
    eff3 = ROOT.TEfficiency(hnum3,hden3)
  if ROOT.TEfficiency.CheckConsistency(hnum4,hden4):
    eff4 = ROOT.TEfficiency(hnum4,hden4)

  eff1.SetLineColor(ROOT.kBlack)
  eff2.SetLineColor(ROOT.kBlue)
  eff3.SetLineColor(ROOT.kRed)
  eff4.SetLineColor(ROOT.kOrange)
  eff1.SetMarkerColor(ROOT.kBlack)
  eff2.SetMarkerColor(ROOT.kBlue)
  eff3.SetMarkerColor(ROOT.kRed)
  eff4.SetMarkerColor(ROOT.kOrange)


  eff1.SetTitle('(generalTracks) x (pat);gen L_{xy} (cm);gen HNL p_{T} (GeV);Vtx reconstruction efficiency')
  eff2.SetTitle('(generalTracks) x (pat+dg);gen L_{xy} (cm);gen HNL p_{T} (GeV);Vtx reconstruction efficiency')
  eff3.SetTitle('(generalTracks) x (pat+dg+dsa);gen L_{xy} (cm);gen HNL p_{T} (GeV);Vtx reconstruction efficiency')
  eff4.SetTitle('(generalTracks+dispTracks) + (muons+dg+dsa);gen L_{xy} (cm);gen HNL p_{T} (GeV);Vtx reconstruction efficiency')

  c = ROOT.TCanvas()
  eff1.Draw('COLZ')
  eff1.Draw('TEXTsame')
  #ROOT.gPad.Update(); 
  #graph = eff1.GetPaintedGraph(); 
  #graph.SetMinimum(0);
  #graph.SetMaximum(0.5); 
  ROOT.gPad.Update(); 
  #eff2.GetXaxis().SetTitle('L_{xy} (cm)')
  #eff2.GetYaxis().SetTitle('Vtx reconstruction efficiency')
  c.SaveAs('recoeff_comparison_2D.pdf')

  c2 = ROOT.TCanvas()
  eff3.Draw('COLZ')
  eff3.Draw('TEXTsame')
  #ROOT.gPad.Update(); 
  #graph = eff3.GetPaintedGraph(); 
  #graph.SetMinimum(0);
  #graph.SetMaximum(0.5); 
  ROOT.gPad.Update(); 

  c2.SaveAs('recoeff_comparison_2D_dgdsa.pdf')

