// usage: root -l PUReweighter.C+

#include<iostream>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;


void PUReweighter(){

  // get the mc pileup profile
  TFile* file_mc = TFile::Open("../../data/pileup/profiles/pileup_mc_2018_Autumn18.root");
  TH1D* hist_mc_fromfile = (TH1D*) file_mc->Get("pileup")->Clone("hist_mc_fromfile");
  hist_mc_fromfile->Scale(1./hist_mc_fromfile->Integral());

  // make sure the binning is the same as the one in data
  TH1D* hist_mc = new TH1D("hist_mc", "hist_mc", 200, 0, 200);
  for(int i(0); i<hist_mc_fromfile->GetNbinsX(); ++i){
    hist_mc->SetBinContent(i, hist_mc_fromfile->GetBinContent(i));
  }
  for(int i(hist_mc_fromfile->GetNbinsX()); i<201; ++i){
    hist_mc->SetBinContent(i, 0.);
  }

  // get the data pileup profile
  TFile* file_data = TFile::Open("../../data/pileup/profiles/pileup_data_2018.root");
  TH1D* hist_data = (TH1D*) file_data->Get("pileup_2018A")->Clone("hist_data");
  hist_data->Scale(1./hist_data->Integral());

  bool do_check = false;
  if(do_check){
    for(int i(0); i<hist_data->GetNbinsX(); ++i){
      float ratio = 0;
      if(hist_mc->GetBinContent(i) != 0.) ratio = hist_data->GetBinContent(i) / hist_mc->GetBinContent(i); 
      std::cout << i << " data " << hist_data->GetBinContent(i) << " mc " << hist_mc->GetBinContent(i) << " ratio " << ratio << std::endl;
    }
  }

  // get the pileup weight distribution
  TH1D* hist_weight = (TH1D*) hist_data->Clone("hist_weight");
  hist_weight->Divide(hist_mc);

  hist_weight->SetTitle("");
  hist_weight->SetMarkerSize(20);
  hist_weight->GetXaxis()->SetTitle("pileup data/MC");

  // store the weight in a root file
  TFile* root_file = TFile::Open("pileup_weight_dataA_mcAutumn18.root", "RECREATE");
  hist_weight->Write();
  root_file->Close();

  // plot some distributions
  TCanvas* c = new TCanvas("c", "c", 1200, 1000);
  TPad* pad = new TPad("pad", "pad", 0.01, 0.01, 0.99, 0.99);
  gStyle->SetPadLeftMargin(5);
  gStyle->SetPadRightMargin(5);
  pad->Draw();
  pad->cd();

  gStyle->SetOptStat(0);

  TH1D* hist_data_A = (TH1D*) file_data->Get("pileup_2018A")->Clone("hist_data_A");
  hist_data_A->Scale(1./hist_data_A->Integral());
  hist_data_A->SetMarkerStyle(41);
  hist_data_A->SetMarkerSize(2);
  hist_data_A->SetMarkerColor(kBlue+2);
  hist_data_A->SetFillColor(0);
  hist_data_A->SetLineWidth(0);

  TH1D* hist_data_B_5streams = (TH1D*) file_data->Get("pileup_2018B_5streams")->Clone("hist_data_B_5streams");
  hist_data_B_5streams->Scale(1./hist_data_B_5streams->Integral());
  hist_data_B_5streams->SetMarkerStyle(41);
  hist_data_B_5streams->SetMarkerSize(2);
  hist_data_B_5streams->SetMarkerColor(kBlue-4);
  hist_data_B_5streams->SetFillColor(0);
  hist_data_B_5streams->SetLineWidth(0);

  TH1D* hist_data_B_6streams = (TH1D*) file_data->Get("pileup_2018B_6streams")->Clone("hist_data_B_6streams");
  hist_data_B_6streams->Scale(1./hist_data_B_6streams->Integral());
  hist_data_B_6streams->SetMarkerStyle(41);
  hist_data_B_6streams->SetMarkerSize(2);
  hist_data_B_6streams->SetMarkerColor(kRed-10);
  hist_data_B_6streams->SetFillColor(0);
  hist_data_B_6streams->SetLineWidth(0);

  TH1D* hist_data_C = (TH1D*) file_data->Get("pileup_2018C")->Clone("hist_data_C");
  hist_data_C->Scale(1./hist_data_C->Integral());
  hist_data_C->SetMarkerStyle(41);
  hist_data_C->SetMarkerSize(2);
  hist_data_C->SetMarkerColor(kOrange+6);
  hist_data_C->SetFillColor(0);
  hist_data_C->SetLineWidth(0);

  TH1D* hist_data_D = (TH1D*) file_data->Get("pileup_2018D")->Clone("hist_data_D");
  hist_data_D->Scale(1./hist_data_D->Integral());
  hist_data_D->SetMarkerStyle(41);
  hist_data_D->SetMarkerSize(2);
  hist_data_D->SetMarkerColor(kOrange+8);
  hist_data_D->SetFillColor(0);
  hist_data_D->SetLineWidth(0);

  TH1D* hist_data_tot = (TH1D*) file_data->Get("pileup_2018total")->Clone("hist_data_tot");
  hist_data_tot->Scale(1./hist_data_tot->Integral());
  hist_data_tot->GetXaxis()->SetTitle("pile-up");
  hist_data_tot->GetXaxis()->SetRangeUser(0, 70);
  hist_data_tot->GetYaxis()->SetTitle("");
  hist_data_tot->GetYaxis()->SetRangeUser(0, 0.08);
  hist_data_tot->SetMarkerStyle(0);
  hist_data_tot->SetMarkerSize(0);
  hist_data_tot->SetFillColor(0);
  hist_data_tot->SetLineWidth(9);
  hist_data_tot->SetLineColor(1);

  hist_mc->SetMarkerStyle(0);
  hist_mc->SetMarkerSize(0);
  hist_mc->SetFillColor(0);
  hist_mc->SetLineWidth(9);
  hist_mc->SetLineColor(2);

  TH1D* hist_mc_weighted = (TH1D*) hist_mc->Clone("hist_mc_weighted");
  for(int i(0); i<hist_mc->GetNbinsX(); ++i){
    hist_mc_weighted->SetBinContent(i, hist_mc->GetBinContent(i)*hist_weight->GetBinContent(i));
  }
  hist_mc_weighted->SetFillColor(kOrange+4);
  hist_mc_weighted->SetFillStyle(3335);
  hist_mc_weighted->SetLineWidth(0);

  TLegend* leg = new TLegend(0.6, 0.6, 0.85, 0.85);
  leg->AddEntry(hist_data_A, "Data - A");
  leg->AddEntry(hist_data_B_5streams, "Data - B (5 streams)");
  leg->AddEntry(hist_data_B_6streams, "Data - B (6 streams)");
  leg->AddEntry(hist_data_C, "Data - C");
  leg->AddEntry(hist_data_D, "Data - D");
  leg->AddEntry(hist_data_tot, "Data - total");
  leg->AddEntry(hist_mc, "MC (Autumn 18)");
  leg->AddEntry(hist_mc_weighted, "MC reweighted");
  leg->SetTextSize(0.04);
  leg->SetLineColor(0);
  leg->SetFillColorAlpha(0, 0);
  leg->SetBorderSize(0);

  hist_data_tot->Draw("hist L");
  hist_data_A->Draw("same");
  hist_data_B_5streams->Draw("same");
  hist_data_B_6streams->Draw("same");
  hist_data_C->Draw("same");
  hist_data_D->Draw("same");
  hist_mc->Draw("hist L same");
  hist_mc_weighted->Draw("same");

  leg->Draw();

  c->SaveAs("pu_data_mc_reweighted.png");

}
