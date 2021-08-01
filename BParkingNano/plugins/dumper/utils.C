#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

// ------------------------------------------ //
//  extra functions needed by the ntupliser
// ------------------------------------------ //


bool sortcansbydesc(const pair<int, float> &a1, const pair<int, float> &a2){
  return a1.second > a2.second;
}


bool sortcansbydesc_opp(const pair<int, float> &a1, const pair<int, float> &a2){
  return a1.second < a2.second;
}


float deltaR(float eta1, float eta2, float phi1, float phi2){
  float pi = 3.14159265359;
  float delta_eta = eta1 - eta2;
  float delta_phi = phi1 - phi2;
  if(fabs(delta_phi) > pi){
    delta_phi -= 2*pi;
  }
  return sqrt(pow(delta_eta, 2) + pow(delta_phi, 2));
}


bool checkLumi(string lumi_ranges, int lumi, int seed = 1){
  // get boundaries of lumi range
  string range_min = lumi_ranges.substr(lumi_ranges.find("[", seed)+1, int(lumi_ranges.find(",", seed)) - int(lumi_ranges.find("[", seed)+1));
  string range_max = lumi_ranges.substr(lumi_ranges.find(",", seed)+2, int(lumi_ranges.find("]", seed)) - int(lumi_ranges.find(",", seed)+2));
  //cout << range_min << " " << range_max << endl;

  if(lumi >= stoi(range_min) && lumi <= stoi(range_max)){
    //cout << "lumi accepted" << endl;
    return true;
  }
  else if(lumi_ranges.substr(int(lumi_ranges.find("]", seed))+1, 1) != "]"){ // there are further lumi blocks
    return checkLumi(lumi_ranges, lumi, int(lumi_ranges.find("]", seed)+2));
  }
  else{
    return false;
  }
}


bool lumiMask(int run, int lumi){
  //cout << "run " << run << " lumi " << lumi << endl;
  
  // put the content of the json in a tree
  boost::property_tree::ptree the_json_tree;
  boost::property_tree::read_json("golden_2018.json", the_json_tree);

  // checking that the run exists
  std::string lumi_ranges = the_json_tree.get<std::string> (std::to_string(run), "RUN NOT FOUND");

  if(lumi_ranges == "RUN NOT FOUND"){
    return false; 
  }
  else{
    // check that lumi is valid
    return checkLumi(lumi_ranges, lumi);
  }
}


float getTriggerScaleFactor(float pt, float eta){
  // get trigger scale factor file
  TString filename_sf = "/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/data/trigger_scale_factors/scaleFactor_results_cat_pt_eta_fit_A1.root";
  TFile* file_sf = TFile::Open(filename_sf);
  file_sf->cd();

  // get histogram
  TH2D* hist_sf = (TH2D*) file_sf->Get("hist_scale_factor")->Clone("hist_sf");

  pt = std::max(6., std::min(99.9, double(pt)));

  // get bin
  int bin_pt = hist_sf->GetXaxis()->FindBin(pt);
  int bin_eta = hist_sf->GetYaxis()->FindBin(eta);

  // get scale factor
  Float_t scale_factor = hist_sf->GetBinContent(bin_pt, bin_eta);
  
  file_sf->Close();

  return scale_factor;
}


