#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

class BToMuMuPiBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BToMuMuPiBuilder(const edm::ParameterSet &cfg):
    pi_selection_      {cfg.getParameter<std::string>("pionSelection"     )},
    pre_vtx_selection_ {cfg.getParameter<std::string>("preVtxSelection"   )},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection"  )},
    isotrk_selection_  {cfg.getParameter<std::string>("isoTracksSelection")},


    // these two collections are ideally created beforehand by MuonTriggerSelector.cc
    //    * the former are muons that pass the preselection defined there AND match one of the 
    //      BParking triggers
    //    * the latter are all muons that pass the preselection (regardless whether they 
    //      fired the trigger). It's a superset of the previous collection
    trg_muons_         {consumes<pat::MuonCollection>              ( cfg.getParameter<edm::InputTag>("trgMuons"               ) )},
    sel_muons_         {consumes<pat::MuonCollection>              ( cfg.getParameter<edm::InputTag>("selMuons"               ) )},
    sel_muons_ttracks_ {consumes<TransientTrackCollection>         ( cfg.getParameter<edm::InputTag>("selMuonsTransientTracks") )},
    
    pions_             {consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pions"                  ) )},
    pions_ttracks_     {consumes<TransientTrackCollection>         ( cfg.getParameter<edm::InputTag>("pionsTransientTracks"   ) )},
    isotracksToken_    {consumes<pat::PackedCandidateCollection>   ( cfg.getParameter<edm::InputTag>("tracks"                 ) )},
    isolostTracksToken_{consumes<pat::PackedCandidateCollection>   ( cfg.getParameter<edm::InputTag>("lostTracks"             ) )},
    beamspot_          {consumes<reco::BeamSpot>                   ( cfg.getParameter<edm::InputTag>("beamSpot"               ) )} 
    {
      produces<pat::CompositeCandidateCollection>();
    }

  ~BToMuMuPiBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> pi_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  const edm::EDGetTokenT<pat::MuonCollection> trg_muons_;
  const edm::EDGetTokenT<pat::MuonCollection> sel_muons_;
  const edm::EDGetTokenT<TransientTrackCollection> sel_muons_ttracks_;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pions_;
  const edm::EDGetTokenT<TransientTrackCollection> pions_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BToMuMuPiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input

  edm::Handle<pat::MuonCollection> trg_muons;
  evt.getByToken(trg_muons_, trg_muons);

  edm::Handle<pat::MuonCollection> sel_muons;
  evt.getByToken(sel_muons_, sel_muons);
  
  edm::Handle<TransientTrackCollection> sel_muons_ttracks;
  evt.getByToken(sel_muons_ttracks_, sel_muons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> pions;
  evt.getByToken(pions_, pions);
  
  edm::Handle<TransientTrackCollection> pions_ttracks;
  evt.getByToken(pions_ttracks_, pions_ttracks);  

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);
  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  std::vector<int> used_lep1_id, used_trk_id;


  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  for(size_t trg_mu_idx = 0; trg_mu_idx < trg_muons->size(); ++trg_mu_idx) {

    edm::Ptr<pat::Muon> trg_mu_ptr(trg_muons, trg_mu_idx);
    
    for(size_t pi_idx = 0; pi_idx < pions->size(); ++pi_idx) {
      edm::Ptr<pat::CompositeCandidate> pi_ptr(pions, pi_idx);
      if( !pi_selection_(*pi_ptr) ) continue;
      
      math::PtEtaPhiMLorentzVector k_p4(
        pi_ptr->pt(), 
        pi_ptr->eta(),
        pi_ptr->phi(),
        PI_MASS
        );
  
      // loop on selected muons and for a mu-pi candidate 
      // as well as a B candidate, that is HNL + trg mu
      for(size_t sel_mu_idx = 0; sel_mu_idx < sel_muons->size(); ++sel_mu_idx) {
        edm::Ptr<pat::Muon> sel_mu_ptr(sel_muons, sel_mu_idx);
        
        // the second muon must be _other_ than the trigger muon
        if(sel_mu_ptr==trg_mu_ptr) continue;

        // HNL candidate
        pat::CompositeCandidate hnl_cand;
        hnl_cand.setP4(sel_mu_ptr->p4() + k_p4);
        hnl_cand.setCharge(sel_mu_ptr->charge() + pi_ptr->charge());

        hnl_cand.addUserCand("mu", sel_mu_ptr);
        hnl_cand.addUserCand("pi", pi_ptr);

        // fit the mu-pi vertex
        KinVtxFitter fitter(
          {sel_muons_ttracks->at(sel_mu_idx), pions_ttracks->at(pi_idx)},
          {sel_mu_ptr->mass(), PI_MASS},
          {LEP_SIGMA, PI_SIGMA} //some small sigma for the lepton mass
        );
        if(!fitter.success()) continue; // hardcoded, but do we need otherwise?
        hnl_cand.setVertex( 
          reco::Candidate::Point( 
            fitter.fitted_vtx().x(),
            fitter.fitted_vtx().y(),
            fitter.fitted_vtx().z()
          )  
        );

        hnl_cand.addUserInt("sv_OK" , fitter.success());
        hnl_cand.addUserFloat("sv_chi2", fitter.chi2());
        hnl_cand.addUserFloat("sv_ndof", fitter.dof()); // float??
        hnl_cand.addUserFloat("sv_prob", fitter.prob());
        hnl_cand.addUserFloat("fitted_mll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass());
        auto fit_p4 = fitter.fitted_p4();
        hnl_cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
        hnl_cand.addUserFloat("fitted_eta" , fit_p4.eta());
        hnl_cand.addUserFloat("fitted_phi" , fit_p4.phi());
        hnl_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());      
        hnl_cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
        hnl_cand.addUserFloat(
          "cos_theta_2D", 
          cos_theta_2D(fitter, *beamspot, hnl_cand.p4())
          );
        hnl_cand.addUserFloat(
          "fitted_cos_theta_2D", 
          cos_theta_2D(fitter, *beamspot, fit_p4)
          );
        auto lxy = l_xy(fitter, *beamspot);
        hnl_cand.addUserFloat("l_xy", lxy.value());
        hnl_cand.addUserFloat("l_xy_unc", lxy.error());
        hnl_cand.addUserFloat("vtx_x", hnl_cand.vx());
        hnl_cand.addUserFloat("vtx_y", hnl_cand.vy());
        hnl_cand.addUserFloat("vtx_z", hnl_cand.vz());
        hnl_cand.addUserFloat("vtx_ex", sqrt(fitter.fitted_vtx_uncertainty().cxx()));
        hnl_cand.addUserFloat("vtx_ey", sqrt(fitter.fitted_vtx_uncertainty().cyy()));
        hnl_cand.addUserFloat("vtx_ez", sqrt(fitter.fitted_vtx_uncertainty().czz()));
  
        hnl_cand.addUserFloat("fitted_l1_pt" , fitter.daughter_p4(0).pt()); 
        hnl_cand.addUserFloat("fitted_l1_eta", fitter.daughter_p4(0).eta());
        hnl_cand.addUserFloat("fitted_l1_phi", fitter.daughter_p4(0).phi());
        hnl_cand.addUserFloat("fitted_l2_pt" , fitter.daughter_p4(1).pt()); 
        hnl_cand.addUserFloat("fitted_l2_eta", fitter.daughter_p4(1).eta());
        hnl_cand.addUserFloat("fitted_l2_phi", fitter.daughter_p4(1).phi());
        hnl_cand.addUserFloat("fitted_k_pt"  , fitter.daughter_p4(2).pt()); 
        hnl_cand.addUserFloat("fitted_k_eta" , fitter.daughter_p4(2).eta());
        hnl_cand.addUserFloat("fitted_k_phi" , fitter.daughter_p4(2).phi());

        // B candidate
        pat::CompositeCandidate b_cand;
        b_cand.setP4(hnl_cand.p4() + trg_mu_ptr->p4());
        b_cand.setCharge(hnl_cand.charge() + trg_mu_ptr->charge());

        b_cand.addUserCand("trg_mu", trg_mu_ptr);
// https://cmssdt.cern.ch/lxr/source/DataFormats/Candidate/interface/Candidate.h
//         sourceCandidatePtr()
        b_cand.addUserCand("hnl", hnl_cand.masterClonePtr());
              
      } // for(size_t sel_mu_idx = 0; sel_mu_idx < sel_muons->size(); ++sel_mu_idx)
      
    } // for(size_t pi_idx = 0; pi_idx < kaons->size(); ++pi_idx)

  } // for(size_t trg_mu_idx = 0; trg_mu_idx < trg_muons->size(); ++trg_mu_idx)
  
//   for (auto & b_cand: *ret_val){
//     b_cand.addUserInt("n_pi_used", std::count(used_trk_id.begin() , used_trk_id.end() , cand.userInt("pi_idx")));
//     b_cand.addUserInt("n_l1_used", std::count(used_lep1_id.begin(), used_lep1_id.end(), cand.userInt("l1_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l1_idx")));
//   }

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToMuMuPiBuilder);
