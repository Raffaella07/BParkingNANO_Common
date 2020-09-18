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
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Common/interface/RefToPtr.h"
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

    size_t trg_mu_position = sel_muons->size(); // make it point to just beyond the size of the collection
    
    edm::Ptr<pat::Muon> trg_mu_ptr(trg_muons, trg_mu_idx);
    
    for(size_t pi_idx = 0; pi_idx < pions->size(); ++pi_idx) {
      edm::Ptr<pat::CompositeCandidate> pi_ptr(pions, pi_idx);
      if( !pi_selection_(*pi_ptr) ) continue;
      
      math::PtEtaPhiMLorentzVector pi_p4(
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
        if(sel_mu_ptr->pt()==trg_mu_ptr->pt()) { // lacking of any better idea for a comparison by pointer... 
            // save anyways the position in the collection
            // trigger muons are a subset of selected muons and selected muons are those that 
            // are saved in the tress eventually (see muonsBPark_cff.py), so
            // find the position of the trigger muon in the collection of selected muons
            trg_mu_position = sel_mu_idx;
//             std::cout << __LINE__ << "]\t selected muon pt\t"     << sel_mu_ptr->pt()  << std::endl
//                                   << "    \t trigger  muon pt\t"  << trg_mu_ptr->pt()  << std::endl
//                                   << "    \t selected muon eta\t" << sel_mu_ptr->eta() << std::endl
//                                   << "    \t trigger  muon eta\t" << trg_mu_ptr->eta() << std::endl
//                                   << "    \t selected muon phi\t" << sel_mu_ptr->phi() << std::endl
//                                   << "    \t trigger  muon phi\t" << trg_mu_ptr->phi() << std::endl
//                                   << std::endl;
            continue;
        }

        // HNL candidate
        pat::CompositeCandidate hnl_cand;
        hnl_cand.setP4(sel_mu_ptr->p4() + pi_p4);
        hnl_cand.setCharge(sel_mu_ptr->charge() + pi_ptr->charge());

        hnl_cand.addUserCand("mu", sel_mu_ptr);
        hnl_cand.addUserCand("pi", pi_ptr);

        // check if pass pre vertex cut
        if( !pre_vtx_selection_(hnl_cand) ) continue;

        // fit the mu-pi vertex
        KinVtxFitter fitter(
          {sel_muons_ttracks->at(sel_mu_idx), pions_ttracks->at(pi_idx)},
          {sel_mu_ptr->mass()               , PI_MASS                  },
          {LEP_SIGMA                        , PI_SIGMA                 } //some small sigma for the lepton mass
        );
        if(!fitter.success()) continue; // hardcoded, but do we need otherwise?
        hnl_cand.setVertex( 
          reco::Candidate::Point( 
            fitter.fitted_vtx().x(),
            fitter.fitted_vtx().y(),
            fitter.fitted_vtx().z()
          )  
        );

        auto fit_p4 = fitter.fitted_p4();
        auto lxy    = l_xy(fitter, *beamspot);

        // B candidate
        pat::CompositeCandidate b_cand;
        b_cand.setP4(hnl_cand.p4() + trg_mu_ptr->p4());
        b_cand.setCharge(hnl_cand.charge() + trg_mu_ptr->charge());

//         b_cand.addUserCand("trg_mu", trg_mu_ptr);
        // https://cmssdt.cern.ch/lxr/source/DataFormats/Candidate/interface/Candidate.h
        
//         std::cout << __LINE__ << "]\t" << hnl_cand.pt() << std::endl;
//         std::cout << __LINE__ << "]\t" << hnl_cand.originalObjectRef().isNull() << std::endl;
//         edm::Ptr<pat::CompositeCandidate> hnl_cand_ptr = edm::refToPtr(hnl_cand.originalObjectRef());
//         b_cand.addUserCand("hnl", hnl_cand.originalObjectRef());
//         b_cand.addUserCand("hnl", hnl_cand.sourceCandidatePtr(0));


        b_cand.addDaughter(*trg_mu_ptr, "trg_mu");
        b_cand.addDaughter( hnl_cand  , "hnl"   );

        b_cand.addUserInt  ("hnlVtxOK"             , fitter.success()                                                        );
        b_cand.addUserFloat("hnlVtxChi2"           , fitter.chi2()                                                           );
        b_cand.addUserFloat("hnlVtxNdof"           , fitter.dof()                                                            ); // float??
        b_cand.addUserFloat("hnlVtxProb"           , fitter.prob()                                                           );
        b_cand.addUserFloat("hnlFittedPt"          , fit_p4.pt()                                                             ); 
        b_cand.addUserFloat("hnlFittedEta"         , fit_p4.eta()                                                            );
        b_cand.addUserFloat("hnlFittedPhi"         , fit_p4.phi()                                                            );
        b_cand.addUserFloat("hnlFittedMass"        , fitter.fitted_candidate().mass()                                        );      
        b_cand.addUserFloat("hnlFittedMassErr"     , sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
        b_cand.addUserFloat("hnlCosTheta2D"        , cos_theta_2D(fitter, *beamspot, hnl_cand.p4())                          );
        b_cand.addUserFloat("hnlFittedCosTheta2D"  , cos_theta_2D(fitter, *beamspot, fit_p4)                                 );
        b_cand.addUserFloat("hnlLxy"               , lxy.value()                                                             );
        b_cand.addUserFloat("hnlLxyUnc"            , lxy.error()                                                             );
        b_cand.addUserFloat("hnlLlsxy"             , lxy.value()/lxy.error()                                                 );
        b_cand.addUserFloat("hnlVtxX"              , hnl_cand.vx()                                                           );
        b_cand.addUserFloat("hnlVtxY"              , hnl_cand.vy()                                                           );
        b_cand.addUserFloat("hnlVtxZ"              , hnl_cand.vz()                                                           );
        b_cand.addUserFloat("hnlVtxXE"             , sqrt(fitter.fitted_vtx_uncertainty().cxx())                             );
        b_cand.addUserFloat("hnlVtxYE"             , sqrt(fitter.fitted_vtx_uncertainty().cyy())                             );
        b_cand.addUserFloat("hnlVtxZE"             , sqrt(fitter.fitted_vtx_uncertainty().czz())                             );
        b_cand.addUserFloat("hnlFittedMuPt"        , fitter.daughter_p4(0).pt()                                              ); 
        b_cand.addUserFloat("hnlFittedMuEta"       , fitter.daughter_p4(0).eta()                                             );
        b_cand.addUserFloat("hnlFittedMuPhi"       , fitter.daughter_p4(0).phi()                                             );
        b_cand.addUserFloat("hnlFittedPiPt"        , fitter.daughter_p4(1).pt()                                              ); 
        b_cand.addUserFloat("hnlFittedPiEta"       , fitter.daughter_p4(1).eta()                                             );
        b_cand.addUserFloat("hnlFittedPiPhi"       , fitter.daughter_p4(1).phi()                                             );
        
        // position of the muons / tracks in their own collections
        b_cand.addUserInt("triggerMuonIdx", trg_mu_position);
        b_cand.addUserInt("muonIdx", sel_mu_idx);
        b_cand.addUserInt("probeTracksIdx", pi_idx    );

        // post fit selection
        if( !post_vtx_selection_(b_cand) ) continue;        

        ret_val->push_back(b_cand);
              
      } // for(size_t sel_mu_idx = 0; sel_mu_idx < sel_muons->size(); ++sel_mu_idx)
      
    } // for(size_t pi_idx = 0; pi_idx < kaons->size(); ++pi_idx)

  } // for(size_t trg_mu_idx = 0; trg_mu_idx < trg_muons->size(); ++trg_mu_idx)
  
  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToMuMuPiBuilder);
