#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "PhysicsTools/BParkingNano/interface/helper.h"
#include <limits>
#include <algorithm>
#include "PhysicsTools/BParkingNano/interface/KinVtxFitter.h"
#include "PhysicsTools/BParkingNano/interface/ETHMuon.h"

class BToMuTrkPiBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
//  typedef std::vector<Lepton> LeptonCollection;
  typedef std::vector<pat::ETHMuon> MuonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BToMuTrkPiBuilder(const edm::ParameterSet &cfg):
    pi_selection_      {cfg.getParameter<std::string>("pionSelection"     )},
    isotrk_selection_  {cfg.getParameter<std::string>("isoTracksSelection")},
    trgmu_selection_   {cfg.getParameter<std::string>("trgMuonSelection"  )},
    lep_selection_     {cfg.getParameter<std::string>("leptonSelection"  )},
    pre_vtx_selection_ {cfg.getParameter<std::string>("preVtxSelection"   )},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection"  )},
    lepton_type_       {cfg.getParameter<std::string>("label")},
    isMC_              {cfg.getParameter<bool>("isMC")},

    // these two collections are ideally created beforehand by MuonTriggerSelector.cc
    //    * the former are muons that pass the preselection defined there AND match one of the 
    //      BParking triggers
    //    * the latter are all muons that pass the preselection (regardless whether they 
    //      fired the trigger). It's a superset of the previous collection
    trg_muons_         {consumes<MuonCollection>              ( cfg.getParameter<edm::InputTag>("trgMuons"               ) )},
    trks_              {consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("leptons"                ) )},
    trk_ttracks_       {consumes<TransientTrackCollection>         ( cfg.getParameter<edm::InputTag>("leptonsTransientTracks" ) )},
    pions_             {consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pions"                  ) )},
    pions_ttracks_     {consumes<TransientTrackCollection>         ( cfg.getParameter<edm::InputTag>("pionsTransientTracks"   ) )},
    isotracksToken_    {consumes<pat::PackedCandidateCollection>   ( cfg.getParameter<edm::InputTag>("tracks"                 ) )},
    isolostTracksToken_{consumes<pat::PackedCandidateCollection>   ( cfg.getParameter<edm::InputTag>("lostTracks"             ) )},
    genParticles_      {consumes<reco::GenParticleCollection>      ( cfg.getParameter<edm::InputTag>("genParticles"           ) )}, 
    beamspot_          {consumes<reco::BeamSpot>                   ( cfg.getParameter<edm::InputTag>("beamSpot"               ) )} 
    {
      produces<pat::CompositeCandidateCollection>();
    }

    // added for fetching the PV
    //vertexSrc_         { consumes<reco::VertexCollection>          ( iConfig.getParameter<edm::InputTag>( "vertexCollection"  ) )},
    //vertexSrc_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) )

  ~BToMuTrkPiBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  // pre-fitter preselection 
  const StringCutObjectSelector<pat::CompositeCandidate> pi_selection_; 
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; // not needed for the moment 
  const StringCutObjectSelector<pat::ETHMuon> trgmu_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> lep_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; 
  // post-fitter preselection 
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; 

  const std::string lepton_type_;
  const bool isMC_;

  const edm::EDGetTokenT<MuonCollection> trg_muons_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> trks_;
  const edm::EDGetTokenT<TransientTrackCollection> trk_ttracks_;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pions_;
  const edm::EDGetTokenT<TransientTrackCollection> pions_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;

  const edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  

  // for PV
  //const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  
};

void BToMuTrkPiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input

  edm::Handle<MuonCollection> trg_muons;
  evt.getByToken(trg_muons_, trg_muons);

/*  edm::Handle<LeptonCollection> leptons;
  evt.getByToken(leptons_, leptons);
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);
*/
  edm::Handle<pat::CompositeCandidateCollection> trks;
  evt.getByToken(trks_, trks);
  
  edm::Handle<TransientTrackCollection> trk_ttracks;
  evt.getByToken(trk_ttracks_, trk_ttracks);  

  edm::Handle<pat::CompositeCandidateCollection> pions;
  evt.getByToken(pions_, pions);
  
  edm::Handle<TransientTrackCollection> pions_ttracks;
  evt.getByToken(pions_ttracks_, pions_ttracks);  

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genParticles_, genParticles);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  


  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);

  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  // PV fetched for getting the trigger muon id (caveat: B is long lived)
  //edm::Handle<reco::VertexCollection> vertexHandle;
  //evt.getByToken(vertexSrc_, vertexHandle);
  //const reco::Vertex & PV = vertexHandle->front();

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  for(size_t trg_mu_idx = 0; trg_mu_idx < trg_muons->size(); ++trg_mu_idx) {

    size_t trg_mu_position = trks->size(); // make it point to just beyond the size of the collection
    
    edm::Ptr<pat::ETHMuon> trg_mu_ptr(trg_muons, trg_mu_idx);

    // only select the trigger muons (are slimmedMuons by construction)
    if(trg_mu_ptr->userInt("isTriggeringBPark") != 1) continue;	
   
  // selection on the trigger muon
    if( !trgmu_selection_(*trg_mu_ptr) ) continue;

      math::PtEtaPhiMLorentzVector trgMu_p4(
        trg_mu_ptr->pt(), 
        trg_mu_ptr->eta(),
        trg_mu_ptr->phi(),
        MUON_MASS
        );
    for(size_t trk_idx = 0; trk_idx < trks->size(); ++trk_idx) {
     
      edm::Ptr<pat::CompositeCandidate> trk_ptr(trks, trk_idx);
      // selection on the pion
      if( !lep_selection_(*trk_ptr) ) continue;
      
      math::PtEtaPhiMLorentzVector trk_p4(
        trk_ptr->pt(), 
        trk_ptr->eta(),
        trk_ptr->phi(),
        0.0005
        );
    for(size_t pi_idx = trk_idx+1; pi_idx < pions->size(); ++pi_idx) {
      edm::Ptr<pat::CompositeCandidate> pi_ptr(pions, pi_idx);
      if(trk_ptr->pt()==pi_ptr->pt() || (pi_idx==trk_idx))continue;

      // selection on the pion
      if( !pi_selection_(*pi_ptr) ) continue;
      
      math::PtEtaPhiMLorentzVector pi_p4(
        pi_ptr->pt(), 
        pi_ptr->eta(),
        pi_ptr->phi(),
        PI_MASS
        );
  
      // loop on selected muons and for a mu-pi candidate 
      // as well as a B candidate, that is HNL + trg mu
      
     /* for(size_t lep_idx = 0; lep_idx < leptons->size(); ++lep_idx) {}
        edm::Ptr<Lepton> lep_ptr(leptons, lep_idx);
        

        // the second muon must be _other_ than the trigger muon
        if(lep_ptr->pt()==trg_mu_ptr->pt()) { // lacking of any better idea for a comparison by pointer... 
            // save anyways the position in the collection
            // trigger muons are a subset of selected muons and selected muons are those that 
            // are saved in the tress eventually (see muonsBPark_cff.py), so
            // find the position of the trigger muon in the collection of selected muons
            trg_mu_position = lep_idx;
//             std::cout << __LINE__ << "]\t selected muon pt\t"     << sel_mu_ptr->pt()  << std::endl
//                                   << "    \t trigger  muon pt\t"  << trg_mu_ptr->pt()  << std::endl
//                                   << "    \t selected muon eta\t" << sel_mu_ptr->eta() << std::endl
//                                   << "    \t trigger  muon eta\t" << trg_mu_ptr->eta() << std::endl
//                                   << "    \t selected muon phi\t" << sel_mu_ptr->phi() << std::endl
//                                   << "    \t trigger  muon phi\t" << trg_mu_ptr->phi() << std::endl
//                                   << std::endl;
            continue;
        }else{
	
	trg_mu_position = trg_mu_idx;
	}
*/
        // selection on the lepton
  //      if( !lep_selection_(*lep_ptr) ) continue;

        // HNL candidate
	trg_mu_position = trg_mu_idx;
        pat::CompositeCandidate hnl_cand;
        hnl_cand.setP4(trk_ptr->p4() + pi_p4);
        hnl_cand.setCharge(trk_ptr->charge() + pi_ptr->charge());

        hnl_cand.addUserCand("trk", trk_ptr);
        hnl_cand.addUserFloat("trk_dxy", trk_ptr->userFloat("dxy"));
        hnl_cand.addUserCand("pi", pi_ptr);
        hnl_cand.addUserFloat("pi_dxy", pi_ptr->userFloat("dxy"));

        // check if pass pre vertex cut
        if( !pre_vtx_selection_(hnl_cand) ) continue;
        // fit the mu-pi vertex
        KinVtxFitter fitter(
          {pions_ttracks->at(trk_idx), pions_ttracks->at(pi_idx)},
          {trk_ptr->mass()             , PI_MASS                  },
          {LEP_SIGMA                   , PI_SIGMA                 } //some small sigma for the lepton mass
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
        b_cand.setP4(hnl_cand.p4() + trgMu_p4);
        b_cand.setCharge(hnl_cand.charge() + trg_mu_ptr->charge());
	//b_cand.addUserCand("hnl",hnl_cand);
//         b_cand.addUserCand("trg_mu", trg_mu_ptr);
        // https://cmssdt.cern.ch/lxr/source/DataFormats/Candidate/interface/Candidate.h
        
//         std::cout << __LINE__ << "]\t" << hnl_cand.pt() << std::endl;
//         std::cout << __LINE__ << "]\t" << hnl_cand.originalObjectRef().isNull() << std::endl;
//         edm::Ptr<pat::CompositeCandidate> hnl_cand_ptr = edm::refToPtr(hnl_cand.originalObjectRef());
//         b_cand.addUserCand("hnl", hnl_cand.originalObjectRef());
//         b_cand.addUserCand("hnl", hnl_cand.sourceCandidatePtr(0));



        b_cand.addDaughter(*trg_mu_ptr, "trg_mu");
        b_cand.addDaughter( hnl_cand  , "hnl"   );
        b_cand.addUserCand("sel_trk", trk_ptr);
        b_cand.addUserCand("pi", pi_ptr);

        b_cand.addUserInt  ("hnl_vtx_OK"             , fitter.success()                                                        );
        b_cand.addUserFloat("hnl_vtx_chi2"           , fitter.chi2()                                                           );
        b_cand.addUserFloat("hnl_vtx_ndof"           , fitter.dof()                                                            ); // float??
        b_cand.addUserFloat("hnl_vtx_prob"           , fitter.prob()                                                           );
        b_cand.addUserFloat("hnl_fitted_pt"          , fit_p4.pt()                                                             ); 
        b_cand.addUserFloat("hnl_fitted_eta"         , fit_p4.eta()                                                            );
        b_cand.addUserFloat("hnl_fitted_phi"         , fit_p4.phi()                                                            );
        b_cand.addUserFloat("hnl_fitted_mass"        , fitter.fitted_candidate().mass()                                        );      
        b_cand.addUserFloat("hnl_fitted_massErr"     , sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
        b_cand.addUserFloat("hnl_cos_theta_2D"       , cos_theta_2D(fitter, *beamspot, hnl_cand.p4())                          );
        b_cand.addUserFloat("hnl_fitted_cos_theta_2D", cos_theta_2D(fitter, *beamspot, fit_p4)                                 );
        b_cand.addUserFloat("hnl_fitted_cos_theta_2D_star", cos_theta_2D_star(fitter, *beamspot, fit_p4,trk_ptr->p4())                   );
        b_cand.addUserFloat("hnl_l_xy"               , lxy.value()                                                             );
        b_cand.addUserFloat("hnl_l_xy_unc"           , lxy.error()                                                             );
        b_cand.addUserFloat("hnl_ls_xy"              , lxy.value()/lxy.error()                                                 );
        b_cand.addUserFloat("hnl_charge"             , hnl_cand.charge()                                                       );
        b_cand.addUserFloat("hnl_vtx_x"              , hnl_cand.vx()                                                           );
        b_cand.addUserFloat("hnl_vtx_y"              , hnl_cand.vy()                                                           );
        b_cand.addUserFloat("hnl_vtx_z"              , hnl_cand.vz()                                                           );
        b_cand.addUserFloat("hnl_vtx_ex"             , sqrt(fitter.fitted_vtx_uncertainty().cxx())                             );
        b_cand.addUserFloat("hnl_vtx_ey"             , sqrt(fitter.fitted_vtx_uncertainty().cyy())                             );
        b_cand.addUserFloat("hnl_vtx_ez"             , sqrt(fitter.fitted_vtx_uncertainty().czz())                             );
        b_cand.addUserFloat("hnl_fitted_lep_pt"       , fitter.daughter_p4(0).pt()                                             ); 
        b_cand.addUserFloat("hnl_fitted_lep_eta"      , fitter.daughter_p4(0).eta()                                            );
        b_cand.addUserFloat("hnl_fitted_lep_phi"      , fitter.daughter_p4(0).phi()                                            );
        b_cand.addUserFloat("hnl_fitted_lep_mass"     , fitter.daughter_p4(0).mass()                                           );
        b_cand.addUserFloat("hnl_fitted_pi_pt"       , fitter.daughter_p4(1).pt()                                              ); 
        b_cand.addUserFloat("hnl_fitted_pi_eta"      , fitter.daughter_p4(1).eta()                                             );
        b_cand.addUserFloat("hnl_fitted_pi_phi"      , fitter.daughter_p4(1).phi()                                             );
        b_cand.addUserFloat("hnl_fitted_pi_mass"     , fitter.daughter_p4(1).mass()                                            );
      
        // adding trigger muon information to the b candidate
        b_cand.addUserFloat("trg_muon_pt"              , trg_mu_ptr->pt()                                                      );
        b_cand.addUserFloat("trg_muon_eta"             , trg_mu_ptr->eta()                                                     );
        b_cand.addUserFloat("trg_muon_phi"             , trg_mu_ptr->phi()                                                     );


        // difference between the z vertex position of the selected muon and tigger muon
        // computed at the prefit stage 
        b_cand.addUserFloat("dilepton_vzdiff"           , fabs(trg_mu_ptr->vz()-trk_ptr->vz())                               );
        b_cand.addUserFloat("dilepton_vxdiff"           , fabs(trg_mu_ptr->vx()-trk_ptr->vx())                               );
        b_cand.addUserFloat("dilepton_vydiff"           , fabs(trg_mu_ptr->vy()-trk_ptr->vy())                               );
        b_cand.addUserFloat("dilepton_Lxy"              , sqrt(pow(trg_mu_ptr->vx()-trk_ptr->vx(), 2) + pow(trg_mu_ptr->vy()-trk_ptr->vy(), 2)));
        b_cand.addUserFloat("dilepton_Lxyz"             , sqrt(pow(trg_mu_ptr->vx()-trk_ptr->vx(), 2) + pow(trg_mu_ptr->vy()-trk_ptr->vy(), 2) + pow(trg_mu_ptr->vz()-trk_ptr->vz(), 2)));
        
        // difference between the z vertex position of the pion and tigger muon
        b_cand.addUserFloat("pion_trgmuon_vzdiff"                , fabs(trg_mu_ptr->vz()-pi_ptr->vz())                            );
      
        // fetch the id of the sel muon at the secondary vertex. Not needed--->use instead info saved in the muonsBPark collection
      /*	 if(lepton_type_ == "muon"){
          float sel_muon_isSoft   = lep_ptr->isSoftMuon  ((const reco::Vertex&) fitter) ? 1. : 0. ;
          float sel_muon_isTight  = lep_ptr->isTightMuon ((const reco::Vertex&) fitter) ? 1. : 0. ;
          float sel_muon_isMedium = lep_ptr->isMediumMuon()                             ? 1. : 0. ;
          float sel_muon_isLoose  = lep_ptr->isLooseMuon ()                             ? 1. : 0. ;

          b_cand.addUserFloat("sel_muon_isSoft"       , sel_muon_isSoft                            );
          b_cand.addUserFloat("sel_muon_isTight"      , sel_muon_isTight                           );
          b_cand.addUserFloat("sel_muon_isMedium"     , sel_muon_isMedium                          );
          b_cand.addUserFloat("sel_muon_isLoose"      , sel_muon_isLoose                           );
        }*/

        
        // adding dR quantities (with fitted quantities)
        float dR_lep_pi = reco::deltaR(fitter.daughter_p4(0), fitter.daughter_p4(1)); 
        float dR_trgmu_hnl = reco::deltaR((*trg_mu_ptr), hnl_cand); 
        b_cand.addUserFloat("dr_lep_pi"              , dR_lep_pi                                 );
        b_cand.addUserFloat("dr_trgmu_hnl"          , dR_trgmu_hnl                               );


        // impact parameter variables (with pre-fit quantities)---> could be included in each particle collection but keeping them here for easy access
      /*  b_cand.addUserFloat("trg_muon_ip3d"   , fabs(trg_mu_ptr->dB(pat::Muon::PV3D))                                    );
        b_cand.addUserFloat("trg_muon_sip3d"  , fabs(trg_mu_ptr->dB(pat::Muon::PV3D) / trg_mu_ptr->edB(pat::Muon::PV3D)) );
        b_cand.addUserFloat("trg_muon_dxy"    , trg_mu_ptr->dB(pat::Muon::PV2D)                                          );
        b_cand.addUserFloat("trg_muon_dz"     , trg_mu_ptr->dB(pat::Muon::PVDZ)                                          );
       */ 
     //   if(lepton_type_ == "muon"){
       /*   b_cand.addUserFloat("sel_lep_ip3d"   , fabs(lep_ptr->dB(Lepton::PV3D))                                    );
          b_cand.addUserFloat("sel_lep_sip3d"  , fabs(lep_ptr->dB(Lepton::PV3D) / lep_ptr->edB(Lepton::PV3D)) );
          b_cand.addUserFloat("sel_lep_dxy"    , trk_ptr->dB(Lepton::PV2D)                                          );
          b_cand.addUserFloat("sel_lep_dz"     , trk_ptr->dB(Lepton::PVDZ)                                          );
       // }*/

        b_cand.addUserFloat("pion_dz"         , pi_ptr->userFloat("dz")                                                  );
        b_cand.addUserFloat("pion_dzS"        , pi_ptr->userFloat("dzS")                                                 );
        b_cand.addUserFloat("pion_dxy"        , pi_ptr->userFloat("dxy")                                                 );
        b_cand.addUserFloat("pion_dxyS"       , pi_ptr->userFloat("dxyS")                                                );
        b_cand.addUserFloat("pion_DCASig"     , pi_ptr->userFloat("DCASig")                                              );

        // isolation
        float trg_mu_iso03 = 0; 
        float trg_mu_iso04 = 0;
        float lep_iso03 = 0; 
        float lep_iso04 = 0;
        float pi_iso03  = 0; 
        float pi_iso04  = 0;
        float hnl_iso03 = 0;
        float hnl_iso04 = 0;

        // with conditions: best track + close to B 
        float lep_iso03_close = 0;
        float lep_iso04_close = 0;
        float trg_mu_iso03_close = 0; 
        float trg_mu_iso04_close = 0;
        float pi_iso03_close  = 0; 
        float pi_iso04_close  = 0;
        float hnl_iso03_close = 0;
        float hnl_iso04_close = 0;
        
        // nTracks     = iso_tracks->size();
        // totalTracks = nTracks + iso_lostTracks->size();

        for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
        
          const pat::PackedCandidate & track = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];

          // same preselection as for tracks
          if( !isotrk_selection_(track) ) continue;

          // check if the track is the pion
          if (pi_ptr->userCand("b_cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
         
          // check if the track is one of the two leptons 
          if (track_to_lepton_match(trg_mu_ptr, iso_tracks.id(), iTrk) || 
              track_to_lepton_match(trk_ptr, iso_tracks.id(), iTrk) ) continue;

          // add to final particle iso if dR < cone
          float dr_to_trgmu = deltaR(b_cand.userFloat("trg_muon_eta")     , b_cand.userFloat("trg_muon_phi")     , track.eta(), track.phi());
          float dr_to_lep = deltaR(b_cand.userFloat("hnl_fitted_lep_eta"), b_cand.userFloat("hnl_fitted_lep_phi"), track.eta(), track.phi());
          float dr_to_pi    = deltaR(b_cand.userFloat("hnl_fitted_pi_eta"), b_cand.userFloat("hnl_fitted_pi_phi"), track.eta(), track.phi());
          float dr_to_hnl   = deltaR(b_cand.userFloat("hnl_fitted_eta")   , b_cand.userFloat("hnl_fitted_phi")   , track.eta(), track.phi());

          if (dr_to_trgmu < 0.4){
            trg_mu_iso04 += track.pt();
            if ( dr_to_trgmu < 0.3) trg_mu_iso03 += track.pt();
          }
          if (dr_to_lep < 0.4){
            lep_iso04 += track.pt();
            if (dr_to_lep < 0.3)  lep_iso03 += track.pt();
          }
          if (dr_to_pi < 0.4){
            pi_iso04 += track.pt();
            if (dr_to_pi < 0.3) pi_iso03 += track.pt();
          }
          if (dr_to_hnl < 0.4){
            hnl_iso04 += track.pt();
            if (dr_to_hnl < 0.3) hnl_iso03 += track.pt();
          }

          // add requirement of the tracks to be close to the B
          if (!trg_mu_ptr->bestTrack() || fabs(track.dz() - trg_mu_ptr->bestTrack()->dz()) > 0.4) continue;
          if (!trk_ptr->bestTrack() || fabs(track.dz() - trk_ptr->bestTrack()->dz()) > 0.4) continue;
          //if (!pi_ptr->bestTrack() || fabs(track.dz() - pi_ptr->bestTrack()->dz()) > 0.4) continue; // pion never passes bestTrack requirement
          if (fabs(track.dz() - pi_ptr->userFloat("dz")) > 0.4) continue; //dropping requirement of best track

          if (dr_to_trgmu < 0.4){
            trg_mu_iso04_close += track.pt();
            if ( dr_to_trgmu < 0.3) trg_mu_iso03_close += track.pt();
          }
          if (dr_to_lep < 0.4){
            lep_iso04_close += track.pt();
            if (dr_to_lep < 0.3)  lep_iso03_close += track.pt();
          }
          if (dr_to_pi < 0.4){
            pi_iso04_close += track.pt();
            if (dr_to_pi < 0.3) pi_iso03_close += track.pt();
          }
          if (dr_to_hnl < 0.4){
            hnl_iso04_close += track.pt();
            if (dr_to_hnl < 0.3) hnl_iso03_close += track.pt();
          }
        }

        //trg_mu_iso03 /= trg_mu_ptr->pt();
        //trg_mu_iso04 /= trg_mu_ptr->pt();
        //sel_mu_iso03 /= fitter.daughter_p4(0).pt();
        //sel_mu_iso04 /= fitter.daughter_p4(0).pt();
        //pi_iso03 /= fitter.daughter_p4(1).pt();
        //pi_iso04 /= fitter.daughter_p4(1).pt();

        b_cand.addUserFloat("trg_mu_iso03", trg_mu_iso03);
        b_cand.addUserFloat("trg_mu_iso04", trg_mu_iso04);
        b_cand.addUserFloat("trk_iso03", lep_iso03);
        b_cand.addUserFloat("trk_iso04", lep_iso04);
        b_cand.addUserFloat("pi_iso03" , pi_iso03 );
        b_cand.addUserFloat("pi_iso04" , pi_iso04 );
        b_cand.addUserFloat("hnl_iso03" , hnl_iso03 );
        b_cand.addUserFloat("hnl_iso04" , hnl_iso04 );

        // add requirement of the tracks to be close to the B
        b_cand.addUserFloat("trg_mu_iso03_close", trg_mu_iso03_close);
        b_cand.addUserFloat("trg_mu_iso04_close", trg_mu_iso04_close);
        b_cand.addUserFloat("sel_trk_iso03_close", lep_iso03_close);
        b_cand.addUserFloat("sel_trk_iso04_close", lep_iso04_close);
        b_cand.addUserFloat("pi_iso03_close", pi_iso03_close);
        b_cand.addUserFloat("pi_iso04_close", pi_iso04_close);
        b_cand.addUserFloat("hnl_iso03_close", hnl_iso03_close);
        b_cand.addUserFloat("hnl_iso04_close", hnl_iso04_close);
        

        // position of the muons / tracks in their own collections
        b_cand.addUserInt("trg_mu_idx", trg_mu_position);
        b_cand.addUserInt("lep_idx", trk_idx);
        b_cand.addUserInt("pi_idx"    , pi_idx    );


        // di-lepton (for control channel)
        float dilepton_mass = (fitter.daughter_p4(0) + trgMu_p4).mass();
//	std::cout << "dilepton mass" <<  dilepton_mass<< std::endl;
        b_cand.addUserFloat("dilepton_mass", dilepton_mass);

        float dilepton_pt = (fitter.daughter_p4(0) + trgMu_p4).pt();
        b_cand.addUserFloat("dilepton_pt", dilepton_pt);

//	std::cout << "post vtx selection" <<  t_vtx_selection_(b_cand)<< std::endl;
//	std::cout << "lxy" << lxy.value()/lxy.error() << std::endl;
        // post fit selection
        if( !post_vtx_selection_(b_cand) ) continue;        


        // gen-matching
        
        int isMatched = 0;
        int trg_mu_genIdx(-1), sel_lep_genIdx(-1), pi_genIdx(-1);
        int genTriggerMuonMother_genPdgId(-1), genLeptonMother_genPdgId(-1), genPionMother_genPdgId(-1);

        // for MC only-- should work for both muon and electron final states now
        if(isMC_ == true){

          // pdgId of the gen particle to which the final-state particles are matched
          int trg_mu_genPdgId = trg_mu_ptr->userInt("mcMatch");
          int sel_lep_genPdgId = trk_ptr->userInt("mcMatch");
          int pi_genPdgId     = pi_ptr->userInt("mcMatch");
          
          // index of the gen particle to which the final-state particles are matched
          trg_mu_genIdx   = trg_mu_ptr->userInt("mcMatchIndex"); 
          sel_lep_genIdx   = trk_ptr->userInt("mcMatchIndex"); 
          pi_genIdx       = pi_ptr->userInt("mcMatchIndex"); 

          if(trg_mu_genIdx != -1 && sel_lep_genIdx != -1 && pi_genIdx != -1){

            // getting the associated gen particles
            edm::Ptr<reco::GenParticle> genTriggerMuon_ptr(genParticles, trg_mu_genIdx);
            edm::Ptr<reco::GenParticle> genLepton_ptr(genParticles, sel_lep_genIdx);
            edm::Ptr<reco::GenParticle> genPion_ptr(genParticles, pi_genIdx);

            // index of the associated mother particle
            int genTriggerMuonMother_genIdx = -1;
            int genLeptonMother_genIdx        = -1;
            int genPionMother_genIdx        = -1;
            if(genTriggerMuon_ptr->numberOfMothers()>0) genTriggerMuonMother_genIdx = genTriggerMuon_ptr->motherRef(0).key();
            if(genLepton_ptr->numberOfMothers()>0) genLeptonMother_genIdx = genLepton_ptr->motherRef(0).key();
            if(genPion_ptr->numberOfMothers()>0) genPionMother_genIdx = genPion_ptr->motherRef(0).key();

            // getting the mother particles
            edm::Ptr<reco::GenParticle> genTriggerMuonMother_ptr(genParticles, genTriggerMuonMother_genIdx);
            edm::Ptr<reco::GenParticle> genLeptonMother_ptr(genParticles, genLeptonMother_genIdx);
            edm::Ptr<reco::GenParticle> genPionMother_ptr(genParticles, genPionMother_genIdx);

            // pdgId of the mother particles
            genTriggerMuonMother_genPdgId = genTriggerMuonMother_ptr->pdgId();
            genLeptonMother_genPdgId        = genLeptonMother_ptr->pdgId();
            genPionMother_genPdgId        = genPionMother_ptr->pdgId();
            if(
               (fabs(sel_lep_genPdgId) == 211 || fabs(sel_lep_genPdgId) == 211) && fabs(genLeptonMother_genPdgId) == 9900015 && 
               fabs(pi_genPdgId) == 211 && fabs(genPionMother_genPdgId) == 9900015 &&
               fabs(trg_mu_genPdgId) == 13 && (fabs(genTriggerMuonMother_genPdgId) == 511 || fabs(genTriggerMuonMother_genPdgId) == 521 
                  || fabs(genTriggerMuonMother_genPdgId) == 531 || fabs(genTriggerMuonMother_genPdgId) == 541)
              ){
                isMatched = 1;
	    std::cout << "sel lep pdgId "<<sel_lep_genPdgId << " Lep mother pdgId " << genLeptonMother_genPdgId <<std::endl; 
	    std::cout << "sel pi pdgId "<<pi_genPdgId << " pi mother pdgId " << genPionMother_genPdgId <<std::endl; 
	    std::cout << "trg mu pdgId "<< trg_mu_genPdgId << " trgmu  mother pdgId " << genTriggerMuonMother_genPdgId <<std::endl; 
            }
          }
        }

        b_cand.addUserInt("isMatched", isMatched);
        b_cand.addUserInt("matching_trg_mu_genIdx", trg_mu_genIdx);
        b_cand.addUserInt("matching_sel_lep_genIdx", sel_lep_genIdx);
        b_cand.addUserInt("matching_pi_genIdx", pi_genIdx);
        b_cand.addUserInt("matching_trg_mu_motherPdgId", genTriggerMuonMother_genPdgId);
        b_cand.addUserInt("matching_sel_lep_motherPdgId", genLeptonMother_genPdgId);
        b_cand.addUserInt("matching_pi_motherPdgId", genPionMother_genPdgId);


        ret_val->push_back(b_cand);
              
      } // for(size_t sel_mu_idx = 0; sel_mu_idx < sel_muons->size(); ++sel_mu_idx)
      
    } // for(size_t pi_idx = 0; pi_idx < kaons->size(); ++pi_idx)

  } // for(size_t trg_mu_idx = 0; trg_mu_idx < trg_muons->size(); ++trg_mu_idx)
  
  evt.put(std::move(ret_val));
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToMuTrkPiBuilder);
/*
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//typedef BToMuLPiHighDisplacementBuilder<pat::Muon> BToMuMuPiHighDisplacementBuilder; //Builder for HNL in muon
typedef BToMuLPiHighDisplacementBuilder<pat::Electron> BToMuEPiHighDisplacementBuilder; //Builder for HNL in electron

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToMuEPiHighDisplacementBuilder);
*/
