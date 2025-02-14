#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 

#include "PhysicsTools/BParkingNano/interface/ETHMuon.h"
#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "PhysicsTools/BParkingNano/interface/helper.h"
#include <limits>
#include <algorithm>
#include "PhysicsTools/BParkingNano/interface/KinVtxFitter.h"

template<typename Lepton>
class BToMuLPiGeneralBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<Lepton> LeptonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BToMuLPiGeneralBuilder(const edm::ParameterSet &cfg):
    pi_selection_      {cfg.getParameter<std::string>("pionSelection"     )},
    isotrk_selection_  {cfg.getParameter<std::string>("isoTracksSelection")},
    trgmu_selection_   {cfg.getParameter<std::string>("trgMuonSelection"  )},
    selLep_selection_  {cfg.getParameter<std::string>("selLeptonSelection")},
    pre_vtx_selection_ {cfg.getParameter<std::string>("preVtxSelection"   )},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection"  )},
    isMC_              {cfg.getParameter<bool>("isMC")},

    // these two collections are ideally created beforehand by MuonTriggerSelector.cc
    //    * the former are muons that pass the preselection defined there AND match one of the 
    //      BParking triggers
    //    * the latter are all muons that pass the preselection (regardless whether they 
    //      fired the trigger). It's a superset of the previous collection
    trg_muons_         {consumes<pat::MuonCollection>              ( cfg.getParameter<edm::InputTag>("trgMuons"               ) )},
    sel_leptons_           {consumes<LeptonCollection>                 ( cfg.getParameter<edm::InputTag>("selLeptons"                ) )},
    sel_leptons_ttracks_ {consumes<TransientTrackCollection>         ( cfg.getParameter<edm::InputTag>("selLeptonsTransientTracks") )},
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

  ~BToMuLPiGeneralBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  // pre-fitter preselection 
  const StringCutObjectSelector<pat::CompositeCandidate> pi_selection_; 
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; // not needed for the moment 
  const StringCutObjectSelector<pat::Muon> trgmu_selection_; 
  const StringCutObjectSelector<Lepton> selLep_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; 
  // post-fitter preselection 
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; 

  const std::string lepton_type_;
  const bool isMC_;

  const edm::EDGetTokenT<pat::MuonCollection> trg_muons_;
  const edm::EDGetTokenT<LeptonCollection> sel_leptons_;
  const edm::EDGetTokenT<TransientTrackCollection> sel_leptons_ttracks_;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pions_;
  const edm::EDGetTokenT<TransientTrackCollection> pions_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;

  const edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  

  // for PV
  //const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  
};

template<typename Lepton>
void BToMuLPiGeneralBuilder<Lepton>::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {


  bool debug = false;
  //input

  edm::Handle<pat::MuonCollection> trg_muons;
  evt.getByToken(trg_muons_, trg_muons);

  edm::Handle<LeptonCollection> leptons;
  evt.getByToken(sel_leptons_, leptons);
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(sel_leptons_ttracks_, leptons_ttracks);

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

    size_t trg_mu_position = leptons->size(); // make it point to just beyond the size of the collection
    
    edm::Ptr<pat::Muon> trg_mu_ptr(trg_muons, trg_mu_idx);

    // selection on the trigger muon
    if( !trgmu_selection_(*trg_mu_ptr) ) continue;

    for(size_t pi_idx = 0; pi_idx < pions->size(); ++pi_idx) {
      edm::Ptr<pat::CompositeCandidate> pi_ptr(pions, pi_idx);

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
      for(size_t lep_idx = 0; lep_idx < leptons->size(); ++lep_idx) {
        edm::Ptr<Lepton> lep_ptr(leptons, lep_idx);
	if( debug) std::cout << "in the loop, before trg mu " << std::endl;
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
	
	if( debug) std::cout << "before lep selection" << std::endl;
	trg_mu_position = trg_mu_idx;
	}
	
	if( debug) std::cout << "before lep selection" << std::endl;
	if( debug) std::cout << "lepton pt " << lep_ptr->pt() << "lepton eta " << fabs(lep_ptr->eta())<< std::endl;
        // selection on the lepton
        //if( lep_ptr->pt()>0.5 && fabs(lep_ptr->eta())<2 ) continue;
        if( !selLep_selection_(*lep_ptr) ) continue;

	if( debug) std::cout << "after lep selection" << std::endl;
	int sel_lep_idx = lep_idx;
        // HNL candidate
        math::PtEtaPhiMLorentzVector lep_p4(
        lep_ptr->pt(), 
        lep_ptr->eta(),
        lep_ptr->phi(),
        lep_ptr->mass()
        );
        pat::CompositeCandidate hnl_cand;
        hnl_cand.setP4(lep_p4 + pi_p4);
        hnl_cand.setCharge(lep_ptr->charge() + pi_ptr->charge());

        hnl_cand.addUserCand("lep", lep_ptr);
        hnl_cand.addUserCand("pi", pi_ptr);

	if( debug) std::cout << "before pre vtx" << std::endl;
	if( debug) std::cout << "lepton pt " << lep_ptr->pt() << "lepton eta " << lep_ptr->eta()<< std::endl;
	if( debug) std::cout << "trgmu pt " << trg_mu_ptr->pt() << "lepton eta " << trg_mu_ptr->eta()<< std::endl;
	if( debug) std::cout << "pi pt " << pi_ptr->pt() << "lepton eta " << pi_ptr->eta()<< std::endl;
	if( debug) std::cout << "hnl pt " << hnl_cand.pt() << "hnl mass " << hnl_cand.mass()<< std::endl;
        // check if pass pre vertex cut
        if( !pre_vtx_selection_(hnl_cand) ) continue;
        // fit the mu-pi vertex
	if( debug) std::cout << "after pre vtx" << std::endl;
        KinVtxFitter fitter(
          {leptons_ttracks->at(lep_idx), pions_ttracks->at(pi_idx)},
          {lep_ptr->mass()             , PI_MASS                  },
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

	if( debug) std::cout << "after fit" << std::endl;
        auto fit_p4 = fitter.fitted_p4();
        auto lxy    = l_xy(fitter, *beamspot);

        // B candidate
        pat::CompositeCandidate b_cand;
        b_cand.setP4(hnl_cand.p4() + trg_mu_ptr->p4());
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
        b_cand.addUserCand("sel_lep", lep_ptr);
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
        b_cand.addUserFloat("hnl_fitted_lep_pt"       , fitter.daughter_p4(0).pt()                                              ); 
        b_cand.addUserFloat("hnl_fitted_lep_eta"      , fitter.daughter_p4(0).eta()                                             );
        b_cand.addUserFloat("hnl_fitted_lep_phi"      , fitter.daughter_p4(0).phi()                                             );
        b_cand.addUserFloat("hnl_fitted_lep_mass"     , fitter.daughter_p4(0).mass()                                            );
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
        b_cand.addUserFloat("dimuon_vzdiff"           , fabs(trg_mu_ptr->vz()-lep_ptr->vz())                               );
        b_cand.addUserFloat("dimuon_vxdiff"           , fabs(trg_mu_ptr->vx()-lep_ptr->vx())                               );
        b_cand.addUserFloat("dimuon_vydiff"           , fabs(trg_mu_ptr->vy()-lep_ptr->vy())                               );
        b_cand.addUserFloat("dimuon_Lxy"              , sqrt(pow(trg_mu_ptr->vx()-lep_ptr->vx(), 2) + pow(trg_mu_ptr->vy()-lep_ptr->vy(), 2)));
        b_cand.addUserFloat("dimuon_Lxyz"             , sqrt(pow(trg_mu_ptr->vx()-lep_ptr->vx(), 2) + pow(trg_mu_ptr->vy()-lep_ptr->vy(), 2) + pow(trg_mu_ptr->vz()-lep_ptr->vz(), 2)));
        
        // difference between the z vertex position of the pion and tigger muon
        b_cand.addUserFloat("pion_muon_vzdiff"                , fabs(trg_mu_ptr->vz()-pi_ptr->vz())                            );
      

        // fetch the id of the sel muon at the secondary vertex (use instead info saved in the muonsBPark collection?)
        
        // adding dR quantities (with fitted quantities)
        float dR_lep_pi = reco::deltaR(fitter.daughter_p4(0), fitter.daughter_p4(1)); 
        float dR_trgmu_hnl = reco::deltaR((*trg_mu_ptr), hnl_cand); 
        b_cand.addUserFloat("dr_lep_pi"              , dR_lep_pi                                   );
        b_cand.addUserFloat("dr_trgmu_hnl"          , dR_trgmu_hnl                               );


        // impact parameter variables (with pre-fit quantities)
        b_cand.addUserFloat("trg_muon_ip3d"   , fabs(trg_mu_ptr->dB(pat::Muon::PV3D))                                    );
        b_cand.addUserFloat("trg_muon_sip3d"  , fabs(trg_mu_ptr->dB(pat::Muon::PV3D) / trg_mu_ptr->edB(pat::Muon::PV3D)) );
        b_cand.addUserFloat("trg_muon_dxy"    , trg_mu_ptr->dB(pat::Muon::PV2D)                                          );
        b_cand.addUserFloat("trg_muon_dz"     , trg_mu_ptr->dB(pat::Muon::PVDZ)                                          );
        
        b_cand.addUserFloat("sel_lep_ip3d"   , fabs(lep_ptr->userFloat("ip3d"))                                          );
        b_cand.addUserFloat("sel_lep_sip3d"  , fabs(lep_ptr->userFloat("sip3d"))                                         );
        b_cand.addUserFloat("sel_lep_dxy"    , lep_ptr->userFloat("dxy")                                                 );
        b_cand.addUserFloat("sel_lep_dz"     , lep_ptr->userFloat("dz")                                                  );

        b_cand.addUserFloat("pion_dz"         , pi_ptr->userFloat("dz")                                                  );
        b_cand.addUserFloat("pion_dxy"        , pi_ptr->userFloat("dxy")                                                 );
        b_cand.addUserFloat("pion_dzS"        , pi_ptr->userFloat("dzS")                                                 );
        b_cand.addUserFloat("pion_dxyS"       , pi_ptr->userFloat("dxyS")                                                );
        b_cand.addUserFloat("pion_DCASig"     , pi_ptr->userFloat("DCASig")                                              );

        // isolation
        float trg_mu_iso03 = 0; 
        float trg_mu_iso04 = 0;
        float sel_lep_iso03 = 0; 
        float sel_lep_iso04 = 0;
        float pi_iso03  = 0; 
        float pi_iso04  = 0;
        float hnl_iso03 = 0;
        float hnl_iso04 = 0;
        // with conditions: best track + close to B 
        float sel_lep_iso03_close = 0; 
        float sel_lep_iso04_close = 0;
        float trg_mu_iso03_close = 0; 
        float trg_mu_iso04_close = 0;
        float pi_iso03_close  = 0; 
        float pi_iso04_close  = 0;
        float hnl_iso03_close = 0;
        float hnl_iso04_close = 0;
        
        // nTracks     = iso_tracks->size();
        // totalTracks = nTracks + iso_lostTracks->size();

        for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
        
          const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];

          // same preselection as for tracks
          if( !isotrk_selection_(trk) ) continue;

          // check if the track is the pion
          if (pi_ptr->userCand("b_cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
         
          // check if the track is one of the two leptons 
          if (track_to_lepton_match(trg_mu_ptr, iso_tracks.id(), iTrk) || 
              track_to_lepton_match(lep_ptr, iso_tracks.id(), iTrk) ) continue;

          // add to final particle iso if dR < cone
          float dr_to_trgmu = deltaR(b_cand.userFloat("trg_muon_eta")     , b_cand.userFloat("trg_muon_phi")     , trk.eta(), trk.phi());
          float dr_to_selLep= deltaR(b_cand.userFloat("hnl_fitted_lep_eta"), b_cand.userFloat("hnl_fitted_lep_phi"), trk.eta(), trk.phi());
          float dr_to_pi    = deltaR(b_cand.userFloat("hnl_fitted_pi_eta"), b_cand.userFloat("hnl_fitted_pi_phi"), trk.eta(), trk.phi());
          float dr_to_hnl   = deltaR(b_cand.userFloat("hnl_fitted_eta")   , b_cand.userFloat("hnl_fitted_phi")   , trk.eta(), trk.phi());

          if (dr_to_trgmu < 0.4){
            trg_mu_iso04 += trk.pt();
            if ( dr_to_trgmu < 0.3) trg_mu_iso03 += trk.pt();
          }
          if (dr_to_selLep < 0.4){
            sel_lep_iso04 += trk.pt();
            if (dr_to_selLep < 0.3)  sel_lep_iso03 += trk.pt();
          }
          if (dr_to_pi < 0.4){
            pi_iso04 += trk.pt();
            if (dr_to_pi < 0.3) pi_iso03 += trk.pt();
          }
          if (dr_to_hnl < 0.4){
            hnl_iso04 += trk.pt();
            if (dr_to_hnl < 0.3) hnl_iso03 += trk.pt();
          }

          // add requirement of the tracks to be close to the B
          if (!trg_mu_ptr->bestTrack() || fabs(trk.dz() - trg_mu_ptr->bestTrack()->dz()) > 0.4) continue;
          if (!lep_ptr->bestTrack() || fabs(trk.dz() - lep_ptr->bestTrack()->dz()) > 0.4) continue;
          //if (!pi_ptr->bestTrack() || fabs(trk.dz() - pi_ptr->bestTrack()->dz()) > 0.4) continue; // pion never passes bestTrack requirement
          if (fabs(trk.dz() - pi_ptr->userFloat("dz")) > 0.4) continue; //dropping requirement of best track

          if (dr_to_trgmu < 0.4){
            trg_mu_iso04_close += trk.pt();
            if ( dr_to_trgmu < 0.3) trg_mu_iso03_close += trk.pt();
          }
          if (dr_to_selLep < 0.4){
            sel_lep_iso04_close += trk.pt();
            if (dr_to_selLep < 0.3)  sel_lep_iso03_close += trk.pt();
          }
          if (dr_to_pi < 0.4){
            pi_iso04_close += trk.pt();
            if (dr_to_pi < 0.3) pi_iso03_close += trk.pt();
          }
          if (dr_to_hnl < 0.4){
            hnl_iso04_close += trk.pt();
            if (dr_to_hnl < 0.3) hnl_iso03_close += trk.pt();
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
        b_cand.addUserFloat("sel_lep_iso03", sel_lep_iso03);
        b_cand.addUserFloat("sel_lep_iso04", sel_lep_iso04);
        b_cand.addUserFloat("pi_iso03" , pi_iso03 );
        b_cand.addUserFloat("pi_iso04" , pi_iso04 );
        b_cand.addUserFloat("hnl_iso03" , hnl_iso03 );
        b_cand.addUserFloat("hnl_iso04" , hnl_iso04 );

        // add requirement of the tracks to be close to the B
        b_cand.addUserFloat("trg_mu_iso03_close", trg_mu_iso03_close);
        b_cand.addUserFloat("trg_mu_iso04_close", trg_mu_iso04_close);
        b_cand.addUserFloat("sel_lep_iso03_close", sel_lep_iso03_close);
        b_cand.addUserFloat("sel_lep_iso04_close", sel_lep_iso04_close);
        b_cand.addUserFloat("pi_iso03_close", pi_iso03_close);
        b_cand.addUserFloat("pi_iso04_close", pi_iso04_close);
        b_cand.addUserFloat("hnl_iso03_close", hnl_iso03_close);
        b_cand.addUserFloat("hnl_iso04_close", hnl_iso04_close);
        

        // position of the muons / tracks in their own collections
        b_cand.addUserInt("trg_mu_idx", trg_mu_position);
        b_cand.addUserInt("sel_lep_idx", sel_lep_idx);
        b_cand.addUserInt("pi_idx"    , pi_idx    );


        // di-lepton (for control channel)
        float dilepton_mass = (fitter.daughter_p4(0) + trg_mu_ptr->p4()).mass();
        b_cand.addUserFloat("dilepton_mass", dilepton_mass);

        float dilepton_pt = (fitter.daughter_p4(0) + trg_mu_ptr->p4()).pt();
        b_cand.addUserFloat("dilepton_pt", dilepton_pt);


        // post fit selection
        if( !post_vtx_selection_(b_cand) ) continue;        

        // gen-matching
        
        int isMatched = 0;
        int trg_mu_genIdx(-1), sel_lep_genIdx(-1), pi_genIdx(-1);
        int genTriggerMuonMother_genPdgId(-1), genLeptonMother_genPdgId(-1), genPionMother_genPdgId(-1);

        // for MC only-- should work for both muon and electron final states now
        if(isMC_ == true){
	std::cout << "in gen matching" << std::endl;
          // pdgId of the gen particle to which the final-state particles are matched
          int trg_mu_genPdgId = trg_mu_ptr->userInt("mcMatch");
          int sel_lep_genPdgId = lep_ptr->userInt("mcMatch");
          int pi_genPdgId     = pi_ptr->userInt("mcMatch");
          // index of the gen particle to which the final-state particles are matched
          trg_mu_genIdx   = trg_mu_ptr->userInt("mcMatchIndex"); 
          sel_lep_genIdx   = lep_ptr->userInt("mcMatchIndex"); 
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
		
	   std::cout << "sel lep pdgId "<<sel_lep_genPdgId << " Lep mother pdgId " << genLeptonMother_genPdgId <<std::endl; 
	    std::cout << "sel pi pdgId "<<pi_genPdgId << " pi mother pdgId " << genPionMother_genPdgId <<std::endl; 
	    std::cout << "trg mu pdgId "<< trg_mu_genPdgId << " trgmu  mother pdgId " << genTriggerMuonMother_genPdgId <<std::endl; 
            if(
               (fabs(sel_lep_genPdgId) == 13 || fabs(sel_lep_genPdgId) == 11) && fabs(genLeptonMother_genPdgId) == 9900015 && 
               fabs(pi_genPdgId) == 211 && fabs(genPionMother_genPdgId) == 9900015 &&
               fabs(trg_mu_genPdgId) == 13 && (fabs(genTriggerMuonMother_genPdgId) == 511 || fabs(genTriggerMuonMother_genPdgId) == 521 
                  || fabs(genTriggerMuonMother_genPdgId) == 531 || fabs(genTriggerMuonMother_genPdgId) == 541)
              ){
                isMatched = 1;
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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
typedef BToMuLPiGeneralBuilder<pat::ETHMuon> BToMuMuPiGeneralBuilder; //Builder for HNL in muon
typedef BToMuLPiGeneralBuilder<pat::Electron> BToMuEPiGeneralBuilder; //Builder for HNL in electron

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToMuMuPiGeneralBuilder);
DEFINE_FWK_MODULE(BToMuEPiGeneralBuilder);
