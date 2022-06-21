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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "../interface/helper.h"
#include <limits>
#include <algorithm>
#include "PhysicsTools/BParkingNano/interface/KinVtxFitter.h"
#include "PhysicsTools/BParkingNano/interface/ETHMuon.h"

class HNLToMuPiBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
  public:
    typedef pat::ETHMuon Lepton;
    typedef std::vector<Lepton> LeptonCollection;
    typedef std::vector<pat::ETHMuon> ETHMuonCollection;
    typedef std::vector<reco::TransientTrack> TransientTrackCollection;

    explicit HNLToMuPiBuilder(const edm::ParameterSet &cfg):
      pi_selection_      {cfg.getParameter<std::string>("pionSelection"     )},
      lep_selection_     {cfg.getParameter<std::string>("leptonSelection"  )},
      pre_vtx_selection_ {cfg.getParameter<std::string>("preVtxSelection"   )},
      post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection"  )},
      lepton_type_       {cfg.getParameter<std::string>("label")},
      isMC_              {cfg.getParameter<bool>("isMC")},

      leptons_           {consumes<LeptonCollection>                 ( cfg.getParameter<edm::InputTag>("leptons"                ) )},
      leptons_ttracks_   {consumes<TransientTrackCollection>         ( cfg.getParameter<edm::InputTag>("leptonsTransientTracks" ) )},
      pions_             {consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pions"                  ) )},
      pions_ttracks_     {consumes<TransientTrackCollection>         ( cfg.getParameter<edm::InputTag>("pionsTransientTracks"   ) )},
      genParticles_      {consumes<reco::GenParticleCollection>      ( cfg.getParameter<edm::InputTag>("genParticles"           ) )}, 
      beamspot_          {consumes<reco::BeamSpot>                   ( cfg.getParameter<edm::InputTag>("beamSpot"               ) )} 
    {
      produces<pat::CompositeCandidateCollection>();
    }

    ~HNLToMuPiBuilder() override {}

    void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

  private:
    // pre-fitter preselection 
    const StringCutObjectSelector<pat::CompositeCandidate> pi_selection_; 
    const StringCutObjectSelector<Lepton> lep_selection_; 
    const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; 
    // post-fitter preselection 
    const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; 

    const std::string lepton_type_;
    const bool isMC_;

    const edm::EDGetTokenT<LeptonCollection> leptons_;
    const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;

    const edm::EDGetTokenT<pat::CompositeCandidateCollection> pions_;
    const edm::EDGetTokenT<TransientTrackCollection> pions_ttracks_;

    const edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;

    const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void HNLToMuPiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<LeptonCollection> leptons;
  evt.getByToken(leptons_, leptons);

  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> pions;
  evt.getByToken(pions_, pions);

  edm::Handle<TransientTrackCollection> pions_ttracks;
  evt.getByToken(pions_ttracks_, pions_ttracks);  

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genParticles_, genParticles);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  //std::cout << std::endl; 

  int gen_pdgId(0);
  float gen_hnl_pt(0.);
  float gen_hnl_eta(0.);
  float gen_hnl_mass(0.);
  float gen_mu_pt(0.);
  float gen_mu_eta(0.);
  float gen_pi_pt(0.);
  float gen_pi_eta(0.);

  /*
  if(isMC_ == true){

    unsigned int hnl_idx(-99);

    for(size_t gen_idx = 0; gen_idx < genParticles->size(); ++gen_idx) {
      edm::Ptr<reco::GenParticle> genPart_ptr(genParticles, gen_idx);
      if(fabs(genPart_ptr->pdgId()) == 9900015){
        //std::cout << "found a hnl pt " << genPart_ptr->pt() << " eta " << genPart_ptr->eta() << std::endl;
        hnl_idx = gen_idx;
        gen_hnl_pt = genPart_ptr->pt();
        gen_hnl_eta = genPart_ptr->eta();
        gen_hnl_mass = genPart_ptr->mass();
        //if(genPart_ptr->numberOfMothers()>0){
        //  mother_idx = genPart_ptr->motherRef(0).key();
        //}
        break;
      }
    }
    int n_muons(0);
    for(size_t gen_idx = 0; gen_idx < genParticles->size(); ++gen_idx) {
      edm::Ptr<reco::GenParticle> genPart_ptr(genParticles, gen_idx);
      //if(genPart_ptr->numberOfMothers()>0 && genPart_ptr->motherRef(0).key() == mother_idx){
        //std::cout << "found a mother " << genPart_ptr->pdgId() << std::endl;
      //}
      //if(genPart_ptr->numberOfMothers()>0) std::cout << genPart_ptr->motherRef(0).key() << " " << hnl_idx << std::endl;
      if(genPart_ptr->numberOfMothers()>0 && genPart_ptr->motherRef(0).key() == hnl_idx){
        //std::cout << "found a daughter" << " pdgId: " << genPart_ptr->pdgId() << " pt " << genPart_ptr->pt() << " eta " << genPart_ptr->eta() << std::endl;
        //if(abs(genPart_ptr->pdgId())==11){
        //  gen_pdgId = 11;
        //}
        if(abs(genPart_ptr->pdgId())==13){
          //std::cout << "found a muon daughter " << genPart_ptr->pt() << " eta " << genPart_ptr->eta() << std::endl;
          gen_mu_pt = genPart_ptr->pt();
          gen_mu_eta = genPart_ptr->eta();
          ++n_muons;
        }
        if(abs(genPart_ptr->pdgId())==211){
          //std::cout << "found a pion daughter " << genPart_ptr->pt() << " eta " << genPart_ptr->eta() << std::endl;
          gen_pi_pt = genPart_ptr->pt();
          gen_pi_eta = genPart_ptr->eta();
        }
      }
    }
  }
  */


  //std::cout << "HNL builder" << std::endl;
  //std::cout << "pion size " << pions->size() << " muon size " << leptons->size() << std::endl;


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

    //if(pi_ptr->pt()<9.5){
    //  continue;
    //}
    //else std::cout << "reco pion pt " << pi_ptr->pt() << " eta " << pi_ptr->eta() << std::endl;

    // loop on selected muons and for a mu-pi candidate 
    for(size_t lep_idx = 0; lep_idx < leptons->size(); ++lep_idx) {
      edm::Ptr<Lepton> lep_ptr(leptons, lep_idx);
      //if(lep_ptr->isDSAMuon()) continue;

      // selection on the lepton
      if( !lep_selection_(*lep_ptr) ) continue;

      math::PtEtaPhiMLorentzVector lep_p4(
        lep_ptr->pt(), 
        lep_ptr->eta(),
        lep_ptr->phi(),
        lep_ptr->mass()
        );

      //if(lep_ptr->pt()<3.3){
      //  continue;
      //}
      //else std::cout << "reco muon pt " << lep_ptr->pt() << " eta " << lep_ptr->eta() << " isDSA " << lep_ptr->isDSAMuon() << std::endl;


      //std::cout << std::endl << std::endl << std::endl;

      // HNL candidate
      pat::CompositeCandidate hnl_cand;
      hnl_cand.setP4(lep_p4 + pi_p4);
      //hnl_cand.setP4(lep_ptr->P4() + pi_p4);
      hnl_cand.setCharge(lep_ptr->charge() + pi_ptr->charge());

      hnl_cand.addUserCand("lep", lep_ptr);
      hnl_cand.addUserCand("pi", pi_ptr);

      // check if pass pre vertex cut
      if( !pre_vtx_selection_(hnl_cand) ) continue;
      //std::cout << "passed prefitter selection" << std::endl;

      //std::cout << "pt " << lep_ptr->pt() - leptons_ttracks->at(lep_idx).track().pt() << " isDSA " << lep_ptr->isDSAMuon() << std::endl;

      //std::cout << "muon pt " << leptons_ttracks->at(lep_idx).track().pt() << " eta " << leptons_ttracks->at(lep_idx).track().eta() << " vx " << leptons_ttracks->at(lep_idx).track().vx() << " vz " << leptons_ttracks->at(lep_idx).track().vz() << std::endl;
      //std::cout << "pion pt " << pions_ttracks->at(lep_idx).track().pt() << " eta " << pions_ttracks->at(lep_idx).track().eta() << " vx " << pions_ttracks->at(lep_idx).track().vx() << " vz " << pions_ttracks->at(lep_idx).track().vz() << std::endl;

      // fit the mu-pi vertex
      KinVtxFitter fitter(
          {leptons_ttracks->at(lep_idx), pions_ttracks->at(pi_idx)},
          {lep_ptr->mass()             , PI_MASS                  },
          {LEP_SIGMA                   , PI_SIGMA                 } //some small sigma for the lepton mass
          );
      // fit with mass constraint
      //KinVtxFitter fitter(
      //    {leptons_ttracks->at(lep_idx), pions_ttracks->at(pi_idx)},
      //    {lep_ptr->mass()             , PI_MASS                  },
      //    {LEP_SIGMA                   , PI_SIGMA                 }, //some small sigma for the lepton mass
      //    3.0
      //    );
      if(!fitter.success()){
        //std::cout << "no fit " << std::endl;
        continue;
      }
      //else std::cout << "fitted" << std::endl;
      hnl_cand.setVertex( 
          reco::Candidate::Point( 
            fitter.fitted_vtx().x(),
            fitter.fitted_vtx().y(),
            fitter.fitted_vtx().z()
            )  
          );

      auto fit_p4 = fitter.fitted_p4();
      auto lxy    = l_xy(fitter, *beamspot);

      hnl_cand.addUserInt  ("hnl_vtx_OK"             , fitter.success()                                                        );
      hnl_cand.addUserFloat("hnl_vtx_chi2"           , fitter.chi2()                                                           );
      hnl_cand.addUserFloat("hnl_vtx_ndof"           , fitter.dof()                                                            ); // float??
      hnl_cand.addUserFloat("hnl_vtx_prob"           , fitter.prob()                                                           );
      hnl_cand.addUserFloat("hnl_fitted_pt"          , fit_p4.pt()                                                             ); 
      hnl_cand.addUserFloat("hnl_fitted_eta"         , fit_p4.eta()                                                            );
      hnl_cand.addUserFloat("hnl_fitted_phi"         , fit_p4.phi()                                                            );
      hnl_cand.addUserFloat("hnl_fitted_mass"        , fitter.fitted_candidate().mass()                                        );      
      hnl_cand.addUserFloat("hnl_fitted_massErr"     , sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
      hnl_cand.addUserFloat("hnl_cos_theta_2D"       , cos_theta_2D(fitter, *beamspot, hnl_cand.p4())                          );
      hnl_cand.addUserFloat("hnl_fitted_cos_theta_2D", cos_theta_2D(fitter, *beamspot, fit_p4)                                 );
      hnl_cand.addUserFloat("hnl_l_xy"               , lxy.value()                                                             );
      hnl_cand.addUserFloat("hnl_l_xy_unc"           , lxy.error()                                                             );
      hnl_cand.addUserFloat("hnl_ls_xy"              , lxy.value()/lxy.error()                                                 );
      hnl_cand.addUserInt  ("hnl_charge"             , hnl_cand.charge()                                                       );
      hnl_cand.addUserFloat("hnl_vtx_x"              , hnl_cand.vx()                                                           );
      hnl_cand.addUserFloat("hnl_vtx_y"              , hnl_cand.vy()                                                           );
      hnl_cand.addUserFloat("hnl_vtx_z"              , hnl_cand.vz()                                                           );
      hnl_cand.addUserFloat("hnl_vtx_ex"             , sqrt(fitter.fitted_vtx_uncertainty().cxx())                             );
      hnl_cand.addUserFloat("hnl_vtx_ey"             , sqrt(fitter.fitted_vtx_uncertainty().cyy())                             );
      hnl_cand.addUserFloat("hnl_vtx_ez"             , sqrt(fitter.fitted_vtx_uncertainty().czz())                             );
      hnl_cand.addUserFloat("hnl_fitted_lep_pt"      , fitter.daughter_p4(0).pt()                                              ); 
      hnl_cand.addUserFloat("hnl_fitted_lep_eta"     , fitter.daughter_p4(0).eta()                                             );
      hnl_cand.addUserFloat("hnl_fitted_lep_phi"     , fitter.daughter_p4(0).phi()                                             );
      hnl_cand.addUserFloat("hnl_fitted_lep_mass"    , fitter.daughter_p4(0).mass()                                            );
      hnl_cand.addUserFloat("hnl_fitted_pi_pt"       , fitter.daughter_p4(1).pt()                                              ); 
      hnl_cand.addUserFloat("hnl_fitted_pi_eta"      , fitter.daughter_p4(1).eta()                                             );
      hnl_cand.addUserFloat("hnl_fitted_pi_phi"      , fitter.daughter_p4(1).phi()                                             );
      hnl_cand.addUserFloat("hnl_fitted_pi_mass"     , fitter.daughter_p4(1).mass()                                            );

      //std::cout << "mass 1 " << fitter.fitted_candidate().mass() << " mass 2 " << (fitter.daughter_p4(0)+daughter_p4(1)).mass() << " mass 3 " << fitter.fitted_p4().mass() << std::endl;
      //std::cout << "mass 1 " << fitter.fitted_candidate().mass() << " mass 2 " << fitter.fitted_p4().mass() << std::endl;

      // post fit selection
      if( !post_vtx_selection_(hnl_cand) ) continue;        

      // position of the muons / tracks in their own collections
      hnl_cand.addUserInt("lep_idx", lep_idx);
      hnl_cand.addUserInt("pi_idx", pi_idx);

      // invariant mass (before fitting)
      float reco_hnl_mass = (lep_ptr->P4() + pi_ptr->p4()).mass();
      hnl_cand.addUserFloat("hnl_mass_unfitted", reco_hnl_mass);
      hnl_cand.addUserFloat("diff_mass", fabs(reco_hnl_mass-fitter.fitted_candidate().mass())/reco_hnl_mass);

      // gen-matching
      int isMatched = 0;
      int sel_mu_isMatched(0), pi_isMatched(0);
      int sel_mu_genIdx(-1), pi_genIdx(-1);
      int genMuonMother_genPdgId(-1), genPionMother_genPdgId(-1);
      float mupi_mass_reldiff(99.), lxy_reldiff(99.);
      float gen_lxy(-1.);

      // for MC only
      if(isMC_ == true){
        // pdgId of the gen particle to which the final-state particles are matched
        int sel_mu_genPdgId = lep_ptr->userInt("mcMatch");
        int pi_genPdgId     = pi_ptr->userInt("mcMatch");

       // index of the gen particle to which the final-state particles are matched
        sel_mu_genIdx   = lep_ptr->userInt("mcMatchIndex"); 
        pi_genIdx       = pi_ptr->userInt("mcMatchIndex"); 

        float mupi_mass_reco = fitter.fitted_candidate().mass(); // taking the fitted mass 
        float mupi_mass_gen = 99.;

        float hnl_vx_gen(99.), hnl_vy_gen(99.), mu_vx_gen(99.), mu_vy_gen(99.);

        //std::cout << "gen mu idx " << sel_mu_genIdx << " gen pi idx " << pi_genIdx << std::endl;

        if(sel_mu_genIdx != -1){
          // getting the associated gen particles
          edm::Ptr<reco::GenParticle> genMuon_ptr(genParticles, sel_mu_genIdx);

          // index of the associated mother particle
          int genMuonMother_genIdx = -1;
          if(genMuon_ptr->numberOfMothers()>0) genMuonMother_genIdx = genMuon_ptr->motherRef(0).key();

          // getting the mother particle
          edm::Ptr<reco::GenParticle> genMuonMother_ptr(genParticles, genMuonMother_genIdx);

          // fetching mass
          mupi_mass_gen = genMuonMother_ptr->mass();

          // pdgId of the mother particle
          genMuonMother_genPdgId = genMuonMother_ptr->pdgId();

          // getting vertices
          mu_vx_gen = genMuon_ptr->vx();
          mu_vy_gen = genMuon_ptr->vy();
          hnl_vx_gen = genMuonMother_ptr->vx();
          hnl_vy_gen = genMuonMother_ptr->vy();

          // matching of the displaced lepton
          if(fabs(sel_mu_genPdgId) == 13 && fabs(genMuonMother_genPdgId) == 9900015){
            sel_mu_isMatched = 1;
          }
        }

        if(pi_genIdx != -1){
          // getting the associated gen particles
          edm::Ptr<reco::GenParticle> genPion_ptr(genParticles, pi_genIdx);

          // index of the associated mother particle
          int genPionMother_genIdx = -1;
          if(genPion_ptr->numberOfMothers()>0) genPionMother_genIdx = genPion_ptr->motherRef(0).key();

          // getting the mother particles
          edm::Ptr<reco::GenParticle> genPionMother_ptr(genParticles, genPionMother_genIdx);

          // pdgId of the mother particle
          genPionMother_genPdgId = genPionMother_ptr->pdgId();

          // matching of the displaced pion
          if(fabs(pi_genPdgId) == 211 && fabs(genPionMother_genPdgId) == 9900015){
            pi_isMatched = 1;
          }
        }

        // computing displacement at gen level
        float gen_hnl_lxy = sqrt(pow(hnl_vx_gen - mu_vx_gen, 2) + pow(hnl_vy_gen - mu_vy_gen, 2));

        // computing relative difference between gen and reco quantities
        //mupi_mass_reldiff = fabs(mupi_mass_reco - mupi_mass_gen) / mupi_mass_gen; // using fitted tracks
        mupi_mass_reldiff = fabs(reco_hnl_mass - mupi_mass_gen) / mupi_mass_gen;
        lxy_reldiff = fabs(lxy.value() - gen_hnl_lxy) / gen_hnl_lxy;

        //std::cout << "mu matched " << sel_mu_isMatched << " pi matched " << pi_isMatched << " mass reldiff " << mupi_mass_reldiff << " lxy reldiff " << lxy_reldiff << std::endl;

        // matching of the full mulpi candidate
        //if(sel_mu_isMatched==1 && pi_isMatched==1 && mupi_mass_reldiff<0.3 && lxy_reldiff<0.5){
        if(sel_mu_isMatched==1 && pi_isMatched==1 && mupi_mass_reldiff<0.15 && lxy_reldiff<0.5){
        //if(sel_mu_isMatched==1 && pi_isMatched==1 && lxy_reldiff<0.5){
          isMatched = 1;
          //if(gen_mu_pt>1 && fabs(gen_mu_eta)<2.5 && gen_pi_pt>1 && fabs(gen_pi_eta)<2.5){
          //  std::cout << "fill h num " << std::endl;
          //  std::cout << "found a hnl pt " << gen_hnl_pt << " eta " << gen_hnl_eta << " mass " << gen_hnl_mass << " reco fitted mass " << mupi_mass_reco << " reco mass " << reco_hnl_mass << " vx " << hnl_cand.vx() << " vz " << hnl_cand.vz() << std::endl;
          //  std::cout << "found a muon pt " << lep_ptr->pt() << " eta " << lep_ptr->eta() << " isDSA " << lep_ptr->isDSAMuon() << std::endl;
          //  std::cout << "found a pion pt " << pi_ptr->pt() << " eta " << pi_ptr->eta() << std::endl;
          //}
        }
      }

      hnl_cand.addUserInt("isMatched", isMatched);
      hnl_cand.addUserInt("sel_mu_isMatched", sel_mu_isMatched);
      hnl_cand.addUserInt("pi_isMatched", pi_isMatched);
      hnl_cand.addUserInt("matching_sel_mu_genIdx", sel_mu_genIdx);
      hnl_cand.addUserInt("matching_pi_genIdx", pi_genIdx);
      hnl_cand.addUserInt("matching_sel_mu_motherPdgId", genMuonMother_genPdgId);
      hnl_cand.addUserInt("matching_pi_motherPdgId", genPionMother_genPdgId);
      hnl_cand.addUserFloat("mupi_mass_reco_gen_reldiff", mupi_mass_reldiff);
      hnl_cand.addUserFloat("lxy_reco_gen_reldiff", lxy_reldiff);
      hnl_cand.addUserFloat("gen_lxy", gen_lxy);

      ret_val->push_back(hnl_cand);

    } // for(size_t sel_mu_idx = 0; sel_mu_idx < sel_muons->size(); ++sel_mu_idx)

  } // for(size_t pi_idx = 0; pi_idx < kaons->size(); ++pi_idx)

  evt.put(std::move(ret_val));
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HNLToMuPiBuilder);
