#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "../interface/helper.h"
#include <limits>
#include <algorithm>
#include "../interface/KinVtxFitter.h"
#include "../interface/ETHMuon.h"

//template<typename Lepton>
class TagAndProbeJPsiToMuMuBuilder : public edm::global::EDProducer<> {

public:
  typedef pat::ETHMuon Lepton;
  typedef std::vector<Lepton> LeptonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit TagAndProbeJPsiToMuMuBuilder(const edm::ParameterSet &cfg):
    l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
    l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src_{consumes<LeptonCollection>( cfg.getParameter<edm::InputTag>("src") )},
    isMC_{cfg.getParameter<bool>("isMC")},
    ttracks_src_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracksSrc") )},
    genParticles_ {consumes<reco::GenParticleCollection>( cfg.getParameter<edm::InputTag>("genParticles") )}, 
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot"               ) )} { 
       produces<pat::CompositeCandidateCollection>();
    }

  ~TagAndProbeJPsiToMuMuBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<Lepton> l1_selection_; // cut on leading lepton
  const StringCutObjectSelector<Lepton> l2_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<LeptonCollection> src_;
  const bool isMC_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

//template<typename Lepton>
void TagAndProbeJPsiToMuMuBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<LeptonCollection> leptons;
  evt.getByToken(src_, leptons);
  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_src_, ttracks);

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genParticles_, genParticles);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  
  // loop on the tag muon
  for(size_t l1_idx = 0; l1_idx < leptons->size(); ++l1_idx) {
    edm::Ptr<Lepton> l1_ptr(leptons, l1_idx);
    if(l1_ptr->isDSAMuon()==1) continue;
    // the next two cuts are removed from the builder and applied in the ntupliser
    ////if(l1_ptr->userInt("HLT_Mu8_v1")!=1) continue;
    ////if(l1_ptr->userInt("HLT_Mu8_v1_prescale")<1) continue;
    if(l1_ptr->userInt("isTriggering")!=1) continue;
    if(l1_ptr->softId()!=1) continue;
    if(!l1_selection_(*l1_ptr)) continue; 
    
    // loop on the probe muon
    for(size_t l2_idx = l1_idx + 1; l2_idx < leptons->size(); ++l2_idx) {
      edm::Ptr<Lepton> l2_ptr(leptons, l2_idx);
      if(l2_ptr->isDSAMuon()==1) continue;
      if(l2_ptr->charge() == l1_ptr->charge()) continue;
      // request at least one BParking trigger to be fired (cut moved to the ntupliser)
      ////if(l2_ptr->userInt("HLT_Mu7_IP4")!=1 && l2_ptr->userInt("HLT_Mu8_IP3")!=1 && l2_ptr->userInt("HLT_Mu8_IP5")!=1 && l2_ptr->userInt("HLT_Mu8_IP6")!=1 && l2_ptr->userInt("HLT_Mu8p5_IP3p5")!=1 && l2_ptr->userInt("HLT_Mu9_IP4")!=1 && l2_ptr->userInt("HLT_Mu9_IP5")!=1 && l2_ptr->userInt("HLT_Mu9_IP6")!=1 && l2_ptr->userInt("HLT_Mu10p5_IP3p5")!=1 && l2_ptr->userInt("HLT_Mu12_IP6")!=1) continue;
      ////if(l2_ptr->userInt("isTriggering")!=1) continue;
      if(!l2_selection_(*l2_ptr)) continue;

      pat::CompositeCandidate lepton_pair;
      lepton_pair.setP4(l1_ptr->polarP4() + l2_ptr->polarP4());
      lepton_pair.setCharge(l1_ptr->charge() + l2_ptr->charge());
      lepton_pair.addUserFloat("lep_deltaR", reco::deltaR(*l1_ptr, *l2_ptr));
      lepton_pair.addUserFloat("lep_vzdiff", l1_ptr->vz() - l2_ptr->vz()); 
      
      // Put the lepton passing the corresponding selection
      lepton_pair.addUserInt("l1_idx", l1_idx );
      lepton_pair.addUserInt("l2_idx", l2_idx );
      lepton_pair.addUserCand("l1", l1_ptr );
      lepton_pair.addUserCand("l2", l2_ptr );

      //if(isMC_){
      //  lepton_pair.addUserInt("l1_mcMatch", l1_ptr->userInt("mcMatch"));
      //  lepton_pair.addUserInt("l2_mcMatch", l2_ptr->userInt("mcMatch"));
      //  lepton_pair.addUserInt("l1_mcMatchIndex", l1_ptr->userInt("mcMatchIndex"));
      //  lepton_pair.addUserInt("l2_mcMatchIndex", l2_ptr->userInt("mcMatchIndex"));
      //}

      if( !pre_vtx_selection_(lepton_pair) ) continue; // before making the SV, cut on the info we have

      KinVtxFitter fitter(
        {ttracks->at(l1_idx), ttracks->at(l2_idx)},
        {l1_ptr->mass(), l2_ptr->mass()},
        {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
        );

      if(!fitter.success()) continue;

      lepton_pair.addUserFloat("sv_chi2", fitter.chi2());
      lepton_pair.addUserFloat("sv_ndof", fitter.dof());
      lepton_pair.addUserFloat("sv_prob", fitter.prob());
      lepton_pair.addUserFloat("fitted_mass", fitter.success() ? fitter.fitted_candidate().mass() : -1);
      lepton_pair.addUserFloat("fitted_massErr", fitter.success() ? sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)) : -1);
      lepton_pair.addUserFloat("charge", lepton_pair.charge());
      lepton_pair.addUserFloat("fitted_lep1_pt", fitter.daughter_p4(0).pt()); 
      lepton_pair.addUserFloat("fitted_lep1_eta", fitter.daughter_p4(0).eta());
      lepton_pair.addUserFloat("fitted_lep1_phi", fitter.daughter_p4(0).phi());
      lepton_pair.addUserFloat("fitted_lep2_pt", fitter.daughter_p4(1).pt()); 
      lepton_pair.addUserFloat("fitted_lep2_eta", fitter.daughter_p4(1).eta());
      lepton_pair.addUserFloat("fitted_lep2_phi", fitter.daughter_p4(1).phi());
        
      auto lxy    = l_xy(fitter, *beamspot);
      lepton_pair.addUserFloat("l_xy" , lxy.value());
      lepton_pair.addUserFloat("l_xy_sig", lxy.value()/lxy.error());

      auto fit_p4 = fitter.fitted_p4();
      lepton_pair.addUserFloat("fitted_pt", fit_p4.pt()); 
      lepton_pair.addUserFloat("fitted_eta", fit_p4.eta());
      lepton_pair.addUserFloat("fitted_phi", fit_p4.phi());
      lepton_pair.addUserFloat("cos_theta_2D", cos_theta_2D(fitter, *beamspot, fit_p4));

      // cut on the SV info
      if( !post_vtx_selection_(lepton_pair) ) continue;
      
      // gen-matching
      int isMatched = 0;
      int l1_genIdx(-1), l2_genIdx(-1);
      int genMuon1Mother_genPdgId(-1), genMuon2Mother_genPdgId(-1);

      // for MC only
      if(isMC_ == true){

        // pdgId of the gen particle to which the final-state particles are matched
        int l1_genPdgId = l1_ptr->userInt("mcMatch");
        int l2_genPdgId = l2_ptr->userInt("mcMatch");
        
        // index of the gen particle to which the final-state particles are matched
        l1_genIdx = l1_ptr->userInt("mcMatchIndex"); 
        l2_genIdx = l2_ptr->userInt("mcMatchIndex"); 

        if(l1_genIdx != -1 && l2_genIdx != -1){

          // getting the associated gen particles
          edm::Ptr<reco::GenParticle> genMuon1_ptr(genParticles, l1_genIdx);
          edm::Ptr<reco::GenParticle> genMuon2_ptr(genParticles, l2_genIdx);

          // index of the associated mother particle
          int genMuon1Mother_genIdx = -1;
          int genMuon2Mother_genIdx = -1;
          if(genMuon1_ptr->numberOfMothers()>0) genMuon1Mother_genIdx = genMuon1_ptr->motherRef(0).key();
          if(genMuon2_ptr->numberOfMothers()>0) genMuon2Mother_genIdx = genMuon2_ptr->motherRef(0).key();

          // getting the mother particles
          edm::Ptr<reco::GenParticle> genMuon1Mother_ptr(genParticles, genMuon1Mother_genIdx);
          edm::Ptr<reco::GenParticle> genMuon2Mother_ptr(genParticles, genMuon2Mother_genIdx);

          // pdgId of the mother particles
          genMuon1Mother_genPdgId = genMuon1Mother_ptr->pdgId();
          genMuon2Mother_genPdgId = genMuon2Mother_ptr->pdgId();

          if(
             fabs(l1_genPdgId) == 13 && fabs(genMuon1Mother_genPdgId) == 443 && 
             fabs(l2_genPdgId) == 13 && fabs(genMuon2Mother_genPdgId) == 443 &&
             genMuon1Mother_genIdx == genMuon2Mother_genIdx
            ){
              isMatched = 1;
          }
        }
      }

      lepton_pair.addUserInt("isMatched", isMatched);

      ret_value->push_back(lepton_pair);
    }
  }
  
  evt.put(std::move(ret_value));
}

//#include "../interface/ETHMuon.h"
//typedef TagAndProbeJPsiToLLBuilder<pat::ETHMuon> TagAndProbeJPsiToMuMuBuilder;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TagAndProbeJPsiToMuMuBuilder);
