// Merges the PF and LowPT collections, sets the isPF and isLowPt 
// UserInt's accordingly

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "ConversionInfo.h"

#include "PhysicsTools/BParkingNano/interface/helper.h"
#include <limits>
#include <algorithm>
#include "PhysicsTools/BParkingNano/interface/KinVtxFitter.h"
#include "PhysicsTools/BParkingNano/interface/ETHMuon.h"

class PFElectronMerger : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
  typedef std::vector<pat::ETHMuon> ETHMuonCollection;


public:
  bool debug=false; 

  explicit PFElectronMerger(const edm::ParameterSet &cfg):
    triggerMuons_{consumes<ETHMuonCollection>( cfg.getParameter<edm::InputTag>("trgMuon") )},
    pf_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("pfSrc") )},
    pf_mvaId_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("pfmvaId") )},
    vertexSrc_{ consumes<reco::VertexCollection> ( cfg.getParameter<edm::InputTag>("vertexCollection") )},
    conversions_{ consumes<edm::View<reco::Conversion> > ( cfg.getParameter<edm::InputTag>("conversions") )},
    beamSpot_{ consumes<reco::BeamSpot> ( cfg.getParameter<edm::InputTag>("beamSpot") )},
    drTrg_cleaning_{cfg.getParameter<double>("drForCleaning_wrtTrgMuon")},
    dzTrg_cleaning_{cfg.getParameter<double>("dzForCleaning_wrtTrgMuon")},
    dr_cleaning_{cfg.getParameter<double>("drForCleaning")},
    dz_cleaning_{cfg.getParameter<double>("dzForCleaning")},
    flagAndclean_{cfg.getParameter<bool>("flagAndclean")},
    pf_ptMin_{cfg.getParameter<double>("pf_ptMin")},
    ptMin_{cfg.getParameter<double>("ptMin")},
    etaMax_{cfg.getParameter<double>("etaMax")},
    bdtMin_{cfg.getParameter<double>("bdtMin")},
    use_gsf_mode_for_p4_{cfg.getParameter<bool>("useGsfModeForP4")},
    use_regression_for_p4_{cfg.getParameter<bool>("useRegressionModeForP4")},
    sortOutputCollections_{cfg.getParameter<bool>("sortOutputCollections")},
    saveLowPtE_{cfg.getParameter<bool>("saveLowPtE")},
    addUserVarsExtra_{cfg.getParameter<bool>("addUserVarsExtra")}
    {
       produces<pat::ElectronCollection>("SelectedElectrons");
       produces<TransientTrackCollection>("SelectedTransientElectrons");  
    }

  ~PFElectronMerger() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const edm::EDGetTokenT<ETHMuonCollection> triggerMuons_;
  const edm::EDGetTokenT<pat::ElectronCollection> pf_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> pf_mvaId_src_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  const edm::EDGetTokenT<edm::View<reco::Conversion> > conversions_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  const double drTrg_cleaning_;
  const double dzTrg_cleaning_;
  const double dr_cleaning_;
  const double dz_cleaning_;
  const bool flagAndclean_;
  const double pf_ptMin_;
  const double ptMin_; //pt min cut
  const double etaMax_; //eta max cut
  const double bdtMin_; //bdt min cut
  const bool use_gsf_mode_for_p4_;
  const bool use_regression_for_p4_;
  const bool sortOutputCollections_;
  const bool saveLowPtE_;
  const bool addUserVarsExtra_;

};

void PFElectronMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const & iSetup) const {

  //input
  edm::Handle<ETHMuonCollection> trgMuon;
  evt.getByToken(triggerMuons_, trgMuon);
  edm::Handle<pat::ElectronCollection> pf;
  evt.getByToken(pf_src_, pf);
  edm::Handle<edm::ValueMap<float> > pfmvaId;  
  evt.getByToken(pf_mvaId_src_, pfmvaId);
  // 
  edm::ESHandle<TransientTrackBuilder> theB ;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  //
  edm::Handle<reco::VertexCollection> vertexHandle;
  evt.getByToken(vertexSrc_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();
  //
  edm::Handle<edm::View<reco::Conversion> > conversions;
  evt.getByToken(conversions_, conversions);
  edm::Handle<reco::BeamSpot> beamSpot;
  evt.getByToken(beamSpot_, beamSpot);

  // output
  std::unique_ptr<pat::ElectronCollection>  ele_out      (new pat::ElectronCollection );
  std::unique_ptr<TransientTrackCollection> trans_ele_out(new TransientTrackCollection);
  std::vector<std::pair<float, float>> pfEtaPhi;
  std::vector<float> pfVz;
  
  // -> changing order of loops ert Arabella's fix this without need for more vectors  
  size_t ipfele=-1;
  for(auto ele : *pf) {
   ipfele++;

   if (debug) std::cout << "ElectronMerger, Event " << (evt.id()).event() 
			<< " => PF: ele.superCluster()->rawEnergy() = " << ele.superCluster()->rawEnergy()
			<< ", ele.correctedEcalEnergy() = " << ele.correctedEcalEnergy()
			<< ", ele gsf track chi2 = " << ele.gsfTrack()->normalizedChi2()
			<< ", ele.p = " << ele.p() << std::endl;

   //cuts
   if (ele.pt()<ptMin_ || ele.pt() < pf_ptMin_) continue;
   if (fabs(ele.eta())>etaMax_) continue;
   // apply conversion veto unless we want conversions
   if (!ele.passConversionVeto()) continue;

   // Fix the mass to the proper one
   reco::Candidate::PolarLorentzVector p4( 
					  ele.pt(),
					  ele.eta(),
					  ele.phi(),
					  ELECTRON_MASS
					   );
   ele.setP4(p4);     

   // skip electrons inside tag's jet or from different PV
   bool skipEle=true;
   float dzTrg = 0.0;
   for(const auto & trg : *trgMuon) {
   if(trg.userInt("isTriggeringBPark") != 1) continue;	
     if(reco::deltaR(ele, trg) < drTrg_cleaning_ && drTrg_cleaning_ > 0)
        continue;
     if(fabs(ele.vz() - trg.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ > 0)
        continue;
     skipEle=false;
     dzTrg = ele.vz() - trg.vz();
     break; // one trg muon to pass is enough :)
   }
   // we skip evts without trg muon
   if (skipEle) continue;

   // for PF e we set BDT outputs to much higher number than the max
   edm::Ref<pat::ElectronCollection> ref(pf,ipfele);
   float pf_mva_id = float((*pfmvaId)[ref]);
   ele.addUserInt("isPF", 1);
   ele.addUserInt("isLowPt", 0);
   ele.addUserFloat("ptBiased", 20.);
   ele.addUserFloat("unBiased", 20.);
   ele.addUserFloat("mvaId", 20.);
   ele.addUserFloat("pfmvaId", pf_mva_id);
   ele.addUserFloat("chargeMode", ele.charge());
   ele.addUserInt("isPFoverlap", 0);
   ele.addUserFloat("dzTrg", dzTrg);
   ele.addUserFloat("ip3d", fabs(ele.dB(pat::Electron::PV3D)));
   ele.addUserFloat("sip3d", fabs(ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D)));
   ele.addUserFloat("dxy", fabs(ele.dB(pat::Electron::PV2D)));
   ele.addUserFloat("dz", fabs(ele.dB(pat::Electron::PVDZ)));

   // Attempt to match electrons to conversions in "gsfTracksOpenConversions" collection (NO MATCHES EXPECTED)
   ConversionInfo info;
   ConversionInfo::match(beamSpot,conversions,ele,info);
   info.addUserVars(ele);
   if ( addUserVarsExtra_ ) { info.addUserVarsExtra(ele); }

   pfEtaPhi.push_back(std::pair<float, float>(ele.eta(), ele.phi()));
   pfVz.push_back(ele.vz());
   ele_out       -> emplace_back(ele);
  }

  unsigned int pfSelectedSize = pfEtaPhi.size();

  if(sortOutputCollections_){

    //sorting increases sligtly the time but improves the code efficiency in the Bcandidate builder
    //easier identification of leading and subleading with smarter loop
    std::sort( ele_out->begin(), ele_out->end(), [] (pat::Electron e1, pat::Electron e2) -> bool {return e1.pt() > e2.pt();}
             );
  }

  // build transient track collection
  for(auto &ele : *ele_out){
    float regErrorRatio = std::abs(ele.corrections().combinedP4Error/ele.p()/ele.gsfTrack()->qoverpModeError()*ele.gsfTrack()->qoverpMode());
    const reco::TransientTrack eleTT = use_regression_for_p4_ ?
      (*theB).buildfromReg(ele.gsfTrack(), math::XYZVector(ele.corrections().combinedP4), regErrorRatio) : (*theB).buildfromGSF( ele.gsfTrack() );
    trans_ele_out -> emplace_back(eleTT);

  } 
  //adding label to be consistent with the muon and track naming
  evt.put(std::move(ele_out),      "SelectedElectrons");
  evt.put(std::move(trans_ele_out),"SelectedTransientElectrons");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFElectronMerger);
