// Code to apply energy regression

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"

#include "PhysicsTools/BParkingNano/interface/helper.h"


class PFElectronRegresser : public edm::global::EDProducer<> {

public:
  bool debug=false; 

  explicit PFElectronRegresser(const edm::ParameterSet &cfg):
    pf_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("pfSrc") )}
    {


      // PF regression                                                                                                           
      if( cfg.existsAs<edm::ParameterSet>("gsfRegressionConfig") ) {
	const edm::ParameterSet& iconf = cfg.getParameterSet("gsfRegressionConfig");
	const std::string& mname = iconf.getParameter<std::string>("modifierName");
      ModifyObjectValueBase* plugin =
        ModifyObjectValueFactory::get()->create(mname,iconf);
      regressionGsf_.reset(plugin);
      edm::ConsumesCollector sumes = consumesCollector();
      regressionGsf_->setConsumes(sumes);
      } else {
	regressionGsf_.reset(nullptr);
      }

      produces<pat::ElectronCollection>("regressedElectrons");
    }

  ~PFElectronRegresser() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const edm::EDGetTokenT<pat::ElectronCollection> pf_src_;

  // regression stuff                                                                                                                                         
  std::unique_ptr<ModifyObjectValueBase> regressionGsf_; // Gsf                                                                                               
};

void PFElectronRegresser::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const & iSetup) const {

  //input
  edm::Handle<pat::ElectronCollection> pf;
  evt.getByToken(pf_src_, pf);

  // regression stuff
  regressionGsf_->setEvent(evt);
  regressionGsf_->setEventContent(iSetup);

  // output
  std::unique_ptr<pat::ElectronCollection>  ele_out_pf      (new pat::ElectronCollection );

  // PF regression
  size_t ipfele=-1;
  for(auto ele : *pf) {
   ipfele++;
   if(debug) {

     std::cout << "PFElectronRegresser, Event " << (evt.id()).event() 

	       << " => Pre regression, PF: ele.superCluster()->rawEnergy() = " << ele.superCluster()->rawEnergy()
	       << ", ele.correctedEcalEnergy() = " << ele.correctedEcalEnergy()
	       << ", ele gsf track chi2 = " << ele.core()->gsfTrack()->normalizedChi2()
	       << ", ele.p = " << ele.p() << std::endl;
   }

   regressionGsf_->modifyObject(ele);

   if(debug) { 

     std::cout << "PFElectronRegresser, Event " << (evt.id()).event() 

	       << " => Post regression, PF: ele.superCluster()->rawEnergy() = " << ele.superCluster()->rawEnergy()
	       << ", ele.correctedEcalEnergy() = " << ele.correctedEcalEnergy()
	       << ", ele gsf track chi2 = " << ele.core()->gsfTrack()->normalizedChi2()
	       << ", ele.p = " << ele.p() << std::endl;
   }
     
   ele_out_pf -> emplace_back(ele);
  }

  // put collections in the event

  evt.put(std::move(ele_out_pf),  "regressedElectrons");
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFElectronRegresser);
