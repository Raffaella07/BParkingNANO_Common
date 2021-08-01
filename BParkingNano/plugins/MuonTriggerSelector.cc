// class to produce 2 pat::MuonCollections
// one matched to the Park triggers
// another fitered wrt the Park triggers


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "ETHMuon.h"

#include <TLorentzVector.h>
#include "helper.h"

using namespace std;

constexpr bool debug = false; //false;

class MuonTriggerSelector : public edm::EDProducer {

  public:

    explicit MuonTriggerSelector(const edm::ParameterSet &iConfig);

    ~MuonTriggerSelector() override {};


  private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<std::vector<reco::Track>> displacedStandaloneMuonSrc_;
    const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

    // added
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

    // trigger muon matching
    const double max_deltaR_;
    const double max_deltaPtRel_;

    // for the sel muon
    const double selmu_ptMin_;          // min pT in all muons for B candidates
    const double selmu_absEtaMax_;      // max eta ""
    std::vector<std::string> HLTPaths_;
    //std::vector<std::string> L1Seeds_;
};


MuonTriggerSelector::MuonTriggerSelector(const edm::ParameterSet &iConfig):
  muonSrc_(consumes<std::vector<pat::Muon>> (iConfig.getParameter<edm::InputTag>("muonCollection"))),
  displacedStandaloneMuonSrc_(consumes<std::vector<reco::Track>> (iConfig.getParameter<edm::InputTag>("displacedStandaloneMuonCollection"))),
  vertexSrc_( consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>( "vertexCollection"))),
  // added
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  //
  max_deltaR_(iConfig.getParameter<double>("max_deltaR")),
  max_deltaPtRel_(iConfig.getParameter<double>("max_deltaPtRel")),
  selmu_ptMin_(iConfig.getParameter<double>("selmu_ptMin")),
  selmu_absEtaMax_(iConfig.getParameter<double>("selmu_absEtaMax")),
  HLTPaths_(iConfig.getParameter<std::vector<std::string>>("HLTPaths"))
  //L1Seeds_(iConfig.getParameter<std::vector<std::string>>("L1seeds"))
{
  // produce 2 collections: trgMuons (tags) and SelectedMuons (probes & tags if survive preselection cuts)
  // trigger muons are slimmedMuons only
  // selected muons are slimmedMuons and displacedStandaloneMuons
  // make sure that the indices between the selected muons and transient tracks are consistent
  produces<pat::MuonCollection>("trgMuons"); 
  produces<std::vector<pat::ETHMuon>>("SelectedMuons");
  produces<TransientTrackCollection>("SelectedTransientMuons");  
}


void MuonTriggerSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::ESHandle<MagneticField> bFieldHandle;
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

    std::unique_ptr<pat::MuonCollection>      trgmuons_out   ( new pat::MuonCollection );
    std::unique_ptr<std::vector<pat::ETHMuon>> ETHmuons_out  ( new std::vector<pat::ETHMuon> );
    std::unique_ptr<TransientTrackCollection> trans_muons_out( new TransientTrackCollection );
    
    edm::Handle<std::vector<pat::Muon>> slimmed_muons;
    iEvent.getByToken(muonSrc_, slimmed_muons);

    edm::Handle<std::vector<reco::Track>> displaced_standalone_muons;
    iEvent.getByToken(displacedStandaloneMuonSrc_, displaced_standalone_muons);

    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByToken(vertexSrc_, vertexHandle);

    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(triggerPrescales_, triggerPrescales);

    std::vector<int> muonIsTrigger(slimmed_muons->size(), 0);
    std::vector<float> muonDR(slimmed_muons->size(),10000.);
    std::vector<float> muonDPT(slimmed_muons->size(),10000.);
    std::vector<int> loose_id(slimmed_muons->size(),0);

    std::vector<int> matched_reco_flag(slimmed_muons->size(),-1);
    std::vector<int> matched_trg_index(slimmed_muons->size(),-1);
    std::vector<float> matched_dr(slimmed_muons->size(),10000.);
    std::vector<float> matched_dpt(slimmed_muons->size(),-10000.);
    std::vector<std::vector<int>> fires;
    std::vector<std::vector<int>> prescales;
    std::vector<std::vector<float>> matcher; 
    std::vector<std::vector<float>> DR;
    std::vector<std::vector<float>> DPT;    

    // fetching the primary vertex
    const reco::Vertex & PV = vertexHandle->front();

    //std::cout << std::endl;

    //////////////////////////////
    /* method 1 */
    // proceeding to the trigger muon matching
    // for that, only use slimmedMuons
    for(const pat::Muon &muon : *slimmed_muons){
        if(debug) std::cout <<"Muon Pt="<< muon.pt() << " Eta=" << muon.eta() << " Phi=" << muon.phi() << endl;
        //std::cout <<"Muon Pt="<< muon.pt() << " Eta=" << muon.eta() << " Phi=" << muon.phi() << " nTriggerObjectMatch " << muon.triggerObjectMatches().size() << endl;

        std::vector<int> frs(HLTPaths_.size(),0); //path fires for each reco muon
        std::vector<int> prescale(HLTPaths_.size(),-1);
        //std::vector<int> sds(L1Seeds_.size(),0);// L1 Seeds for each L1 muon
        std::vector<float> temp_matched_to(HLTPaths_.size(),1000.);
        std::vector<float> temp_DR(HLTPaths_.size(),1000.);
        std::vector<float> temp_DPT(HLTPaths_.size(),1000.);
        int ipath=-1;

/*        int iseed=-1;
        for (const std::string seed: L1Seeds_){
            iseed++;
            char cstr[(seed+"*").size()+1];
            strcpy( cstr,(seed+"*").c_str());
            if(muon.triggerObjectMatches().size()!=0){
                for(size_t i=0;i<muon.triggerObjectMatches().size(); i++){
                    if(muon.triggerObjectMatch(i)!=0 && muon.triggerObjectMatch(i)->hasAlgorithmName(cstr,true)){
                        sds[iseed]=1;
                        std::cout << "L1 Seed="<< cstr <<" fired="<< sds[iseed] << endl;
                    }
                }
            }
        } */

        // for each muon, we study whether it fires a HLT line or not
        for (const std::string path: HLTPaths_){
            if(debug) std::cout << path << " " << muon.triggerObjectMatches().size() << std::endl;
            ipath++;

            // the following vectors are used in order to find the minimum DR between a reco muon and all the HLT objects that is matched to it
            // as a reco muon will be matched with only one HLT object every time, so there is a one-to-one correspondence between the two collections
            // DPt_rel is not used to create this one-to-one correspondence but only to create a few plots, debugging and be sure that everything is working fine. 
            std::vector<float> temp_dr(muon.triggerObjectMatches().size(),1000.);
            std::vector<float> temp_dpt(muon.triggerObjectMatches().size(),1000.);
            std::vector<float> temp_pt(muon.triggerObjectMatches().size(),1000.);

            // put the HLT path name into a wildcard
            char cstr[ (path+"*").size() + 1];
            strcpy( cstr, (path+"*").c_str() );       

            // Here we find all the HLT objects from each HLT path each time that are matched with the reco muon.
            if(muon.triggerObjectMatches().size()!=0){
                for(size_t i=0; i<muon.triggerObjectMatches().size();i++){
                  if(debug) std::cout << "triggerObjMatch " << i << " " << cstr << " " <<  muon.triggerObjectMatch(i)->hasPathName(cstr,true,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,false) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,true,false) << std::endl;
                    //std::cout << "trgObjMatch " << i << " " << cstr << " " <<  muon.triggerObjectMatch(i)->hasPathName(cstr,true,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,false) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,true,false) << std::endl;
                    // first bool is pathLastFilterAccepted, second is pathL3FilterAccepted
                    if(muon.triggerObjectMatch(i)!=0 && muon.triggerObjectMatch(i)->hasPathName(cstr,false,false)){
                    //if(muon.triggerObjectMatch(i)!=0 && muon.triggerObjectMatch(i)->hasPathName(cstr,false,false) && DR[iMuo][path]<max_deltaR_ && fabs(DPT[iMuo][path])<max_deltaPtRel_ && DR[iMuo][path]!=10000){
                        //if(abs(muon.triggerObjectMatch(i)->eta())>1.5) std::cout << "eta=" <<muon.triggerObjectMatch(i)->eta();
                        frs[ipath]=1;
                        float dr = reco::deltaR(muon.triggerObjectMatch(i)->p4(), muon.p4()); 
                        float dpt=(muon.triggerObjectMatch(i)->pt()-muon.pt())/muon.triggerObjectMatch(i)->pt();
                        temp_dr[i]=dr;
                        temp_dpt[i]=dpt;
                        temp_pt[i]=muon.triggerObjectMatch(i)->pt();                   
                        if(debug)std::cout <<"Path=" <<cstr << endl;
                        if(debug)std::cout <<"HLT  Pt="<<muon.triggerObjectMatch(i)->pt() <<" Eta="<<muon.triggerObjectMatch(i)->eta() <<" Phi="<<muon.triggerObjectMatch(i)->phi() << endl;
                        //std::cout <<"HLT " << cstr << " Pt="<<muon.triggerObjectMatch(i)->pt() <<" Eta="<<muon.triggerObjectMatch(i)->eta() <<" Phi="<<muon.triggerObjectMatch(i)->phi() << endl;
                        if(debug)std::cout <<"Muon Pt="<< muon.pt() << " Eta=" << muon.eta() << " Phi=" << muon.phi()  <<endl;
                        if(debug)std::cout <<"DR = " << temp_dr[i] <<endl;
                        //std::cout <<"DR = " << temp_dr[i] << " DPT = " << temp_dpt[i] <<endl;
                    }
                }
                // and now we find the real minimum between the reco muon and all its matched HLT objects. 
                temp_DR[ipath]=*min_element(temp_dr.begin(),temp_dr.end());
                int position=std::min_element(temp_dr.begin(),temp_dr.end()) - temp_dr.begin();
                temp_DPT[ipath]=temp_dpt[position];
                temp_matched_to[ipath]=temp_pt[position];
                }
            }
        //and now since we have found the minimum DR we save a few variables for plots       
        fires.push_back(frs); //This is used in order to see if a reco muon fired a Trigger (1) or not (0).
        //for(unsigned int i(0); i<frs.size(); ++i){
        //  std::cout << "fires " << frs[i] << " deltaR " << temp_DR[i] << " deltaPt " << temp_DPT[i] << std::endl;
        //}
        prescales.push_back(prescale); //This is used for the initialisation of the prescales.
        matcher.push_back(temp_matched_to); //This is used in order to see if a reco muon is matched with a HLT object. PT of the reco muon is saved in this vector. 
        DR.push_back(temp_DR);
        DPT.push_back(temp_DPT);

    }

    //now, check for different reco muons that are matched to the same HLTObject.
    for(unsigned int path=0; path<HLTPaths_.size(); path++){
        for(unsigned int iMuo=0; iMuo<slimmed_muons->size(); iMuo++){
            for(unsigned int im=(iMuo+1); im<slimmed_muons->size(); im++){
                if(matcher[iMuo][path]!=1000. && matcher[iMuo][path]==matcher[im][path]){
                //if(matcher[iMuo][path]!=1000. && matcher[iMuo][path]==matcher[im][path] && DR[iMuo][path]<max_deltaR_ && fabs(DPT[iMuo][path])<max_deltaPtRel_ && DR[im][path]<max_deltaR_ && fabs(DPT[im][path])<max_deltaPtRel_){
                    if(DR[iMuo][path]<DR[im][path]){ //Keep the one that has the minimum DR with the HLT object
                        fires[im][path]=0;
                        matcher[im][path]=1000.;
                        DR[im][path]=1000.;                       
                        DPT[im][path]=1000.;
                    }
                    else{
                        fires[iMuo][path]=0;
                        matcher[iMuo][path]=1000.;
                        DR[iMuo][path]=1000.;                       
                        DPT[iMuo][path]=1000.;
                    }
                }              
            }
            if(matcher[iMuo][path]!=1000. && DR[iMuo][path]<max_deltaR_ && fabs(DPT[iMuo][path])<max_deltaPtRel_ && DR[iMuo][path]!=10000){
              muonIsTrigger[iMuo]=1;
              muonDR[iMuo]=DR[iMuo][path];
              muonDPT[iMuo]=DPT[iMuo][path];                
              
              for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
                //if(names.triggerName(i).find(HLTPaths_[path]) != std::string::npos || names.triggerName(i)==HLTPaths_[path]){
                if(names.triggerName(i).find(HLTPaths_[path]) != std::string::npos){
                  prescales[iMuo][path] = triggerPrescales->getPrescaleForIndex(i);
                }
              }
            }
            else{
              fires[iMuo][path]=0;
              //prescales[iMuo][path] = -1;
            }
        }
    }

    //////////////////////////////
    /* method 2 */

    /*
    // we first get the trigger objects
    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);


    // container for the trigger objects, that will be later on matched to reco muons
    std::vector<pat::TriggerObjectStandAlone> triggeringMuons;


    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackFilterLabels(iEvent, *triggerBits);
      obj.unpackPathNames(names);

      bool isTriggerMuon = false;
      for (unsigned h = 0; h < obj.filterIds().size(); ++h){
        if(obj.filterIds()[h] == 83){ 
          isTriggerMuon = true; 
          break;
        } 
      }

      if(!isTriggerMuon) continue; 
      for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
        std::string filterName = obj.filterLabels()[h];
        if(filterName.find("hltL3") != std::string::npos  && filterName.find("Park") != std::string::npos){
          isTriggerMuon = true;
          if(debug) std::cout << "\t   Filters:   " << filterName; 
          break;
        }
        else{ isTriggerMuon = false; }
      }

      if(!isTriggerMuon) continue;

      std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      // Print trigger object collection and type
      std::cout << "\t   Collection: " << obj.collection() << std::endl;
      std::cout << "\t   Type IDs:   ";

      for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
      std::cout << std::endl;

      for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
      std::cout << std::endl;

      std::vector pathNamesAll = obj.pathNames(false);
      std::vector pathNamesLast = obj.pathNames(true);

      std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";

      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
            bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
            bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
            bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
            bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
            std::cout << "   " << pathNamesAll[h];
            if (isBoth) std::cout << "(L,3)";
            if (isL3 && !isBoth) std::cout << "(*,3)";
            if (isLF && !isBoth) std::cout << "(L,*)";
            if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
      }
      std::cout << std::endl;

      triggeringMuons.push_back(obj);
      if(debug){ std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
        // Print trigger object collection and type
        std::cout << "\t   Collection: " << obj.collection() << std::endl;
      }
    }//trigger objects
    */


    if(debug)std::cout << "number of Muons=" <<slimmed_muons->size() << endl;
    //And now create a collection with all trg muons
    for(const pat::Muon & muon : *slimmed_muons){
        unsigned int iMuo(&muon -&(slimmed_muons->at(0)));
        if(muonIsTrigger[iMuo]==1){
            pat::Muon recoTriggerMuonCand(muon);
            trgmuons_out->emplace_back(recoTriggerMuonCand);
        }
    }

    // add the displaced standalone muons to the collection
    for(const reco::Track & muon : *displaced_standalone_muons){
      if(muon.pt()<selmu_ptMin_) continue;
      if(fabs(muon.eta())>selmu_absEtaMax_) continue;

      // add the muon to the transient track collection
      // one has to make sure that the indices in the muon and transient track collections match
      const reco::TransientTrack muonTT((*(&muon)),&(*bFieldHandle)); 
      if(!muonTT.isValid()) continue; 

      pat::ETHMuon the_muon(muon);
      ETHmuons_out->emplace_back(the_muon);

      for(unsigned int i=0; i<HLTPaths_.size(); i++){
        ETHmuons_out->back().addUserInt(HLTPaths_[i], -1);
        ETHmuons_out->back().addUserInt(HLTPaths_[i] + "_prescale", -1);
      }
      ETHmuons_out->back().addUserInt("isTriggering", -1);
      ETHmuons_out->back().addUserFloat("DR", -1.);
      ETHmuons_out->back().addUserFloat("DPT", -1.);
      ETHmuons_out->back().addUserFloat("dz", muon.dz(PV.position()));
      ETHmuons_out->back().addUserFloat("dzS", muon.dz(PV.position())/muon.dzError());
      ETHmuons_out->back().addUserFloat("dxy", muon.dxy(PV.position()));
      ETHmuons_out->back().addUserFloat("dxyS", muon.dxy(PV.position())/muon.dxyError());

      trans_muons_out->emplace_back(muonTT);
    }

    // add the slimmed muons to the collection 
    for(const pat::Muon & muon : *slimmed_muons){
      unsigned int iMuo(&muon - &(slimmed_muons->at(0)) );
      if(muon.pt()<selmu_ptMin_) continue;
      if(fabs(muon.eta())>selmu_absEtaMax_) continue;

      const reco::TransientTrack muonTT((*(muon.bestTrack())),&(*bFieldHandle)); //sara:check,why not using inner track for muons? GM: What is this and why do we need this???
      if(!muonTT.isValid()) continue; // GM: and why do we skip this muon if muonTT is invalid? This seems to have no effect so I kept it.

      //std::cout << "muon " << iMuo << " pt " << muon.pt()  << " eta " << muon.eta() << std::endl;
      
      pat::ETHMuon the_muon(muon);
      ETHmuons_out->emplace_back(the_muon);

      std::cout << "is triggering " << muonIsTrigger[iMuo] << std::endl;
      for(unsigned int i=0; i<HLTPaths_.size(); i++){
        ETHmuons_out->back().addUserInt(HLTPaths_[i], fires[iMuo][i]);
        ETHmuons_out->back().addUserInt(HLTPaths_[i] + "_prescale", prescales[iMuo][i]);
      }
      ETHmuons_out->back().addUserInt("isTriggering", muonIsTrigger[iMuo]);
      ETHmuons_out->back().addUserFloat("DR", muonDR[iMuo]);
      ETHmuons_out->back().addUserFloat("DPT" ,muonDPT[iMuo]);
      ETHmuons_out->back().addUserFloat("dz", muon.dB(muon.PVDZ));
      ETHmuons_out->back().addUserFloat("dzS", muon.dB(muon.PVDZ)/muon.edB(muon.PVDZ));
      ETHmuons_out->back().addUserFloat("dxy", muon.dB(muon.PV2D));
      ETHmuons_out->back().addUserFloat("dxyS", muon.dB(muon.PV2D)/muon.edB(muon.PV2D));

      trans_muons_out->emplace_back(muonTT);
    }

    iEvent.put(std::move(trgmuons_out), "trgMuons"); 
    iEvent.put(std::move(ETHmuons_out), "SelectedMuons");
    iEvent.put(std::move(trans_muons_out), "SelectedTransientMuons");
}



DEFINE_FWK_MODULE(MuonTriggerSelector);
