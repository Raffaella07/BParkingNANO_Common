#ifndef PhysicsTools_BParkingNano_ETHMuon_h
#define PhysicsTools_BParkingNano_ETHMuon_h

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVect;

namespace pat {
  class ETHMuon : public PATObject<reco::CompositeCandidate>{
    public:
      // default constructor
      ETHMuon() {}

      // constructor for slimmedMuons
      ETHMuon(const pat::Muon& muon)
      {
        isSlimmedMuon_ = true;
        isDSAMuon_ = false;
        pt_ = muon.pt(); 
        eta_ = muon.eta(); 
        phi_ = muon.phi();
        charge_ = muon.charge();
        mass_ = muon.mass();
        vx_ = muon.vx();
        vy_ = muon.vy();
        vz_ = muon.vz();
        dz_ = muon.dB(muon.PVDZ);
        dzS_ = muon.dB(muon.PVDZ)/muon.edB(muon.PVDZ);
        dxy_ = muon.dB(muon.PV2D);
        dxyS_= muon.dB(muon.PV2D)/muon.edB(muon.PV2D);
        ip3d_ = muon.dB(muon.PV3D);
        sip3d_ = muon.dB(muon.PV3D)/muon.edB(muon.PV3D); 
        pfiso03_all_ = muon.pfIsolationR03().sumChargedHadronPt + std::max(muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - muon.pfIsolationR03().sumPUPt/2,static_cast<float>(0.0));
        pfiso03Rel_all_ = (muon.pfIsolationR03().sumChargedHadronPt + std::max(muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - muon.pfIsolationR03().sumPUPt/2,static_cast<float>(0.0)))/muon.pt();
        pfiso04_all_ = muon.pfIsolationR04().sumChargedHadronPt + std::max(muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - muon.pfIsolationR04().sumPUPt/2,static_cast<float>(0.0));
        pfiso04Rel_all_ = (muon.pfIsolationR04().sumChargedHadronPt + std::max(muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - muon.pfIsolationR04().sumPUPt/2,static_cast<float>(0.0)))/muon.pt();
        isPFSlimmedMuon_ = muon.isPFMuon();
        isGlobalSlimmedMuon_ = muon.isGlobalMuon();
        isTrackerSlimmedMuon_ = muon.isTrackerMuon();
        isGlobalOrTrackerSlimmedMuon_ = muon.isGlobalMuon() || muon.isTrackerMuon();
        isGlobalNotTrackerSlimmedMuon_ = muon.isGlobalMuon() && !muon.isTrackerMuon();
        isTrackerNotGlobalSlimmedMuon_ = !muon.isGlobalMuon() && muon.isTrackerMuon();
        looseId_ = muon.isLooseMuon();
        mediumId_ = muon.reco::Muon::passed(reco::Muon::CutBasedIdMedium) ? 1 : 0;
        tightId_ = muon.reco::Muon::passed(reco::Muon::CutBasedIdTight) ? 1 : 0;
        softId_ = muon.reco::Muon::passed(reco::Muon::SoftCutBasedId) ? 1 : 0;
        pfIsoId_ = (muon.reco::Muon::passed(reco::Muon::PFIsoVeryLoose) + muon.reco::Muon::passed(reco::Muon::PFIsoLoose) + muon.reco::Muon::passed(reco::Muon::PFIsoMedium) + muon.reco::Muon::passed(reco::Muon::PFIsoTight) + muon.reco::Muon::passed(reco::Muon::PFIsoVeryTight) + muon.reco::Muon::passed(reco::Muon::PFIsoVeryVeryTight)) ? 1 : 0;
        tkIsoId_ = muon.reco::Muon::passed(reco::Muon::TkIsoTight) || muon.reco::Muon::passed(reco::Muon::TkIsoLoose) ? 1 : 0;
        triggerIdLoose_ = muon.reco::Muon::passed(reco::Muon::TriggerIdLoose);
        inTimeMuon_ = muon.reco::Muon::passed(reco::Muon::InTimeMuon);
        segmentCompatibility_ = muon.segmentCompatibility();
        caloCompatibility_ = muon.caloCompatibility();
        validHitFraction_ = muon.isGlobalMuon() || muon.isTrackerMuon() ? muon.innerTrack()->validFraction(): -1.;
        kinkFinderChi2_ = muon.combinedQuality().trkKink;
        globalNormalisedChi2_ = muon.isGlobalMuon() ? muon.globalTrack()->normalizedChi2(): -1.;
        localPositionChi2_ = muon.combinedQuality().chi2LocalPosition;
        trackerHighPurityFlag_ = muon.isGlobalMuon() || muon.isTrackerMuon() ? muon.innerTrack()->quality(reco::TrackBase::highPurity): -1;
        numberOfValidMuonHits_ = muon.isGlobalMuon() ? muon.globalTrack()->hitPattern().numberOfValidMuonHits(): -1;
        numberOfValidPixelHits_ = muon.isGlobalMuon() || muon.isTrackerMuon() ? muon.innerTrack()->hitPattern().numberOfValidPixelHits(): -1;
        numberOfTrackerLayers_ = muon.isGlobalMuon() || muon.isTrackerMuon() ? muon.innerTrack()->hitPattern().trackerLayersWithMeasurement(): -1;
        numberOfPixelLayers_ = muon.isGlobalMuon() || muon.isTrackerMuon() ? muon.innerTrack()->hitPattern().pixelLayersWithMeasurement(): -1;
        numberOfStations_ = muon.numberOfMatchedStations();
      }

      // constructor for DSAMuons
      ETHMuon(const reco::Track& muon)
      {
        isSlimmedMuon_ = false;
        isDSAMuon_ = true;
        pt_ = muon.pt(); 
        eta_ = muon.eta(); 
        phi_ = muon.phi();
        charge_ = muon.charge();
        mass_ = 0.105658;
        vx_ = muon.vx();
        vy_ = muon.vy();
        vz_ = muon.vz();
        dz_ = muon.dz(); //0.001; //99.; // needs to be studied, see https://cmssdt.cern.ch/lxr/source/DataFormats/TrackReco/interface/TrackBase.h
        dzS_ = 1.; //99.; // needs to be studied, see https://cmssdt.cern.ch/lxr/source/DataFormats/TrackReco/interface/TrackBase.h
        dxy_ = 0.001; // needs to be studied, see https://cmssdt.cern.ch/lxr/source/DataFormats/TrackReco/interface/TrackBase.h
        dxyS_ = 1.; //99.; // needs to be studied, see https://cmssdt.cern.ch/lxr/source/DataFormats/TrackReco/interface/TrackBase.h
        ip3d_ = 99; // to modify
        sip3d_ = 99; // to modify
        pfiso03_all_ = -99.; // to modify
        pfiso03Rel_all_ = -99.; // to modify
        pfiso04_all_ = -99.; // to modify
        pfiso04Rel_all_ = -99.; // to modify
        isPFSlimmedMuon_ = -1;
        isGlobalSlimmedMuon_ = -1;
        isTrackerSlimmedMuon_ = -1;
        isGlobalOrTrackerSlimmedMuon_ = -1;
        isGlobalNotTrackerSlimmedMuon_ = -1;
        isTrackerNotGlobalSlimmedMuon_ = -1;
        looseId_ = -1;
        mediumId_ = -1;
        tightId_ = -1;
        softId_ = -1;
        pfIsoId_ = -1;
        tkIsoId_ = -1;
        triggerIdLoose_ = -1;
        inTimeMuon_ = -1;
        segmentCompatibility_ = -1.;
        caloCompatibility_ = -1.;
        validHitFraction_ = muon.validFraction();
        kinkFinderChi2_ =  -1.;
        globalNormalisedChi2_ =  muon.normalizedChi2();
        localPositionChi2_ = -1.;
        trackerHighPurityFlag_ =  muon.quality(reco::TrackBase::highPurity);
        numberOfValidMuonHits_ =  muon.hitPattern().numberOfValidMuonHits();
        numberOfValidPixelHits_ = muon.hitPattern().numberOfValidPixelHits();
        numberOfTrackerLayers_ =  muon.hitPattern().trackerLayersWithMeasurement();
        numberOfPixelLayers_ = muon.hitPattern().pixelLayersWithMeasurement();
        numberOfStations_ = -1;
      }

      // destructor
      ~ETHMuon(){};

      // methods
      bool isSlimmedMuon() const { return (*this).isSlimmedMuon_; }
      bool isDSAMuon() const { return (*this).isDSAMuon_; }
      double pt() const override { return (*this).pt_; }
      double eta() const override { return (*this).eta_; }
      double phi() const override { return (*this).phi_; }
      int charge() const override { return (*this).charge_; }
      double mass() const override { return (*this).mass_; } 
      double vx() const override { return (*this).vx_; }
      double vy() const override { return (*this).vy_; }
      double vz() const override { return (*this).vz_; }
      double dz() const { return (*this).dz_; }
      double dzS() const { return (*this).dzS_; }
      double dxy() const { return (*this).dxy_; }
      double dxyS() const { return (*this).dxyS_; }
      double ip3d() const { return (*this).ip3d_; }
      double sip3d() const { return (*this).sip3d_; }
      double pfiso03_all() const { return (*this).pfiso03_all_; }
      double pfiso03Rel_all() const { return (*this).pfiso03Rel_all_; }
      double pfiso04_all() const { return (*this).pfiso04_all_; }
      double pfiso04Rel_all() const { return (*this).pfiso04Rel_all_; }
      int isPFSlimmedMuon() const { return (*this).isPFSlimmedMuon_; }
      int isGlobalSlimmedMuon() const { return (*this).isGlobalSlimmedMuon_; }
      int isTrackerSlimmedMuon() const { return (*this).isTrackerSlimmedMuon_; }
      int isGlobalOrTrackerSlimmedMuon() const { return (*this).isGlobalOrTrackerSlimmedMuon_; }
      int isGlobalNotTrackerSlimmedMuon() const { return (*this).isGlobalNotTrackerSlimmedMuon_; }
      int isTrackerNotGlobalSlimmedMuon() const { return (*this).isTrackerNotGlobalSlimmedMuon_; }
      int looseId() const { return (*this).looseId_; }
      int mediumId() const { return (*this).mediumId_; }
      int tightId() const { return (*this).tightId_; }
      int softId() const { return (*this).softId_; }
      int pfIsoId() const { return (*this).pfIsoId_; }
      int tkIsoId() const { return (*this).tkIsoId_; }
      int triggerIdLoose() const { return (*this).triggerIdLoose_; }
      int inTimeMuon() const { return (*this).inTimeMuon_; }
      double segmentCompatibility() const { return (*this).segmentCompatibility_; }
      double caloCompatibility() const { return (*this).caloCompatibility_; }
      double validHitFraction() const { return (*this).validHitFraction_; }
      double kinkFinderChi2() const { return (*this).kinkFinderChi2_; }
      double globalNormalisedChi2() const { return (*this).globalNormalisedChi2_; }
      double localPositionChi2() const { return (*this).localPositionChi2_; }
      int trackerHighPurityFlag() const { return (*this).trackerHighPurityFlag_; }
      int numberOfValidMuonHits() const { return (*this).numberOfValidMuonHits_; }
      int numberOfValidPixelHits() const { return (*this).numberOfValidPixelHits_; }
      int numberOfTrackerLayers() const { return (*this).numberOfTrackerLayers_; }
      int numberOfPixelLayers() const { return (*this).numberOfPixelLayers_; }
      int numberOfStations() const { return (*this).numberOfStations_; }

      LorentzVect P4() const{
        LorentzVect p4((*this).pt(), (*this).eta(), (*this).phi(), (*this).mass());
        return p4;
      };

      const PolarLorentzVector & polarP4() const override{
        static PolarLorentzVector p4((*this).pt(), (*this).eta(), (*this).phi(), (*this).mass());
        return p4;
      };

    private:
      bool isSlimmedMuon_;
      bool isDSAMuon_;
      double pt_;
      double eta_;
      double phi_;
      int charge_;
      double mass_;
      double vx_;
      double vy_;
      double vz_;
      double dz_;
      double dzS_;
      double dxy_;
      double dxyS_;
      double ip3d_;
      double sip3d_;
      double pfiso03_all_;
      double pfiso03Rel_all_;
      double pfiso04_all_;
      double pfiso04Rel_all_;
      int isPFSlimmedMuon_;
      int isGlobalSlimmedMuon_;
      int isTrackerSlimmedMuon_;
      int isGlobalOrTrackerSlimmedMuon_;
      int isGlobalNotTrackerSlimmedMuon_;
      int isTrackerNotGlobalSlimmedMuon_;
      int looseId_;
      int mediumId_;
      int tightId_;
      int softId_;
      int pfIsoId_;
      int tkIsoId_;
      int triggerIdLoose_;
      int inTimeMuon_;
      double segmentCompatibility_;
      double caloCompatibility_;
      double validHitFraction_;
      double kinkFinderChi2_;
      double globalNormalisedChi2_;
      double localPositionChi2_;
      int trackerHighPurityFlag_;
      int numberOfValidMuonHits_;
      int numberOfValidPixelHits_;
      int numberOfTrackerLayers_;
      int numberOfPixelLayers_;
      int numberOfStations_;
  };
}


#endif
