#ifndef BPHNANO_CLASSES_H
#define BPHNANO_CLASSES_H

#include "DataFormats/Common/interface/Wrapper.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "PhysicsTools/BParkingNano/src/KalmanVertexFitter.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"

namespace {
  struct dictionary {
      std::vector<reco::TransientTrack> ttv;
      edm::Wrapper<std::vector<reco::TransientTrack> > wttv; 
      HNLKalmanVertexFitter KinVtx_;
  };
}

#endif  
