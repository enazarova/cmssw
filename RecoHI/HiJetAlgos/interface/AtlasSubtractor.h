#ifndef __AtlasSubtractor_h_
#define __AtlasSubtractor_h_

#include "RecoJets/JetProducers/interface/PileUpSubtractor.h"

class AtlasSubtractor : public PileUpSubtractor {
 public:
 AtlasSubtractor(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC) : PileUpSubtractor(iConfig, std::move(iC)),
     sumRecHits_(iConfig.getParameter<bool>("sumRecHits")),
     dropZeroTowers_(iConfig.getUntrackedParameter<bool>("dropZeroTowers",true))
       {;}
    virtual void offsetCorrectJets();
    void rescaleRMS(double s);
    double getEt(const reco::CandidatePtr & in) const;
    double getEta(const reco::CandidatePtr & in) const;
    virtual void calculatePedestal(std::vector<fastjet::PseudoJet> const & coll);
    virtual void subtractPedestal(std::vector<fastjet::PseudoJet> & coll);

    bool sumRecHits_;
    bool dropZeroTowers_;
    ~AtlasSubtractor(){;}
    
};


#endif
