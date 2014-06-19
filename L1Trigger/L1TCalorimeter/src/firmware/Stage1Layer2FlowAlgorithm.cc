///
/// \class l1t::Stage1Layer2FlowAlgorithm
///
/// \authors: Maxime Guilbaud 
///
/// Description: Flow Algorithm HI

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2HFRingSumAlgorithmImp.h"
#include "L1Trigger/L1TCalorimeter/interface/PUSubtractionMethods.h"
#include "L1Trigger/L1TCalorimeter/interface/legacyGtHelper.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2EtSumAlgorithmImp.h"

l1t::Stage1Layer2FlowAlgorithm::Stage1Layer2FlowAlgorithm(CaloParamsStage1* params) : params_(params)
{
  //now do what ever initialization is needed
  for(unsigned int i = 0; i < L1CaloRegionDetId::N_PHI; i++) {
    sinPhi.push_back(sin(2. * 3.1415927 * i * 1.0 / L1CaloRegionDetId::N_PHI));
    cosPhi.push_back(cos(2. * 3.1415927 * i * 1.0 / L1CaloRegionDetId::N_PHI));
  }
}


l1t::Stage1Layer2FlowAlgorithm::~Stage1Layer2FlowAlgorithm() {


}


void l1t::Stage1Layer2FlowAlgorithm::processEvent(const std::vector<l1t::CaloRegion> & regions,
						  const std::vector<l1t::CaloEmCand> & EMCands,
						  std::vector<l1t::HFRingSum> * counts) {


  double q2x = 0;
  double q2y = 0;
  double regionET=0.;  
    
  for(std::vector<CaloRegion>::const_iterator region = regions.begin(); region != regions.end(); region++) {
   
    int ieta=region->hwEta();    
    if (ieta > 3 && ieta < 18) {
      continue;
    }

    int iphi=region->hwPhi();    
    regionET=region->hwPt();
    
    q2x+= regionET * cosPhi[iphi];
    q2y+= regionET * sinPhi[iphi];
  }

  double HFq2 = q2x*q2x+q2y*q2y;
  //double psi2 = 0.5 * atan(q2y/q2x);
  //std::cout << HFq2 << std::endl;
  
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > etLorentz(0,0,0,0);

  // convert back to hardware ET
  l1t::HFRingSum etTot (*&etLorentz,HFRingSum::HFRingSumType::V2,HFq2,0,0,0);

  std::vector<l1t::HFRingSum> *preGtHFRingSums = new std::vector<l1t::HFRingSum>();
  preGtHFRingSums->push_back(etTot);

  // All algorithms
  //EtSumToGtScales(params_, preGtHFRingSums, counts);

  delete preGtHFRingSums;

}
