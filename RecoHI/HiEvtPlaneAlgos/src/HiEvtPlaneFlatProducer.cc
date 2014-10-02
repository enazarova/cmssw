// -*- C++ -*-
//
// Package:    HiEvtPlaneFlatProducer
// Class:      HiEvtPlaneFlatProducer
// 
/**\class HiEvtPlaneFlatProducer HiEvtPlaneFlatProducer.cc HiEvtPlaneFlatten/HiEvtPlaneFlatProducer/src/HiEvtPlaneFlatProducer.cc


 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Stephen Sanders
//         Created:  Sat Jun 26 16:04:04 EDT 2010
// $Id: HiEvtPlaneFlatProducer.cc,v 1.12 2011/12/04 05:13:03 ssanders Exp $
//
//
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"

#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/UpdatablePSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "TList.h"
#include "TString.h"
#include <time.h>
#include <cstdlib>

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"

using namespace std;
using namespace hi;

#include <vector>
using std::vector;


//
// class declaration
//

class HiEvtPlaneFlatProducer : public edm::EDProducer {
   public:
      explicit HiEvtPlaneFlatProducer(const edm::ParameterSet&);
      ~HiEvtPlaneFlatProducer();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

  edm::InputTag vtxCollection_;
  edm::InputTag inputPlanes_;
  CentralityProvider * centProvider;

  bool FirstEvent;
  HiEvtPlaneFlatten * flat[NumEPNames];
  RPFlatParams * rpFlat;
  int nRP;
  bool storeNames_;
  bool useOffsetPsi_;
  int Hbins;
  int Obins;
};

//
// constants, enums and typedefs
//
typedef std::vector<TrackingParticle>                   TrackingParticleCollection;
typedef TrackingParticleRefVector::iterator               tp_iterator;


//
// static data member definitions
//

//
// constructors and destructor
//
HiEvtPlaneFlatProducer::HiEvtPlaneFlatProducer(const edm::ParameterSet& iConfig)
{
  centProvider = 0;
  vtxCollection_  = iConfig.getParameter<edm::InputTag>("vtxCollection_");
  inputPlanes_ = iConfig.getParameter<edm::InputTag>("inputPlanes_");
  storeNames_ = 1;
  useOffsetPsi_ = iConfig.getUntrackedParameter<bool>("useOffsetPsi_",true);
  FirstEvent = kTRUE;
   //register your products
  produces<reco::EvtPlaneCollection>();
   //now do what ever other initialization is needed
  Int_t FlatOrder = 9;
  for(int i = 0; i<NumEPNames; i++) {
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->Init(FlatOrder,NumFlatCentBins,CentBinCompression,EPNames[i],EPOrder[i]);
  }
  Hbins = flat[0]->GetHBins();
  Obins = flat[0]->GetOBins();

  
}


HiEvtPlaneFlatProducer::~HiEvtPlaneFlatProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HiEvtPlaneFlatProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  //
  //Get Flattening Parameters
  //
  if(FirstEvent) {
    FirstEvent = kFALSE;
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    int flatTableSize = flatparmsDB_->m_table.size();
    for(int i = 0; i<flatTableSize; i++) {
      const RPFlatParams::EP* thisBin = &(flatparmsDB_->m_table[i]);
      for(int j = 0; j<NumEPNames; j++) {
	int indx = thisBin->RPNameIndx[j];
	if(indx>=0) {
	    if(i<Hbins) {
	      flat[indx]->SetXDB(i, thisBin->x[j]);
	      flat[indx]->SetYDB(i, thisBin->y[j]);
	    } else if(i>=Hbins && i<Hbins+Obins) {
	      flat[indx]->SetXoffDB(i - Hbins, thisBin->x[j]);
	      flat[indx]->SetYoffDB(i - Hbins, thisBin->y[j]);
	    } else if (i>=Hbins+Obins) {
	      flat[indx]->SetPtDB(i - Hbins- Obins, thisBin->x[j]);
	      flat[indx]->SetPt2DB(i - Hbins- Obins, thisBin->y[j]);
	    }
	}
      }
    }
    
  }
  //
  //Get Centrality
  //
  if(!centProvider) centProvider = new CentralityProvider(iSetup);
  centProvider->newEvent(iEvent,iSetup);
  centProvider->raw();
  int bin = centProvider->getBin();
  if(centProvider->getNbins()<=100) bin=2*bin;
 
  //
  //Get Vertex
  //
  int vs_sell;   // vertex collection size
  float vzr_sell;
  edm::Handle<reco::VertexCollection> vertexCollection3;
  iEvent.getByLabel(vtxCollection_,vertexCollection3);
  const reco::VertexCollection * vertices3 = vertexCollection3.product();
  vs_sell = vertices3->size();
  if(vs_sell>0) {
    vzr_sell = vertices3->begin()->z();
  } else
    vzr_sell = -999.9;
  
  //
  //Get Event Planes
  //
  
  Handle<reco::EvtPlaneCollection> evtPlanes;
  iEvent.getByLabel(inputPlanes_,evtPlanes);
  
  if(!evtPlanes.isValid()){
    //    cout << "Error! Can't get hiEvtPlane product!" << endl;
    return ;
  }

  std::auto_ptr<EvtPlaneCollection> evtplaneOutput(new EvtPlaneCollection);
  EvtPlane * ep[NumEPNames];
  for(int i = 0; i<NumEPNames; i++) {
    ep[i]=0;
  }
  for (EvtPlaneCollection::const_iterator rp = evtPlanes->begin();rp !=evtPlanes->end(); rp++) {
      string baseName = rp->label();
      for(int i = 0; i< NumEPNames; i++) {
	if(EPNames[i].compare(baseName)==0) {
	  double angorig = rp->angle();
	  double c = rp->sumCos();
	  double s = rp->sumSin();
	 
	  double psiOffset = angorig;
	  if(useOffsetPsi_) psiOffset = flat[i]->GetOffsetPsi(s,c,vzr_sell,bin);
	  double psiFlat = flat[i]->GetFlatPsi(psiOffset,vzr_sell,bin);
;
	  if(EPNames[i].compare(rp->label())==0) {	    
	    if(storeNames_) ep[i]= new EvtPlane(psiFlat, rp->sumSin(), rp->sumCos(),rp->mult(),rp->label().data());
	    else ep[i]= new EvtPlane(psiFlat, rp->sumSin(), rp->sumCos(),rp->mult(),"");
	  } 
	}
      }
  }
  
  for(int i = 0; i< NumEPNames; i++) {
    if(ep[i]!=0) evtplaneOutput->push_back(*ep[i]);
    
  }
  iEvent.put(evtplaneOutput);
  //storeNames_ = 0;  
}

// ------------ method called once each job just before starting event loop  ------------
void 
HiEvtPlaneFlatProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiEvtPlaneFlatProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiEvtPlaneFlatProducer);
