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
#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"

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


  edm::InputTag centralityTag_;  
  //edm::EDGetTokenT<reco::Centrality> centralityToken;
  edm::Handle<reco::Centrality> centrality_;
  CentralityProvider * centProvider;
  edm::InputTag vertexTag_;
  //edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken;
  edm::Handle<std::vector<reco::Vertex>> vertex_;

  edm::InputTag inputPlanesTag_;


  edm::InputTag trackTag_;
  //edm::EDGetTokenT<reco::TrackCollection> trackToken;
  edm::Handle<reco::TrackCollection> trackCollection_;

  int FlatOrder_;
  double caloCentRef_;
  double caloCentRefWidth_;
  int caloCentRefMinBin_;
  int caloCentRefMaxBin_;
  int NumFlatBins_;
  int CentBinCompression_;
  int HFEtScale_;
  int Noffmin_;
  int Noffmax_;
  double ntrkval;
  unsigned int runno_;
  bool FirstEvent;
  HiEvtPlaneFlatten * flat[NumEPNames];
  RPFlatParams * rpFlat;
  int nRP;
  bool useOffsetPsi_;
  int Hbins;
  int Obins;
  bool UseEtHF;

  int getNoff(const edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
    int Noff = 0;
    using namespace edm;
    using namespace reco;
  
    iEvent.getByLabel(vertexTag_, vertex_);
    const VertexCollection * recoVertices = vertex_.product();
  
    int primaryvtx = 0;
    math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
    double vxError = (*recoVertices)[primaryvtx].xError();
    double vyError = (*recoVertices)[primaryvtx].yError();
    double vzError = (*recoVertices)[primaryvtx].zError();
    iEvent.getByLabel(trackTag_,trackCollection_);
    for(TrackCollection::const_iterator itTrack = trackCollection_->begin();
	itTrack != trackCollection_->end();                      
	++itTrack) {    
      if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
      if ( itTrack->charge() == 0 ) continue;
      if ( itTrack->pt() < 0.4 ) continue;
      double d0 = -1.* itTrack->dxy(v1);
      double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
      double dz=itTrack->dz(v1);
      double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
      if ( fabs(itTrack->eta()) > 2.4 ) continue;
      if ( fabs( dz/dzerror ) > 3. ) continue;
      if ( fabs( d0/derror ) > 3. ) continue;
      if ( itTrack->ptError()/itTrack->pt() > 0.1 ) continue;
      Noff++;
    }
    return Noff;
  }


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
HiEvtPlaneFlatProducer::HiEvtPlaneFlatProducer(const edm::ParameterSet& iConfig):runno_(0)
{
  UseEtHF = kFALSE;
  centralityTag_ = iConfig.getParameter<edm::InputTag>("centralityTag_");
  vertexTag_  = iConfig.getParameter<edm::InputTag>("vertexTag_");
  trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag_");
  inputPlanesTag_ = iConfig.getParameter<edm::InputTag>("inputPlanesTag_");
  FlatOrder_ = iConfig.getUntrackedParameter<int>("FlatOrder_", 9);
  NumFlatBins_ = iConfig.getUntrackedParameter<int>("NumFlatBins_",40);
  CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",5);
  caloCentRef_ = iConfig.getUntrackedParameter<double>("caloCentRef_",80.);
  caloCentRefWidth_ = iConfig.getUntrackedParameter<double>("caloCentRefWidth_",5.);
  HFEtScale_ = iConfig.getUntrackedParameter<int>("HFEtScale_",3800);
  Noffmin_ = iConfig.getUntrackedParameter<int>("Noffmin_", 0);
  Noffmax_ = iConfig.getUntrackedParameter<int>("Noffmax_", 50000);	
  useOffsetPsi_ = iConfig.getUntrackedParameter<bool>("useOffsetPsi_",true);
  FirstEvent = kTRUE;
   //register your products
  produces<reco::EvtPlaneCollection>();
   //now do what ever other initialization is needed
  for(int i = 0; i<NumEPNames; i++) {
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->Init(FlatOrder_,NumFlatBins_,HFEtScale_,EPNames[i],EPOrder[i]);
  }
  Hbins = flat[0]->GetHBins();
  Obins = flat[0]->GetOBins();
  centProvider = 0;
  
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

  bool newrun = false;
  if(runno_ != iEvent.id().run()) newrun = true;
  runno_ = iEvent.id().run();
  //
  //Get Flattening Parameters
  //
  if(Noffmin_>=0) {
    int Noff = getNoff( iEvent, iSetup);
    ntrkval = Noff;
    if ( (Noff < Noffmin_) or (Noff >= Noffmax_) ) {
      return;
    }
  }
  if(FirstEvent || newrun) {
    FirstEvent = false;
    newrun = false;
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    LoadEPDB * db = new LoadEPDB(flatparmsDB_,flat);
    if(!db->IsSuccess()) return;
  }
  //
  //Get Centrality
  //

  int bin = 0;
  try{iEvent.getByLabel(centralityTag_, centrality_);} catch(...){;}
  if(centrality_.isValid()) { 
    double hfetval = centrality_->EtHFtowerSum();
    bin = flat[0]->GetHFbin(hfetval);
  }
  
  if(!centProvider) {
    for(int i = 0; i<NumEPNames; i++) flat[i]->SetCaloCentRefBins(-1,-1);
    centProvider = new CentralityProvider(iSetup);
    int minbin = (caloCentRef_-caloCentRefWidth_/2.)*centProvider->getNbins()/100.;
    int maxbin = (caloCentRef_+caloCentRefWidth_/2.)*centProvider->getNbins()/100.;
    minbin/=CentBinCompression_;
    maxbin/=CentBinCompression_;
    if(minbin>0 && maxbin>=minbin) {
      for(int i = 0; i<NumEPNames; i++) {
	if(EPDet[i]==HF || EPDet[i]==Castor) flat[i]->SetCaloCentRefBins(minbin,maxbin);
      }
    }
  }
  centProvider->newEvent(iEvent,iSetup);
  int cbin = centProvider->getBin();
  if(!UseEtHF) bin = cbin/CentBinCompression_; 

  //
  //Get Vertex
  //
  int vs_sell;   // vertex collection size
  float vzr_sell;
  iEvent.getByLabel(vertexTag_,vertex_);
  const reco::VertexCollection * vertices3 = vertex_.product();
  vs_sell = vertices3->size();
  if(vs_sell>0) {
    vzr_sell = vertices3->begin()->z();
  } else
    vzr_sell = -999.9;
  
  //
  //Get Event Planes
  //
  
  Handle<reco::EvtPlaneCollection> evtPlanes;
  iEvent.getByLabel(inputPlanesTag_,evtPlanes);
  
  if(!evtPlanes.isValid()){
    //    cout << "Error! Can't get hiEvtPlane product!" << endl;
    return ;
  }

  std::auto_ptr<EvtPlaneCollection> evtplaneOutput(new EvtPlaneCollection);
  EvtPlane * ep[NumEPNames];
  for(int i = 0; i<NumEPNames; i++) {
    ep[i]=0;
  }
  int i = 0;
  for (EvtPlaneCollection::const_iterator rp = evtPlanes->begin();rp !=evtPlanes->end(); rp++) {
	double angorig = rp->angle();
	double s = rp->sumSin();
	double c = rp->sumCos();
	double w = rp->sumw();
	uint m = rp->mult();
	
	double psiOffset = angorig;
	if(useOffsetPsi_) psiOffset = flat[i]->GetOffsetPsi(s,c,w,m,vzr_sell,bin);
	double psiFlat = flat[i]->GetFlatPsi(psiOffset,vzr_sell,bin);
	;
	ep[i]= new EvtPlane(i, 2, psiFlat, flat[i]->sumSin(), flat[i]->sumCos(),rp->sumw(), rp->sumw2(), rp->sumPtOrEt(), rp->sumPtOrEt2(),  rp->mult());
	ep[i]->AddLevel(0,rp->angle(), rp->sumSin(), rp->sumCos());
	ep[i]->AddLevel(3,0., rp->sumSin(3), rp->sumCos(3));
	if(useOffsetPsi_) ep[i]->AddLevel(1, psiOffset, s, c);
	++i;
    
  }
  
  for(int i = 0; i< NumEPNames; i++) {
    if(ep[i]!=0) evtplaneOutput->push_back(*ep[i]);
    
  }
  iEvent.put(evtplaneOutput);
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
