// -*- C++ -*-
//
// Package:    EvtPlaneProducer
// Class:      EvtPlaneProducer
// 
/**\class EvtPlaneProducer EvtPlaneProducer.cc RecoHI/EvtPlaneProducer/src/EvtPlaneProducer.cc
   SS
Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Sergey Petrushanko
//         Created:  Fri Jul 11 10:05:00 2008
// $Id: EvtPlaneProducer.cc,v 1.18 2011/10/07 09:41:29 yilmaz Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <time.h>
#include "TMath.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <cstdlib>
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"

using namespace std;
using namespace hi;

//
// class decleration
//

class EvtPlaneProducer : public edm::EDProducer {
public:
  explicit EvtPlaneProducer(const edm::ParameterSet&);
  ~EvtPlaneProducer();
  
private:
  class GenPlane {
  public: 
    GenPlane(string name,double etaminval1,double etamaxval1,double etaminval2,double etamaxval2,int orderval){
      epname=name;
      etamin1=etaminval1;
      etamax1=etamaxval1;
      etamin2=etaminval2;
      etamax2=etamaxval2;
      sumsin=0;
      sumcos=0;
      mult = 0;
      order = (double) orderval;
    }
    ~GenPlane(){;}
    void addParticle(double w, double s, double c, double eta) {
      if(fabs(w) < 0.0001) return;
      if((eta>=etamin1 && eta<etamax1) || 
	 (etamin2!= etamax2 && eta>=etamin2 && eta<etamax2 )) {
	sumsin+=w*s;
	sumcos+=w*c;
	sumw+=w;
	++mult;
      }
    }
    
    double getAngle(double &ang, double &sv, double &cv, double &w, uint &epmult){
      ang = -10;
      sv = 0;
      cv = 0;
      sv = sumsin;
      cv = sumcos;
      w = sumw;
      epmult = mult;
      double q = sv*sv+cv*cv;
      if(q>0) ang = atan2(sv,cv)/order;
      return ang;
    }
    void reset() {
      sumsin=0;
      sumcos=0;
      sumw = 0;
      mult = 0;
    }
  private:
    string epname;
    double etamin1;
    double etamax1;

    double etamin2;
    double etamax2;
    double sumsin;
    double sumcos;
    uint mult;
    double sumw;
    double order;
  };
  

  GenPlane *rp[NumEPNames];

  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------

  edm::InputTag centralityTag_;  
  edm::EDGetTokenT<reco::Centrality> centralityToken;
  edm::Handle<reco::Centrality> centrality_;

  edm::InputTag vertexTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken;
  edm::Handle<std::vector<reco::Vertex>> vertex_;

  edm::InputTag caloTag_;
  edm::EDGetTokenT<CaloTowerCollection> caloToken;
  edm::Handle<CaloTowerCollection> caloCollection_;

  edm::InputTag trackTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken;
  edm::Handle<reco::TrackCollection> trackCollection_;




  bool useECAL_;
  bool useHCAL_;
  bool useTrack_;
  bool loadDB_;
  double minet_;
  double maxet_;
  double effm_;
  double minpt_;
  double maxpt_;
  double minvtx_;
  double maxvtx_;
  double dzerr_;
  double chi2_;
  bool FirstEvent;
  int FlatOrder_;
  int NumFlatCentBins_;
  int CentBinCompression_;
  CentralityProvider * centProvider;
  HiEvtPlaneFlatten * flat[NumEPNames];
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EvtPlaneProducer::EvtPlaneProducer(const edm::ParameterSet& iConfig) {
  centProvider = 0;

  centralityTag_ = iConfig.getParameter<edm::InputTag>("hiCentrality_");
  centralityToken = consumes<reco::Centrality>(centralityTag_);

  vertexTag_  = iConfig.getParameter<edm::InputTag>("vtxCollection_");
  vertexToken = consumes<std::vector<reco::Vertex>>(vertexTag_);

  caloTag_ = iConfig.getParameter<edm::InputTag>("caloCollection_");
  caloToken = consumes<CaloTowerCollection>(caloTag_);

  trackTag_ = iConfig.getParameter<edm::InputTag>("trackCollection_");
  trackToken = consumes<reco::TrackCollection>(trackTag_);

  FlatOrder_ = iConfig.getUntrackedParameter<int>("FlatOrder_", 9);
  NumFlatCentBins_ = iConfig.getUntrackedParameter<int>("NumFlatCentBins_",50);
  CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",4);
  loadDB_ = iConfig.getUntrackedParameter<bool>("loadDB_",true);

  minet_ = iConfig.getUntrackedParameter<double>("minet_",0.5);
  maxet_ = iConfig.getUntrackedParameter<double>("maxet_",80.);
  minpt_ = iConfig.getUntrackedParameter<double>("minpt_",0.3);

  maxpt_ = iConfig.getUntrackedParameter<double>("maxpt_",2.6);
  minvtx_ = iConfig.getUntrackedParameter<double>("minvtx_",-25.);
  maxvtx_ = iConfig.getUntrackedParameter<double>("maxvtx_",25.);
  dzerr_ = iConfig.getUntrackedParameter<double>("dzerr_",10.);
  chi2_  = iConfig.getUntrackedParameter<double>("chi2_",40.);
  FirstEvent = kTRUE;
  produces<reco::EvtPlaneCollection>("recoLevel");
  for(int i = 0; i<NumEPNames; i++ ) {
    rp[i] = new GenPlane(EPNames[i].data(),EPEtaMin1[i],EPEtaMax1[i],EPEtaMin2[i],EPEtaMax2[i],EPOrder[i]);
  }
  for(int i = 0; i<NumEPNames; i++) {
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->Init(FlatOrder_,NumFlatCentBins_,CentBinCompression_,EPNames[i],EPOrder[i]);
  }
  cout<<"====================="<<endl;
  cout<<"EvtPlaneProducer: "<<endl;
  cout<<"  minet_:            "<<minet_<<endl;
  cout<<"  maxet_:            "<<maxet_<<endl;
  cout<<"  minpt_:            "<<minpt_<<endl;
  cout<<"  maxpt_:            "<<maxpt_<<endl; 
  cout<<"  minvtx_:           "<<minvtx_<<endl;
  cout<<"  maxvtx_:           "<<maxvtx_<<endl;
  cout<<"  dzerr_             "<<dzerr_<<endl;
  cout<<"  chi2_              "<<chi2_<<endl;
  cout<<"====================="<<endl;

}


EvtPlaneProducer::~EvtPlaneProducer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
EvtPlaneProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
 
  if(FirstEvent && loadDB_) {
    FirstEvent = kFALSE;
    //
    //Get flattening parameter file.  
    //
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    int flatTableSize = flatparmsDB_->m_table.size();
    for(int i = 0; i<flatTableSize; i++) {
      const RPFlatParams::EP* thisBin = &(flatparmsDB_->m_table[i]);
      for(int j = 0; j<NumEPNames; j++) {
	int indx = thisBin->RPNameIndx[j];
	int    Hbins = flat[j]->GetHBins();
	int    Obins = flat[j]->GetOBins();
	if(indx>=0) {
	  if(i<Hbins) {
	    flat[indx]->SetXDB(i, thisBin->x[j]);
	    flat[indx]->SetYDB(i, thisBin->y[j]);
	  } else if(i>=Hbins && i<Hbins+Obins) {
	    flat[indx]->SetXoffDB(i - Hbins, thisBin->x[j]);
	    flat[indx]->SetYoffDB(i - Hbins, thisBin->y[j]);
	  } else if (i>=Hbins+Obins) {
	    flat[indx]->SetPtDB(i - Hbins - Obins, thisBin->x[j]);
	    flat[indx]->SetPt2DB(i - Hbins - Obins, thisBin->y[j]);
	  }
	}
      }
    }
  } //First event

  //
  //Get Centrality
  //
  //if(!centProvider) centProvider = new CentralityProvider(iSetup);
  //centProvider->newEvent(iEvent,iSetup);
  //centProvider->raw();
  iEvent.getByToken(centralityToken, centrality_); 

  //int bin = centProvider->getBin();
  int bin = 0;
  int vs_sell;
  float vzr_sell;

  //
  //Get Vertex
  //
  iEvent.getByToken(vertexToken,vertex_);
  const reco::VertexCollection * vertices3 = vertex_.product();
  vs_sell = vertices3->size();
  if(vs_sell>0) {
    vzr_sell = vertices3->begin()->z();
  } else
    vzr_sell = -999.9;
  //
  for(int i = 0; i<NumEPNames; i++) rp[i]->reset();
  if(vzr_sell>minvtx_ && vzr_sell<maxvtx_) {

  

    //calorimetry part
    
    double tower_eta, tower_phi;
    double tower_energyet, tower_energyet_e, tower_energyet_h;

    iEvent.getByToken(caloToken,caloCollection_);
    
    if(caloCollection_.isValid()){
      for (CaloTowerCollection::const_iterator j = caloCollection_->begin();j !=caloCollection_->end(); j++) {   
	tower_eta        = j->eta();
	tower_phi        = j->phi();
	tower_energyet_e   = j->emEt();
	tower_energyet_h   = j->hadEt();
	tower_energyet     = tower_energyet_e + tower_energyet_h;
	if(tower_energyet<minet_) continue;
	if(tower_energyet>maxet_) continue;
	for(int i = 0; i<NumEPNames; i++) {
	  if(EPDet[i]==HF) {
	    double w = tower_energyet;
	    if(EPOrder[i]==1 ) {
	      if(MomConsWeight[i][0]=='y' && loadDB_ ) {
		w = flat[i]->GetW(tower_energyet, vzr_sell, bin);
	      }
	      if(tower_eta<0 ) w=-w;
	    }
	    rp[i]->addParticle(w,sin(EPOrder[i]*tower_phi),cos(EPOrder[i]*tower_phi),tower_eta);
	  }
	}
      } 
    }

    //Tracking part
    
    double track_eta;
    double track_phi;
    double track_pt;
   
    iEvent.getByToken(trackToken, trackCollection_);
    if(trackCollection_.isValid()){
      for(reco::TrackCollection::const_iterator j = trackCollection_->begin(); j != trackCollection_->end(); j++){	
	iEvent.getByToken(vertexToken, vertex_);	    
	math::XYZPoint vtxPoint(0.0,0.0,0.0);
	double vzErr =0.0, vxErr=0.0, vyErr=0.0;
	if(vertex_->size()>0) {
	  vtxPoint=vertex_->begin()->position();
	  vzErr=vertex_->begin()->zError();
	  vxErr=vertex_->begin()->xError();
	  vyErr=vertex_->begin()->yError();
	}
	bool accepted = true;
	bool isPixel = false;
	// determine if the track is a pixel track
	if ( j->numberOfValidHits() < 7 ) isPixel = true;
	
	// determine the vertex significance 
	double d0=0.0, dz=0.0, d0sigma=0.0, dzsigma=0.0;
	d0 = -1.*j->dxy(vtxPoint);
	dz = j->dz(vtxPoint);
	d0sigma = sqrt(j->d0Error()*j->d0Error()+vxErr*vyErr);
	dzsigma = sqrt(j->dzError()*j->dzError()+vzErr*vzErr);
	
	// cuts for pixel tracks
	if( isPixel ){
	  // dz significance cut 
	  if ( fabs(dz/dzsigma) > dzerr_ ) accepted = false;
	  // chi2/ndof cut 
	  if ( j->normalizedChi2() > chi2_ ) accepted = false;
	}   
	// cuts for full tracks
	if ( ! isPixel) {
	  // dz and d0 significance cuts 
	  if ( fabs(dz/dzsigma) > 3 ) accepted = false;
	  if ( fabs(d0/d0sigma) > 3 ) accepted = false;
	  // pt resolution cut
	  if ( j->ptError()/j->pt() > 0.1 ) accepted = false;
	  // number of valid hits cut
	  if ( j->numberOfValidHits() < 12 ) accepted = false;
	}
	if( accepted ) {
	  track_eta = j->eta();
	  track_phi = j->phi();
	  track_pt = j->pt();
	  if(track_pt<minpt_) continue;
	  if(track_pt>maxpt_) continue;
	  for(int i = 0; i<NumEPNames; i++) {
	    if(EPDet[i]==Tracker) {
	      double w = track_pt;
	      if(w>2.5) w=2.0;   //v2 starts decreasing above ~2.5 GeV/c
	      if(EPOrder[i]==1) {
		if(MomConsWeight[i][0]=='y' && loadDB_) {
		  w = flat[i]->GetW(track_pt, vzr_sell, bin);
		}
		if(track_eta<0) w=-w;
	      }
	      rp[i]->addParticle(w,sin(EPOrder[i]*track_phi),cos(EPOrder[i]*track_phi),track_eta);
	    }
	  }
	}  
      } //end for
    }
    
    std::auto_ptr<EvtPlaneCollection> evtplaneOutput(new EvtPlaneCollection);
    EvtPlane *ep[NumEPNames];
    
    double ang=-10;
    double sv = 0;
    double cv = 0;
    double wv = 0;
    uint epmult = 0;

    for(int i = 0; i<NumEPNames; i++) {
      rp[i]->getAngle(ang,sv,cv,wv,epmult);
      ep[i] = new EvtPlane(ang,sv,cv,wv,epmult,EPNames[i]);
    }

    for(int i = 0; i<NumEPNames; i++) {
      evtplaneOutput->push_back(*ep[i]);
    }  
    
    
    iEvent.put(evtplaneOutput, "recoLevel");
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
EvtPlaneProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EvtPlaneProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EvtPlaneProducer);
