
// -*- C++ -*-
//
// Package:    CheckFlattening
// Class:      CheckFlattening
// 

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/HeavyIon.h"
#include "HepMC/SimpleVector.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"


#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <time.h>
#include <cstdlib>
	
#include <vector>
using std::vector;
using std::rand;
using namespace std;
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
using namespace hi;

//
// class declaration
//

class CheckFlattening : public edm::EDAnalyzer {
public:
  explicit CheckFlattening(const edm::ParameterSet&);
  ~CheckFlattening();
  
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  // ----------member data ---------------------------
  int eporder_;
  CentralityProvider * centProvider;

  edm::InputTag centralityTag_;  
  edm::EDGetTokenT<reco::Centrality> centralityToken;
  edm::Handle<reco::Centrality> centrality_;

  edm::InputTag vertexTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken;
  edm::Handle<std::vector<reco::Vertex>> vertex_;


  edm::InputTag inputPlanesTag_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> inputPlanesToken;
  edm::Handle<reco::EvtPlaneCollection> inputPlanes_;

  edm::InputTag inputPlanesFlatTag_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> inputPlanesFlatToken;
  edm::Handle<reco::EvtPlaneCollection> inputPlanesFlat_;


  edm::Service<TFileService> fs;
  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  TH1D * hcent;
  TH1D * hcentbins;
  double centval;
  double vtx;
  Double_t epang[NumEPNames];
  Double_t epang_orig[NumEPNames];
  Double_t epsin[NumEPNames];
  Double_t epcos[NumEPNames];
  Double_t epw[NumEPNames];
  Double_t epqx[NumEPNames];
  Double_t epqy[NumEPNames];
  Double_t epq[NumEPNames];
  Double_t epmult[NumEPNames];
  unsigned int runno_;

  int nEtaBins;
  TH1I * hrun;
  string rpnames[NumEPNames];
  TH1D * PsiRaw[50];
  TH1D * Psi[50];
  TTree * tree;
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
CheckFlattening::CheckFlattening(const edm::ParameterSet& iConfig):runno_(0)
  
{
  centProvider = 0;
  centralityTag_ = iConfig.getParameter<edm::InputTag>("hiCentrality_");
  centralityToken = consumes<reco::Centrality>(centralityTag_);

  vertexTag_  = iConfig.getParameter<edm::InputTag>("vtxCollection_");
  vertexToken = consumes<std::vector<reco::Vertex>>(vertexTag_);

  inputPlanesTag_ = iConfig.getParameter<edm::InputTag>("inputPlanes_");
  inputPlanesToken = consumes<reco::EvtPlaneCollection>(inputPlanesTag_);

  inputPlanesFlatTag_ = iConfig.getParameter<edm::InputTag>("inputPlanesFlat_");
  inputPlanesFlatToken = consumes<reco::EvtPlaneCollection>(inputPlanesFlatTag_);

  
  hcent = fs->make<TH1D>("cent","cent",220,-10,110);
  hcentbins = fs->make<TH1D>("centbins","centbins",201,0,200);
  hrun = fs->make<TH1I>("runs","runs",100000,150001,250000);

  TString epnames = EPNames[0].data();
  epnames = epnames+"/D";
  for(int i = 0; i<NumEPNames; i++) {
    if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/D";
    TFileDirectory subdir = fs->mkdir(Form("%s",EPNames[i].data()));
    Double_t psirange = 4;
    if(EPOrder[i]==2 ) psirange = 2;
    if(EPOrder[i]==3 ) psirange = 1.5;
    if(EPOrder[i]==4 ) psirange = 1;
    if(EPOrder[i]==5) psirange = 0.8;
    if(EPOrder[i]==6) psirange = 0.6;

    PsiRaw[i] = subdir.make<TH1D>(Form("PsiRaw_%s",EPNames[i].data()),Form("PsiRaw_%s",EPNames[i].data()),800,-psirange,psirange);
    Psi[i] = subdir.make<TH1D>(Form("Psi_%s",EPNames[i].data()),Form("Psi_%s",EPNames[i].data()),800,-psirange,psirange);

  }
  tree = fs->make<TTree>("tree","EP tree");
  tree->Branch("Cent",&centval,"cent/D");
  tree->Branch("Vtx",&vtx,"vtx/D");
  tree->Branch("EP",&epang, epnames.Data());
  tree->Branch("EP_orig",&epang_orig, epnames.Data());
  tree->Branch("Sin",     &epsin,      epnames.Data());
  tree->Branch("Cos",     &epcos,      epnames.Data());
  tree->Branch("Weight",  &epw,        epnames.Data());
  tree->Branch("qx",      &epqx,       epnames.Data());
  tree->Branch("qy",      &epqy,       epnames.Data());
  tree->Branch("Mult",    &epmult,     epnames.Data());
  tree->Branch("Run",&runno_,"run/i");
}


CheckFlattening::~CheckFlattening()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
CheckFlattening::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace HepMC;
  runno_ = iEvent.id().run();
  hrun->Fill(runno_);
  //
  //Get Centrality
  //
  int bin = 0;
  if(!centProvider) centProvider = new CentralityProvider(iSetup);
  centProvider->newEvent(iEvent,iSetup);
  centProvider->raw();
  bin = centProvider->getBin();
  centval = (100./centProvider->getNbins())*(bin+0.5);
  edm::Handle<reco::Centrality> cent;
  iEvent.getByToken(centralityToken, cent); 
  hcent->Fill(centval);
  hcentbins->Fill(bin);
  //
  //Get Vertex
  //
  iEvent.getByToken(vertexToken,vertex_);
  const reco::VertexCollection * vertices3 = vertex_.product();
  vs_sell = vertices3->size();
  if(vs_sell>0) {
    vzr_sell = vertices3->begin()->z();
    vzErr_sell = vertices3->begin()->zError();
  } else {
    vzr_sell = -999.9;
  }
  vtx=vzr_sell;
  //
  //Get Event Planes
  //
  iEvent.getByToken(inputPlanesFlatToken,inputPlanesFlat_);
    if(!inputPlanesFlat_.isValid()){
    cout << "Error! Can't get hiEvtPlaneFlat product!" << endl;
    return ;
  }
  iEvent.getByToken(inputPlanesToken, inputPlanes_);
  if(!inputPlanes_.isValid()){
    cout << "Error! Can't get hiEvtPlane product!" << endl;
    return ;
  }
 
  Int_t indx = 0;
  for (EvtPlaneCollection::const_iterator rp = inputPlanesFlat_->begin();rp !=inputPlanesFlat_->end(); rp++) {
    epang[indx]=rp->angle();
    epsin[indx] = rp->sumSin();
    epcos[indx] = rp->sumCos();
    epqx[indx]  = rp->qx(); 
    epqy[indx]  = rp->qy();
    epq[indx]   = rp->q();
    epw[indx]   = rp->sumw();
    epmult[indx] = (double) rp->mult();
    if(centval<=80) Psi[indx]->Fill( rp->angle() );
    ++indx; 
  }

  indx = 0;
  for (EvtPlaneCollection::const_iterator rp = inputPlanes_->begin();rp !=inputPlanes_->end(); rp++) {
    epang_orig[indx] = rp->angle();
    if(centval<=80) PsiRaw[indx]->Fill( rp->angle() );
    ++indx; 
  }

  tree->Fill(); 
}



// ------------ method called once each job just before starting event loop  ------------
void 
CheckFlattening::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CheckFlattening::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CheckFlattening);

