
// -*- C++ -*-
//
// Package:    CheckFlattening
// Class:      CheckFlattening
// 


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"
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
  // edm::InputTag trackCollection_;
  //edm::InputTag vtxCollection_;
  //edm::InputTag inputPlanes_;
  //  edm::InputTag centrality_;

  //edm::InputTag centralityTag_;  
  //edm::EDGetTokenT<reco::Centrality> centralityToken;
  //edm::Handle<reco::Centrality> centrality_;

  edm::InputTag vertexTag_;
  //edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken;
  edm::Handle<std::vector<reco::Vertex>> vertex_;

  edm::InputTag inputPlanesTag_;
  //edm::EDGetTokenT<reco::EvtPlaneCollection> inputPlanesToken;
  edm::Handle<reco::EvtPlaneCollection> inputPlanes_;

  edm::InputTag trackTag_;
  //edm::EDGetTokenT<reco::TrackCollection> trackToken;
  edm::Handle<reco::TrackCollection> trackCollection_;

  CentralityProvider * centProvider;
  edm::Service<TFileService> fs;
  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  TH1D * hcent;
  TH1D * hcentbins;
  double centval;
  double ntrkval;
  double vtx;
  int noff;

  Double_t epang[NumEPNames];
  Double_t epsin[NumEPNames];
  Double_t epcos[NumEPNames];

  Double_t epang_orig[NumEPNames];
  Double_t epsin_orig[NumEPNames];
  Double_t epcos_orig[NumEPNames];

  Double_t epang_RecenterOnly[NumEPNames];
  Double_t epsin_RecenterOnly[NumEPNames];
  Double_t epcos_RecenterOnly[NumEPNames];


  Double_t epang_NoWgt[NumEPNames];
  Double_t epsin_NoWgt[NumEPNames];
  Double_t epcos_NoWgt[NumEPNames];

  Double_t sumw[NumEPNames];
  Double_t sumw2[NumEPNames];

  Double_t sumPtOrEt[NumEPNames];
  Double_t sumPtOrEt2[NumEPNames];

  Double_t qx[NumEPNames];
  Double_t qy[NumEPNames];
  Double_t q[NumEPNames];
  Double_t epmult[NumEPNames];

  Double_t rescor[NumEPNames];
  Double_t rescorErr[NumEPNames];
  unsigned int runno_;

  TH1D * hNtrkoff;
  int nEtaBins;
  TH1I * hrun;
  string rpnames[NumEPNames];
  TH1D * PsiRaw[50];
  TH1D * Psi[50];
  TTree * tree;

  int FlatOrder_;
  int NumFlatBins_;
  int CentBinCompression_;
  int HFEtScale_;
  int Noffmin_;
  int Noffmax_;

  HiEvtPlaneFlatten * flat[NumEPNames];
  bool loadDB_;
  bool FirstEvent_;

  bool Branch_Cent;
  bool Branch_NtrkOff;
  bool Branch_Vtx;
  bool Branch_epang;
  bool Branch_epang_orig;
  bool Branch_epang_RecenterOnly;
  bool Branch_NoWgt;
  bool Branch_epsin;
  bool Branch_epcos;
  bool Branch_epsin_orig;
  bool Branch_epcos_orig;
  bool Branch_epsin_RecenterOnly;
  bool Branch_epcos_RecenterOnly;
  bool Branch_epsin_NoWgt;
  bool Branch_epcos_NoWgt;
  bool Branch_sumw;
  bool Branch_sumw2;
  bool Branch_sumPtOrEt;
  bool Branch_sumPtOrEt2;
  bool Branch_qx;
  bool Branch_qy;
  bool Branch_q;
  bool Branch_mult;
  bool Branch_Run;
  bool Branch_Rescor;
  bool Branch_RescorErr;


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
// constructors and destructor
//
CheckFlattening::CheckFlattening(const edm::ParameterSet& iConfig):runno_(0)
  
{
  runno_ = 0;
  centProvider = 0;
  loadDB_ = kTRUE;
  FirstEvent_ = kTRUE;
  //centralityTag_ = iConfig.getParameter<edm::InputTag>("centralityTag_");
  //centralityToken = consumes<reco::Centrality>(centralityTag_);

  vertexTag_  = iConfig.getParameter<edm::InputTag>("vertexTag_");
  //vertexToken = consumes<std::vector<reco::Vertex>>(vertexTag_);

  trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag_");
  //trackToken = consumes<reco::TrackCollection>(trackTag_);

  inputPlanesTag_ = iConfig.getParameter<edm::InputTag>("inputPlanesTag_");
  //inputPlanesToken = consumes<reco::EvtPlaneCollection>(inputPlanesTag_);

  FlatOrder_ = iConfig.getUntrackedParameter<int>("FlatOrder_", 9);
  NumFlatBins_ = iConfig.getUntrackedParameter<int>("NumFlatBins_",20);
  CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",5);
  HFEtScale_ = iConfig.getUntrackedParameter<int>("HFEtScale_",3800);
  Noffmin_ = iConfig.getUntrackedParameter<int>("Noffmin_", 0);
  Noffmax_ = iConfig.getUntrackedParameter<int>("Noffmax_", 50000);	
  hNtrkoff = fs->make<TH1D>("Ntrkoff","Ntrkoff",1001,0,3000);

  Branch_Cent = iConfig.getUntrackedParameter<bool>("Branch_Cent",true);
  Branch_NtrkOff = iConfig.getUntrackedParameter<bool>("Branch_NtrkOff",true);
  Branch_Vtx = iConfig.getUntrackedParameter<bool>("Branch_Vtx",true);
  Branch_epang = iConfig.getUntrackedParameter<bool>("Branch_epang",true);
  Branch_epang_orig = iConfig.getUntrackedParameter<bool>("Branch_epang_orig",true);
  Branch_epang_RecenterOnly = iConfig.getUntrackedParameter<bool>("Branch_epang_RecenterOnly",true);
  Branch_NoWgt = iConfig.getUntrackedParameter<bool>("Branch_NoWgt",true);
  Branch_epsin = iConfig.getUntrackedParameter<bool>("Branch_epsin",true);
  Branch_epcos = iConfig.getUntrackedParameter<bool>("Branch_epcos",true);
  Branch_epsin_orig = iConfig.getUntrackedParameter<bool>("Branch_epsin_orig",true);
  Branch_epcos_orig = iConfig.getUntrackedParameter<bool>("Branch_epcos_orig",true);
  Branch_epsin_RecenterOnly = iConfig.getUntrackedParameter<bool>("Branch_epsin_RecenterOnly",true);
  Branch_epcos_RecenterOnly = iConfig.getUntrackedParameter<bool>("Branch_epcos_RecenterOnly",true);
  Branch_epsin_NoWgt = iConfig.getUntrackedParameter<bool>("Branch_epsin_NoWgt",true);
  Branch_epcos_NoWgt = iConfig.getUntrackedParameter<bool>("Branch_epcos_NoWgt",true);
  Branch_sumw = iConfig.getUntrackedParameter<bool>("Branch_sumw",true);
  Branch_sumw2 = iConfig.getUntrackedParameter<bool>("Branch_sumw2",true);
  Branch_sumPtOrEt = iConfig.getUntrackedParameter<bool>("Branch_sumPtOrEt",true);
  Branch_sumPtOrEt2 = iConfig.getUntrackedParameter<bool>("Branch_sumPtOrEt2",true);
  Branch_qx = iConfig.getUntrackedParameter<bool>("Branch_qx",true);
  Branch_qy = iConfig.getUntrackedParameter<bool>("Branch_qy",true);
  Branch_q = iConfig.getUntrackedParameter<bool>("Branch_q",true);
  Branch_mult = iConfig.getUntrackedParameter<bool>("Branch_mult",true);
  Branch_Run = iConfig.getUntrackedParameter<bool>("Branch_Run",true);
  Branch_Rescor = iConfig.getUntrackedParameter<bool>("Branch_Rescor",true);
  Branch_RescorErr = iConfig.getUntrackedParameter<bool>("Branch_RescorErr",true);

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

    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->Init(FlatOrder_,NumFlatBins_,HFEtScale_,EPNames[i],EPOrder[i]);

  }
  tree = fs->make<TTree>("tree","EP tree");

  if(Branch_Cent)              tree->Branch("Cent",&centval,"cent/D");
  if(Branch_NtrkOff)           tree->Branch("NtrkOff",&noff,"noff/D");
  if(Branch_Vtx)               tree->Branch("Vtx",&vtx,"vtx/D");
  if(Branch_epang)             tree->Branch("epang",&epang, epnames.Data());
  if(Branch_epang_orig)        tree->Branch("epang_orig",&epang_orig, epnames.Data());
  if(Branch_epang_RecenterOnly)  tree->Branch("epang_RecenterOnly", &epang_RecenterOnly, epnames.Data());
  if(Branch_NoWgt)             tree->Branch("epang_NoWgt", &epang_NoWgt, epnames.Data());

  if(Branch_epsin)             tree->Branch("epsin",     &epsin,      epnames.Data());
  if(Branch_epcos)             tree->Branch("epcos",     &epcos,      epnames.Data());
  if(Branch_epsin_orig)        tree->Branch("epsin_orig",     &epsin_orig,      epnames.Data());
  if(Branch_epcos_orig)        tree->Branch("epcos_orig",     &epcos_orig,      epnames.Data());
  if(Branch_epsin_RecenterOnly)tree->Branch("epsin_RecenterOnly",     &epsin_RecenterOnly,      epnames.Data());
  if(Branch_epcos_RecenterOnly)tree->Branch("epcos_RecenterOnly",     &epcos_RecenterOnly,      epnames.Data());
  if(Branch_epsin_NoWgt)       tree->Branch("epsin_NoWgt",     &epsin_NoWgt,      epnames.Data());
  if(Branch_epcos_NoWgt)       tree->Branch("epcos_NoWgt",     &epcos_NoWgt,      epnames.Data());
  if(Branch_sumw)              tree->Branch("sumw",  &sumw,        epnames.Data());
  if(Branch_sumw2)             tree->Branch("sumw2",  &sumw2,        epnames.Data());
  if(Branch_sumPtOrEt)         tree->Branch("sumPtOrEt",  &sumPtOrEt,        epnames.Data());
  if(Branch_sumPtOrEt2)        tree->Branch("sumPtOrEt2",  &sumPtOrEt2,        epnames.Data());

  if(Branch_qx)                tree->Branch("qx",      &qx,       epnames.Data());
  if(Branch_qy)                tree->Branch("qy",      &qy,       epnames.Data());
  if(Branch_q)                 tree->Branch("q",       &q,       epnames.Data());
  if(Branch_mult)              tree->Branch("mult",    &epmult,  epnames.Data());
  if(Branch_Run)               tree->Branch("Run",     &runno_,   "run/i");
  if(Branch_Rescor)            tree->Branch("Rescor",  &rescor,   epnames.Data());
  if(Branch_RescorErr)         tree->Branch("RescorErr",  &rescorErr,   epnames.Data());

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
  Bool_t newrun = kFALSE;
  if(runno_ != iEvent.id().run()) newrun = kTRUE;
  runno_ = iEvent.id().run();
  hrun->Fill(runno_);

  if(FirstEvent_ || newrun) {
    FirstEvent_ = kFALSE;
    newrun = kFALSE;
    cout<<"Load flattening parameters"<<endl;
    //
    //Get flattening parameter file.  
    //
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    LoadEPDB * db = new LoadEPDB(flatparmsDB_,flat);
    if(!db->IsSuccess()) {
      loadDB_ = kFALSE;
    }        
  } //First event

  //
  //Get Centrality
  //
  int bin = 0;
  int Noff = getNoff( iEvent, iSetup);
  ntrkval = Noff;
  if ( (Noff < Noffmin_) or (Noff >= Noffmax_) ) {
    return;
  }
  hNtrkoff->Fill(Noff);


  if(!centProvider) centProvider = new CentralityProvider(iSetup);
  centProvider->newEvent(iEvent,iSetup);
  centProvider->raw();
  bin = centProvider->getBin();
  centval = (100./centProvider->getNbins())*(bin+0.5);
  hcent->Fill(centval);
  hcentbins->Fill(bin);
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
  
  vtx = vzr_sell;
  //
  //Get Event Planes
  //
  iEvent.getByLabel(inputPlanesTag_,inputPlanes_);
  
  if(!inputPlanes_.isValid()){
    cout << "Error! Can't get hiEvtPlane product!" << endl;
    return ;
  }
  
  Int_t indx = 0;
  for (EvtPlaneCollection::const_iterator rp = inputPlanes_->begin();rp !=inputPlanes_->end(); rp++) {
    if(indx != rp->indx() ) {
      cout<<"EP inconsistency found. Abort."<<endl;
      return;
    }
    epang[indx]=rp->angle();
    epsin[indx] = rp->sumSin();
    epcos[indx] = rp->sumCos();

    epang_orig[indx]=rp->angle(0);
    epsin_orig[indx] = rp->sumSin(0);
    epcos_orig[indx] = rp->sumCos(0);

    epang_RecenterOnly[indx]=rp->angle(1);
    epsin_RecenterOnly[indx] = rp->sumSin(1);
    epcos_RecenterOnly[indx] = rp->sumCos(1);

    epang_NoWgt[indx]=rp->angle(3);
    epsin_NoWgt[indx] = rp->sumSin(3);
    epcos_NoWgt[indx] = rp->sumCos(3);

    qx[indx]  = rp->qx(); 
    qy[indx]  = rp->qy();
    q[indx]   = rp->q();
    sumw[indx]   = rp->sumw();
    sumw2[indx] = rp->sumw2();
    sumPtOrEt[indx] = rp->sumPtOrEt();
    sumPtOrEt2[indx] = rp->sumPtOrEt2();

    epmult[indx] = (double) rp->mult();

    if(Branch_Rescor || Branch_RescorErr) {
      rescor[indx] = flat[indx]->GetCentRes1((int) centval);
      rescorErr[indx] = flat[indx]->GetCentResErr1((int) centval);
    }
    if(centval<=80) {
      Psi[indx]->Fill( rp->angle() );
      PsiRaw[indx]->Fill( rp->angle(0) );
    }
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

