// -*- C++ -*-
//
// Package:    HiEvtPlaneFlatCalib
// Class:      HiEvtPlaneFlatCalib
// 
/**\class HiEvtPlaneFlatCalib HiEvtPlaneFlatCalib.cc HiEvtPlaneFlatten/HiEvtPlaneFlatCalib/src/HiEvtPlaneFlatCalib.cc


 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Stephen Sanders
//         Created:  Sat Jun 26 16:04:04 EDT 2010
// $Id: HiEvtPlaneFlatCalib.cc,v 1.7 2012/02/15 11:04:09 eulisse Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"

#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"

#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TList.h"
#include "TString.h"
#include <time.h>
#include <cstdlib>
#include <vector>

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"

using namespace std;
using namespace hi;

//
// class declaration
//

class HiEvtPlaneFlatCalib : public edm::EDAnalyzer {
   public:
      explicit HiEvtPlaneFlatCalib(const edm::ParameterSet&);
      ~HiEvtPlaneFlatCalib();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  edm::InputTag vtxCollection_;
  edm::InputTag inputPlanes_;

  //  const CentralityBins * cbins_;
  CentralityProvider * centrality_;
  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  TTree * tree;
  double centval;
  double vtx;
  Double_t epang[NumEPNames];
  Double_t epangorig[NumEPNames];
  unsigned int runno_;
 
  int FlatOrder_;
  int NumFlatCentBins_;
  int CentBinCompression_;

  static const  int NumCentBins=9;
  double wcent[10];

  TH1D * hcent;
  TH1D * hvtx;
  TH1D * flatXhist[NumEPNames];
  TH1D * flatYhist[NumEPNames];
  TH1D * flatCnthist[NumEPNames];

  TH1D * flatXDBhist[NumEPNames];
  TH1D * flatYDBhist[NumEPNames];
  TH1D * hPsi[NumEPNames];
  TH1D * hPsiFlat[NumEPNames];
  TH1D * hPsiFlatCent[NumEPNames][NumCentBins];
  TH1D * hPsiFlatSub1[NumEPNames];
  TH1D * hPsiFlatSub2[NumEPNames];
  HiEvtPlaneFlatten * flat[NumEPNames];
  RPFlatParams * rpFlat;
  int nRP;
  bool genFlatPsi_;

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
HiEvtPlaneFlatCalib::HiEvtPlaneFlatCalib(const edm::ParameterSet& iConfig)
{
  genFlatPsi_ = iConfig.getUntrackedParameter<bool>("genFlatPsi_",true);
  vtxCollection_  = iConfig.getParameter<edm::InputTag>("vtxCollection_");
  inputPlanes_ = iConfig.getParameter<edm::InputTag>("inputPlanes_");
  FlatOrder_ = iConfig.getUntrackedParameter<int>("FlatOrder_", 9);
  NumFlatCentBins_ = iConfig.getUntrackedParameter<int>("NumFlatCentBins_",50);
  CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",4);
  //  NumCentBins=9;
  wcent[0] = 0;
  wcent[1] = 5;
  wcent[2] = 10;
  wcent[3] = 20;
  wcent[4] = 30;
  wcent[5] = 40;
  wcent[6] = 50;
  wcent[7] = 60;
  wcent[8] = 70;
  wcent[9] = 100;

  //now do what ever other initialization is needed
  //  cbins_ = 0;
  centrality_ = 0;
  hcent = fs->make<TH1D>("cent","cent",201,0,200);
  hvtx = fs->make<TH1D>("vtx","vtx",1000,-50,50);
  //setting before 8:35 CDT on 11Jan2011

  TString epnames = EPNames[0].data();
  epnames = epnames+"/D";
  for(int i = 0; i<NumEPNames; i++) {
    if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/D";
    TFileDirectory subdir = fs->mkdir(Form("%s",EPNames[i].data()));
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->Init(FlatOrder_,NumFlatCentBins_,CentBinCompression_,EPNames[i],EPOrder[i]);
    int nbins = flat[i]->GetHBins();
    flatXhist[i] = subdir.make<TH1D>(Form("x_%s",EPNames[i].data()),Form("x_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    flatYhist[i] = subdir.make<TH1D>(Form("y_%s",EPNames[i].data()),Form("y_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    flatXDBhist[i] = subdir.make<TH1D>(Form("xDB_%s",EPNames[i].data()),Form("xDB_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    flatYDBhist[i] = subdir.make<TH1D>(Form("yDB_%s",EPNames[i].data()),Form("yDB_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    flatCnthist[i] = subdir.make<TH1D>(Form("cnt_%s",EPNames[i].data()),Form("cnt_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    Double_t psirange = 4;
    if(EPOrder[i]==2 ) psirange = 2;
    if(EPOrder[i]==3 ) psirange = 1.5;
    if(EPOrder[i]==4 ) psirange = 1;
    if(EPOrder[i]==5) psirange = 0.8;
    if(EPOrder[i]==6) psirange = 0.6;
    hPsi[i] = subdir.make<TH1D>("psi","psi",800,-psirange,psirange);
    hPsi[i]->SetXTitle("#Psi");
    hPsi[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    hPsiFlat[i] = subdir.make<TH1D>("psiFlat","psiFlat",800,-psirange,psirange);
    hPsiFlat[i]->SetXTitle("#Psi");
    hPsiFlat[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    for(int j = 0; j<NumCentBins; j++) {
      TString hname = Form("psiFlat_%d_%d",(int) wcent[j],(int) wcent[j+1]);
      hPsiFlatCent[i][j] = subdir.make<TH1D>(hname.Data(),hname.Data(),800,-psirange,psirange);
      hPsiFlatCent[i][j]->SetXTitle("#Psi");
      hPsiFlatCent[i][j]->SetYTitle(Form("Counts (%d<cent#leq%d%c)",(int) wcent[j],(int) wcent[j+1],'%'));
    }  

  }
  tree = fs->make<TTree>("tree","EP tree");
  tree->Branch("Cent",&centval,"cent/D");
  tree->Branch("Vtx",&vtx,"vtx/D");
  tree->Branch("EP",&epang, epnames.Data());
  tree->Branch("EP_orig",&epangorig, epnames.Data());
  tree->Branch("Run",&runno_,"run/i");
  
}


HiEvtPlaneFlatCalib::~HiEvtPlaneFlatCalib()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HiEvtPlaneFlatCalib::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  runno_ = iEvent.id().run(); 
  //
  //Get Centrality
  //
  if(!centrality_) centrality_ = new CentralityProvider(iSetup);

   centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
   int bin = centrality_->getBin();
   centval = (100./centrality_->getNbins())*(bin+0.5);
   hcent->Fill(bin);
  //
  //Get Vertex
  //
  edm::Handle<reco::VertexCollection> vertexCollection3;
  iEvent.getByLabel(vtxCollection_,vertexCollection3);
  const reco::VertexCollection * vertices3 = vertexCollection3.product();
  vs_sell = vertices3->size();
  if(vs_sell>0) {
    vzr_sell = vertices3->begin()->z();
    vzErr_sell = vertices3->begin()->zError();
  } else {
    vzr_sell = -999.9;
  }
  vtx=vzr_sell;
  hvtx->Fill(vzr_sell);
  //
  //Get Flattening Parameters
  //
  if(genFlatPsi_) {
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    int flatTableSize = flatparmsDB_->m_table.size();
    for(int i = 0; i<flatTableSize; i++) {
      const RPFlatParams::EP* thisBin = &(flatparmsDB_->m_table[i]);
      for(int j = 0; j<NumEPNames; j++) {
	int indx = thisBin->RPNameIndx[j];
	if(indx>=0) {
	    flat[indx]->SetXDB(i, thisBin->x[j]);
	    flat[indx]->SetYDB(i, thisBin->y[j]);
	    flatXDBhist[j]->SetBinContent(i+1,thisBin->x[j]);
	    flatYDBhist[j]->SetBinContent(i+1,thisBin->y[j]);
	}
      }
    }
  }
  //
  //Get Event Planes
  //

  Handle<reco::EvtPlaneCollection> evtPlanes;
  iEvent.getByLabel(inputPlanes_,evtPlanes);

  if(!evtPlanes.isValid()){
    //cout << "Error! Can't get hiEvtPlane product!" << endl;
    return ;
  }

  for (EvtPlaneCollection::const_iterator rp = evtPlanes->begin();rp !=evtPlanes->end(); rp++) {
    string baseName = rp->label();
    if(rp->angle() > -5) {
      for(int i = 0; i< NumEPNames; i++) {
	if(EPNames[i].compare(baseName)==0) {
	  double psiFlat = -5;
	  epang[i] = psiFlat;
	  if(genFlatPsi_){
	    psiFlat = flat[i]->GetFlatPsi(rp->angle(),vzr_sell,bin);
	    epang[i]=psiFlat;
	  }
	  epangorig[i]=rp->angle();
	  flat[i]->Fill(rp->angle(),vzr_sell,bin);
	  
	  if(centval<=80) hPsi[i]->Fill(rp->angle());
	  
	  if(genFlatPsi_) {
	    if(centval<=80) hPsiFlat[i]->Fill(psiFlat);
	    for(int j = 0; j<NumCentBins; j++) {
	      if(centval>wcent[j]&&centval<=wcent[j+1]) hPsiFlatCent[i][j]->Fill(psiFlat);
	    }
	  } 
	}
      }
    }    
  }
  tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
HiEvtPlaneFlatCalib::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiEvtPlaneFlatCalib::endJob() {
  for(int i = 0; i<NumEPNames; i++) {
    for(int j = 0; j<flat[i]->GetHBins();j++) {
      flatXhist[i]->SetBinContent(j+1,flat[i]->GetX(j));
      flatYhist[i]->SetBinContent(j+1,flat[i]->GetY(j));
      flatCnthist[i]->SetBinContent(j+1,flat[i]->GetCnt(j));
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiEvtPlaneFlatCalib);
