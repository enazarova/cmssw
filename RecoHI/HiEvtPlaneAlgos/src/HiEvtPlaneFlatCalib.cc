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
// $Id: HiEvtPlaneFlatCalib.cc,v 1.5 2011/11/06 23:17:46 ssanders Exp $
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TList.h"
#include "TString.h"
//#include "TRandom3.h"
#include <time.h>
#include <cstdlib>
#include <vector>

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"

using namespace std;
using namespace hi;

//
// class declaration
//

static const  int NumCentBins=14;
static double wcent[]={0,5,10,15,20,25,30,35,40,50,60,70,80,90,100};

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
  edm::InputTag caloCollection_;
  edm::InputTag trackCollection_;
  edm::InputTag inputPlanes_;

  TTree * tree;
  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  double centval;
  double vtx;
  unsigned int runno_;

  bool FirstEvent;

  TH1D * hcent;
  TH1D * hcentbins;
  TH1D * hvtx;
  TH1I * hruns;

  TH1D * flatXhist[NumEPNames];
  TH1D * flatYhist[NumEPNames];
  TH2D * hPsiPsiFlat[NumEPNames];	
  TH2D * hPsiPsiOffset[NumEPNames];	
  TH1D * flatCnthist[NumEPNames];

  TH1D * xoffhist[NumEPNames];
  TH1D * yoffhist[NumEPNames];
  TH1D * xyoffcnthist[NumEPNames];
  TH1D * xyoffmulthist[NumEPNames];

  TH1D * pthist[NumEPNames];
  TH1D * pt2hist[NumEPNames];
  TH1D * ptcnthist[NumEPNames];
  TH1D * avwpt[NumEPNames];
  TH1D * avwptev[NumEPNames];
  TH1D * avwptnorm[NumEPNames];
  TH1D * avwptcnt[NumEPNames];

  TH1D * flatXDBhist[NumEPNames];
  TH1D * flatYDBhist[NumEPNames];
  TH1D * xoffDBhist[NumEPNames];
  TH1D * yoffDBhist[NumEPNames];
  TH1D * cntoffDBhist[NumEPNames];
  TH1D * ptDBhist[NumEPNames];
  TH1D * cntptDBhist[NumEPNames];
  TH1D * pt2DBhist[NumEPNames];
  TH1D * ptMerge[NumEPNames];
  TH1D * hPsi[NumEPNames];
  TH1D * hPsiFlat[NumEPNames];
  TH1D * hPsiOffset[NumEPNames];
  TH1D * hPsiFlatCent[NumEPNames][NumCentBins];
  TH1D * hPsiFlatSub1[NumEPNames];
  TH1D * hPsiFlatSub2[NumEPNames];
  bool evtcnt[NumEPNames];
  Double_t epang[NumEPNames];
  CentralityProvider * centProvider;
  HiEvtPlaneFlatten * flat[NumEPNames];
  int nRP;
  int Hbins;
  int Obins;
  bool genFlatPsi_;
  bool useOffsetPsi_;
  bool useECAL_;
  bool useHCAL_;
  bool useTrack_;
  bool useMomentumCorrV1_;
  bool useTrackPtWeight_;
  double minet_;
  double maxet_;
  double effm_;
  double minpt_;
  double maxpt_;
  double minvtx_;
  double maxvtx_;
  double dzerr_;
  double chi2_;
  bool storeNames_;
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
HiEvtPlaneFlatCalib::HiEvtPlaneFlatCalib(const edm::ParameterSet& iConfig):runno_(0)
{
  vtxCollection_  = iConfig.getParameter<edm::InputTag>("vtxCollection_");
  caloCollection_  = iConfig.getParameter<edm::InputTag>("caloCollection_");
  trackCollection_  = iConfig.getParameter<edm::InputTag>("trackCollection_");
  inputPlanes_ = iConfig.getParameter<edm::InputTag>("inputPlanes_");
  genFlatPsi_ = iConfig.getUntrackedParameter<bool>("genFlatPsi_",true);
  useOffsetPsi_ = iConfig.getUntrackedParameter<bool>("useOffsetPsi_",true);
  useECAL_ = iConfig.getUntrackedParameter<bool>("useECAL_",true);
  useHCAL_ = iConfig.getUntrackedParameter<bool>("useHCAL_",true);
  useTrack_ = iConfig.getUntrackedParameter<bool>("useTrack",true);
  minet_ = iConfig.getUntrackedParameter<double>("minet_",0.3);
  maxet_ = iConfig.getUntrackedParameter<double>("maxet_",500.);
  effm_ = iConfig.getUntrackedParameter<double>("effm_",0.0);
  minpt_ = iConfig.getUntrackedParameter<double>("minpt_",0.3);
  maxpt_ = iConfig.getUntrackedParameter<double>("maxpt_",2.5);
  minvtx_ = iConfig.getUntrackedParameter<double>("minvtx_",-50.);
  maxvtx_ = iConfig.getUntrackedParameter<double>("maxvtx_",50.);
  dzerr_ = iConfig.getUntrackedParameter<double>("dzerr_",10.);
  chi2_  = iConfig.getUntrackedParameter<double>("chi2_",40.);
  storeNames_ = 1;
  FirstEvent = kTRUE;
  //now do what ever other initialization is needed
  hcent = fs->make<TH1D>("cent","cent",NumCentBins,wcent);
  hcentbins = fs->make<TH1D>("centbins","centbins",201,0,200);
  hvtx = fs->make<TH1D>("vtx","vtx",1000,-50,50);  
  hruns = fs->make<TH1I>("runs","runs",3001,181000,184000);
  cout<<"====================="<<endl;
  cout<<"HiEvtPlaneFlatCalib: "<<endl;
  cout<<"  minet_:            "<<minet_<<endl;
  cout<<"  maxet_:            "<<maxet_<<endl;
  cout<<"  minpt_:            "<<minpt_<<endl;
  cout<<"  maxpt_:            "<<maxpt_<<endl;
  cout<<"  minvtx_:           "<<minvtx_<<endl;
  cout<<"  maxvtx_:           "<<maxvtx_<<endl;
  cout<<"  dzerr_             "<<dzerr_<<endl;
  cout<<"  chi2_              "<<chi2_<<endl;
  cout<<"====================="<<endl;
  TString epnames = EPNames[0].data();
  epnames = epnames+"/D";
  for(int i = 0; i<NumEPNames; i++) {
    if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/D";
    TFileDirectory subdir = fs->mkdir(Form("%s",EPNames[i].data()));
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->Init(FlatOrder,NumFlatCentBins,CentBinCompression,EPNames[i],EPOrder[i]);
    Hbins = flat[i]->GetHBins();
    Obins = flat[i]->GetOBins();
    int nbins = flat[i]->GetHBins() + 2*flat[i]->GetOBins();
    flatXhist[i] = subdir.make<TH1D>(Form("x_%s",EPNames[i].data()),Form("x_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    flatYhist[i] = subdir.make<TH1D>(Form("y_%s",EPNames[i].data()),Form("y_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    flatXDBhist[i] = subdir.make<TH1D>(Form("xDB_%s",EPNames[i].data()),Form("x_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    flatYDBhist[i] = subdir.make<TH1D>(Form("yDB_%s",EPNames[i].data()),Form("y_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);
    flatCnthist[i] = subdir.make<TH1D>(Form("cnt_%s",EPNames[i].data()),Form("cnt_%s",EPNames[i].data()),nbins,-0.5,nbins-0.5);

    avwpt[i] = subdir.make<TH1D>(Form("avwpt_%s",EPNames[i].data()),Form("avwpt_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    avwpt[i]->SetXTitle("bin");
    avwpt[i]->SetYTitle("<w*p_{T}>");
    avwptcnt[i] = subdir.make<TH1D>(Form("avwptcnt_%s",EPNames[i].data()),Form("avwptcnt_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    avwptnorm[i] = subdir.make<TH1D>(Form("avwptnorm_%s",EPNames[i].data()),Form("avwptnorm_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    avwptev[i] = subdir.make<TH1D>(Form("avwptev_%s",EPNames[i].data()),Form("avwptev_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);

    xoffhist[i] = subdir.make<TH1D>(Form("xoff_%s",EPNames[i].data()),Form("xoff_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    yoffhist[i] = subdir.make<TH1D>(Form("yoff_%s",EPNames[i].data()),Form("yoff_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    xyoffcnthist[i] = subdir.make<TH1D>(Form("xyoffcnt_%s",EPNames[i].data()),Form("xyoffcnt_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    xyoffmulthist[i] = subdir.make<TH1D>(Form("xyoffmult_%s",EPNames[i].data()),Form("xyoffmult_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    xoffDBhist[i] = subdir.make<TH1D>(Form("xoffDB_%s",EPNames[i].data()),Form("xoffDB_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    yoffDBhist[i] = subdir.make<TH1D>(Form("yoffDB_%s",EPNames[i].data()),Form("yoffDB_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    cntoffDBhist[i] = subdir.make<TH1D>(Form("cntoffDB_%s",EPNames[i].data()),Form("cntoffDB_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);

    pthist[i] = subdir.make<TH1D>(Form("pt_%s",EPNames[i].data()),Form("pt_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    ptMerge[i] = subdir.make<TH1D>(Form("ptMerge_%s",EPNames[i].data()),Form("ptMerge_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    for(int j = 0; j< ptMerge[i]->GetNbinsX(); j++) ptMerge[i]->SetBinContent(j+1,1);
    pt2hist[i] = subdir.make<TH1D>(Form("pt2_%s",EPNames[i].data()),Form("pt2_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    ptcnthist[i] = subdir.make<TH1D>(Form("ptcnt_%s",EPNames[i].data()),Form("ptcnt_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    ptDBhist[i] = subdir.make<TH1D>(Form("ptDB_%s",EPNames[i].data()),Form("ptDB_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    cntptDBhist[i] = subdir.make<TH1D>(Form("cntptDB_%s",EPNames[i].data()),Form("cntptDB_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);
    pt2DBhist[i] = subdir.make<TH1D>(Form("pt2DB_%s",EPNames[i].data()),Form("pt2DB_%s",EPNames[i].data()),Obins,-0.5,Obins-0.5);


    Double_t psirange = 4;
    if(EPOrder[i]==1 ) psirange = 3.5;
    if(EPOrder[i]==2 ) psirange = 2;
    if(EPOrder[i]==3 ) psirange = 1.5;
    if(EPOrder[i]==4 ) psirange = 1;
    if(EPOrder[i]==5) psirange = 0.8;
    if(EPOrder[i]==6) psirange = 0.6;
    //std::cout<<EPNames[i].data()<<"\t"<<EPOrder[i]<<"\t"<<psirange<<endl;
    hPsi[i] = subdir.make<TH1D>("psi","psi",800,-psirange,psirange);
    hPsi[i]->SetXTitle("#Psi");
    hPsi[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    hPsiFlat[i] = subdir.make<TH1D>("psiFlat","psiFlat",800,-psirange,psirange);
    hPsiFlat[i]->SetXTitle("#Psi");
    hPsiFlat[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    hPsiOffset[i] = subdir.make<TH1D>("psiOffset","psiOffset",800,-psirange,psirange);
    hPsiOffset[i]->SetXTitle("#Psi");
    hPsiOffset[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    hPsiPsiFlat[i] = subdir.make<TH2D>("PsiPsiFlat","PsiFlat vs. Psi",80,-psirange,psirange,80,-psirange,psirange);
    hPsiPsiFlat[i]->SetOption("colz");
    hPsiPsiFlat[i]->SetXTitle("#Psi_{flat}");
    hPsiPsiFlat[i]->SetYTitle("#Psi");

    hPsiPsiOffset[i] = subdir.make<TH2D>("PsiPsiOffset","PsiOfset vs. PsiOrig",80,-psirange,psirange,80,-psirange,psirange);
    hPsiPsiOffset[i]->SetOption("colz");
    hPsiPsiOffset[i]->SetXTitle("#Psi_{flat}");
    hPsiPsiOffset[i]->SetYTitle("#Psi");
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
  //cout<<"====NEW EVENT===="<<endl;
  if(FirstEvent) {
    FirstEvent = kFALSE;
    //
    //Check if new run
    //

    //
    //Get Flattening Parameters
    //
    if(genFlatPsi_ ) {
      edm::ESHandle<RPFlatParams> flatparmsDB_;
      iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
      int flatTableSize = flatparmsDB_->m_table.size();
      //cout<<"flatTableSize: "<<flatTableSize<<endl;
      for(int i = 0; i<flatTableSize; i++) {
	const RPFlatParams::EP* thisBin = &(flatparmsDB_->m_table[i]);
	for(int j = 0; j<NumEPNames; j++) {
	  int indx = thisBin->RPNameIndx[j];
	  if(indx>=0) {
	    if(i<Hbins) {
	      flat[indx]->SetXDB(i, thisBin->x[j]);
	      flat[indx]->SetYDB(i, thisBin->y[j]);
	      flatXDBhist[indx]->SetBinContent(i+1,thisBin->x[j]);
	      flatYDBhist[indx]->SetBinContent(i+1,thisBin->y[j]);
	    } else if(i>=Hbins && i<Hbins+Obins) {
	      flat[indx]->SetXoffDB(i - Hbins, thisBin->x[j]);
	      flat[indx]->SetYoffDB(i - Hbins, thisBin->y[j]);
	      xoffDBhist[indx]->SetBinContent(i+1-Hbins,thisBin->x[j]);
	      yoffDBhist[indx]->SetBinContent(i+1-Hbins,thisBin->y[j]);
	      cntoffDBhist[indx]->SetBinContent(i+1-Hbins,1);

	    } else if (i>=Hbins+Obins) {
	      flat[indx]->SetPtDB(i - Hbins- Obins, thisBin->x[j]);
	      flat[indx]->SetPt2DB(i - Hbins- Obins, thisBin->y[j]);
	      ptDBhist[indx]->SetBinContent(i+1-Hbins-Obins,thisBin->x[j]);
	      cntptDBhist[indx]->SetBinContent(i+1-Hbins-Obins,1);
	      pt2DBhist[indx]->SetBinContent(i+1-Hbins-Obins,thisBin->y[j]);
	    }
	  }
	}
      }
    }
  } //First event

  //
  //Get Centrality
  //
  centProvider = new CentralityProvider(iSetup);
  centProvider->newEvent(iEvent,iSetup);
  centProvider->raw();
  int bin = centProvider->getBin();
  centval = (100./centProvider->getNbins())*(bin+0.5);
  //  if(hcentbins->Integral(1,1000)>0 || centval>50) return;  //For one-event test...
  hcent->Fill(centval);
  hcentbins->Fill(bin);

  //
  //Get Vertex
  //
  edm::Handle<reco::VertexCollection> vertexCollection3;
  iEvent.getByLabel("hiSelectedVertex",vertexCollection3);
  const reco::VertexCollection * vertices3 = vertexCollection3.product();
  vs_sell = vertices3->size();
  if(vs_sell>0) {
    vzr_sell = vertices3->begin()->z();
    vzErr_sell = vertices3->begin()->zError();
  } else
    vzr_sell = -999.9;

  vtx = vzr_sell;

  hruns->Fill(runno_);
  //Set up momentum correction weights
  //
    
    double tower_eta;
    double tower_energyet, tower_energyet_e, tower_energyet_h;
    Handle<CaloTowerCollection> calotower;
    iEvent.getByLabel(caloCollection_,calotower);
    
    if(calotower.isValid()){

	for(int i = 0; i<NumEPNames; i++) {
	  evtcnt[i] = kTRUE;
	  avwptev[i]->Reset();
	  avwptnorm[i]->Reset();
	}

	for (CaloTowerCollection::const_iterator j = calotower->begin();j !=calotower->end(); j++) {   
	  tower_eta        = j->eta();
	  tower_energyet_e   = j->emEt();
	  tower_energyet_h   = j->hadEt();
	  tower_energyet     = tower_energyet_e + tower_energyet_h;
	  if(tower_energyet<minet_) continue;
	  if(tower_energyet>maxet_) continue;
	  double ept = sqrt(pow(tower_energyet,2.)-pow(effm_,2.));
	  for(int i = 0; i<NumEPNames; i++) {
	    if(EPDet[i]==HF && EPOrder[i]==1) {
	      if((tower_eta>=EPEtaMin1[i]&&tower_eta<EPEtaMax1[i] )||
		 (EPEtaMin2[i]!=EPEtaMax2[i]&&(tower_eta>=EPEtaMin2[i]&&tower_eta<EPEtaMax2[i]))) {
		flat[i]->FillPt(ept,vtx,bin);
		double w = flat[i]->GetW(ept,vtx,bin);
		if(tower_eta<0) w=-w;
		int evbin = flat[i]->GetOffsetIndx(bin,vtx);
		avwptev[i]->Fill(evbin, w*ept);
	       
		avwptnorm[i]->Fill(evbin);
		if(evtcnt[i]) {
		  evtcnt[i]=kFALSE;
		  avwptcnt[i]->Fill(evbin);
		}
	      }
	    }
	  }
	} 
	for(int i = 0; i<NumEPNames; i++) {
	  if(EPDet[i]==HF && avwptnorm[i]->Integral(1,avwptnorm[i]->GetNbinsX())>0) {
	    avwptev[i]->Divide(avwptnorm[i]);
	    avwpt[i]->Add(avwptev[i]);
	  }
	}  


    }
    double track_eta;
    double track_pt;
    double track_charge;
    Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel(trackCollection_, tracks);
    if(tracks.isValid()){
	for(int i = 0; i<NumEPNames; i++) {
	  avwptnorm[i]->Reset();
	  avwptev[i]->Reset();
	  evtcnt[i] = kTRUE;
	}
	for(reco::TrackCollection::const_iterator j = tracks->begin(); j != tracks->end(); j++){
	  edm::Handle<reco::VertexCollection> vertex;
	  iEvent.getByLabel(vtxCollection_, vertex);
	  math::XYZPoint vtxPoint(0.0,0.0,0.0);
	  double vzErr =0.0, vxErr=0.0, vyErr=0.0;
	  if(vertex->size()>0) {
	    vtxPoint=vertex->begin()->position();
	    vzErr=vertex->begin()->zError();
	    vxErr=vertex->begin()->xError();
	    vyErr=vertex->begin()->yError();
	  }
	  bool accepted = true;
	  bool isPixel = false;
	  if ( j->numberOfValidHits() < 7 ) isPixel = true;
	  double d0=0.0, dz=0.0, d0sigma=0.0, dzsigma=0.0;
	  d0 = -1.*j->dxy(vtxPoint);
	  dz = j->dz(vtxPoint);
	  d0sigma = sqrt(j->d0Error()*j->d0Error()+vxErr*vyErr);
	  dzsigma = sqrt(j->dzError()*j->dzError()+vzErr*vzErr);
	  if( isPixel ) {
	    if ( fabs(dz/dzsigma) > dzerr_ ) accepted = false;
	    if ( j->normalizedChi2() > chi2_ ) accepted = false;
	  }
	  if ( ! isPixel){
	    if ( fabs(dz/dzsigma) > 3 ) accepted = false;
	    if ( fabs(d0/d0sigma) > 3 ) accepted = false;
	    if ( j->ptError()/j->pt() > 0.05 ) accepted = false;
	    if ( j->numberOfValidHits() < 12 ) accepted = false;
	  }
	  if( accepted ) {
	    track_eta = j->eta();
	    track_pt = j->pt();
	    track_charge = j->charge();
	    if(track_pt<minpt_) continue;
	    if(track_pt>maxpt_) continue;
	    if(fabs(track_eta)<=2.4 && track_charge!=0) {
	      for(int i = 0; i<NumEPNames; i++) {
		if(EPDet[i]==Tracker && EPOrder[i]==1) {
		  if((track_eta>=EPEtaMin1[i]&&track_eta<EPEtaMax1[i]) ||
		     (EPEtaMin2[i]!=EPEtaMax2[i]&&(track_eta>=EPEtaMin2[i]&&track_eta<EPEtaMax2[i]))) {
		    flat[i]->FillPt(track_pt,vtx,bin);
		    double w = flat[i]->GetW(track_pt,vtx,bin);
		    if(track_eta < 0) w=-w;

		    int evbin = flat[i]->GetOffsetIndx(bin,vtx);
		    avwptev[i]->Fill(evbin, w*track_pt);
		    avwptnorm[i]->Fill(flat[i]->GetOffsetIndx(bin,vtx));
		    if(evtcnt[i]) {
		      avwptcnt[i]->Fill(flat[i]->GetOffsetIndx(bin,vtx));
		      evtcnt[i] = kFALSE;
		    }
		  }
		}
	      }
	    }
	  } 
	}  
	for(int i = 0; i<NumEPNames; i++) {
	  if(EPDet[i]==Tracker && avwptnorm[i]->Integral(1,avwptnorm[i]->GetNbinsX())>0) {

	    avwptev[i]->Divide(avwptnorm[i]);
	    avwpt[i]->Add(avwptev[i]);
	  }
	}  
    }

  //
  //Get Event Planes
  //
  Handle<reco::EvtPlaneCollection> evtPlanes;
  iEvent.getByLabel(inputPlanes_,evtPlanes);

  if(!evtPlanes.isValid()){
    cout << "Error! Can't get hiEvtPlane product!" << endl;
    return ;
  }

  double plotZrange = 25;
  for (EvtPlaneCollection::const_iterator rp = evtPlanes->begin();rp !=evtPlanes->end(); rp++) {
    string tmpName = rp->label();
    if(rp->angle() > -5) {
      string baseName = rp->label();
      for(int i = 0; i< NumEPNames; i++) {
	if(EPNames[i].compare(baseName)==0) {
	  double angorig = rp->angle();
	  double c = rp->sumCos();
	  double s = rp->sumSin();
	  uint m = rp->mult();
	  double psiOffset = angorig;
	  if(useOffsetPsi_) psiOffset = flat[i]->GetOffsetPsi(s,c,vtx,bin);
	  double psiFlat = flat[i]->GetFlatPsi(psiOffset,vtx,bin);
	  epang[i]=psiFlat;
	  flat[i]->Fill(psiOffset,vtx,bin);
	  flat[i]->FillOffset(s,c,m,vtx,bin);
	  if(i==0)  hvtx->Fill(vtx);	  
	  if(centval<=80&&abs(vtx)<plotZrange) hPsi[i]->Fill(angorig);
	  if(genFlatPsi_) {
	    if(centval<=80&&abs(vtx)<plotZrange) {
	      hPsiFlat[i]->Fill(psiFlat);
	      hPsiOffset[i]->Fill(psiOffset);
	      hPsiPsiFlat[i]->Fill(angorig,psiFlat);
	      hPsiPsiOffset[i]->Fill(angorig,psiOffset);
	    }
	    for(int j = 0; j<NumCentBins; j++) {
	      if(centval>wcent[j]&&centval<=wcent[j+1]&&abs(vtx)<plotZrange) hPsiFlatCent[i][j]->Fill(psiFlat);
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
    for(int j = 0; j<flat[i]->GetOBins();j++) {
      flatXhist[i]->SetBinContent(j+1+Hbins,flat[i]->GetXoff(j));
      flatYhist[i]->SetBinContent(j+1+Hbins,flat[i]->GetYoff(j));
      flatCnthist[i]->SetBinContent(j+1+Hbins,flat[i]->GetXYoffcnt(j));
      xoffhist[i]->SetBinContent(j+1,flat[i]->GetXoff(j));
      yoffhist[i]->SetBinContent(j+1,flat[i]->GetYoff(j));
      xyoffcnthist[i]->SetBinContent(j+1,flat[i]->GetXYoffcnt(j));
      xyoffmulthist[i]->SetBinContent(j+1,flat[i]->GetXYoffmult(j));

    }
    for(int j = 0; j<flat[i]->GetOBins();j++) {
      flatXhist[i]->SetBinContent(j+1+Hbins+Obins,flat[i]->GetPt(j));
      flatYhist[i]->SetBinContent(j+1+Hbins+Obins,flat[i]->GetPt2(j));
      flatCnthist[i]->SetBinContent(j+1+Hbins+Obins,flat[i]->GetPtcnt(j));
      pthist[i]->SetBinContent(j+1,flat[i]->GetPt(j));
      pt2hist[i]->SetBinContent(j+1,flat[i]->GetPt2(j));
      ptcnthist[i]->SetBinContent(j+1,flat[i]->GetPtcnt(j));
    }
  }
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiEvtPlaneFlatCalib);
