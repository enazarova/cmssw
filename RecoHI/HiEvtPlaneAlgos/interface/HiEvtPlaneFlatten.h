#ifndef __HiEvtPlaneFlatten__
#define __HiEvtPlaneFlatten__
// -*- C++ -*-
//
// Package:    HiEvtPlaneFlatten
// Class:      HiEvtPlaneFlatten
// 

//
//
// Original Author:  Stephen Sanders
//         Created:  Mon Jun  7 14:40:12 EDT 2010
// $Id: HiEvtPlaneFlatten.h,v 1.4 2011/11/06 23:17:27 ssanders Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

#include "TMath.h"
#include <vector>

#define MAXCUT 10000
#define MAXCUTOFF 1000

//
// class declaration
//

class HiEvtPlaneFlatten {
public:

  explicit HiEvtPlaneFlatten()
  {
    pi = TMath::Pi();
    hcentbins = 1;
    hOrder = 9;
    vorder = 2;    //sets default order of event plane
    minvtx = -25;
    delvtx = 5;
    nvtxbins = 10;
  }
  void Init(int order, int ncentbins,const int centbinCompression, std::string tag, int vord)
  {
    hOrder = order;  //order of flattening
    vorder = vord;   //1(v1), 2(v2), 3(v3), 4(v4)	
    hcentbins = ncentbins;
    centbinComp = centbinCompression;
    if(hcentbins<=0) hcentbins = 1;
    hbins = hcentbins*nvtxbins*hOrder;
    obins = hcentbins*nvtxbins;
    if(hbins>MAXCUT) {
      std::cout<<"Too many cuts for flattening calculation.  RESET to deaults"<<std::endl;
      hcentbins = 1;
      hOrder = 9;
    }
    for(int i = 0; i<hbins; i++) {
      flatX[i]=0;
      flatY[i]=0;
      flatXDB[i]=0;
      flatYDB[i]=0;
      flatCnt[i]=0;
    } 
    for(int i = 0; i<obins; i++) {
      xoff[i]=0;
      yoff[i]=0;
      xoffDB[i]=0;
      yoffDB[i]=0;
      xyoffcnt[i]=0;
      xyoffmult[i]=0;
      pt[i]=0;
      pt2[i]=0;
      ptDB[i]=0;
      pt2DB[i]=0;
      ptcnt[i]=0;
    }
  }

  int GetCutIndx(int centbin, double vtx, int iord)
  {
    int cut;
    int icent = centbin/centbinComp;
    if(icent < 0 || icent > hcentbins) return -1;
    int ivtx = (vtx-minvtx)/delvtx;
    if(ivtx < 0 || ivtx > nvtxbins) return -1;
    cut = hOrder*nvtxbins*icent + hOrder*ivtx + iord;
    if(cut<0 || cut>hbins) return -1;
    return cut;
  }

  int GetOffsetIndx(int centbin, double vtx)
  {
    int cut;
    int icent = centbin/centbinComp;
    if(icent < 0 || icent > hcentbins) return -1;
    int ivtx = (vtx-minvtx)/delvtx;
    if(ivtx < 0 || ivtx > nvtxbins) return -1;
    cut = nvtxbins*icent + ivtx ;
    if(cut<0 || cut>hbins) return -1;
    return cut;
  }
  
  void Fill(double psi, double vtx, int centbin)
  {
    if(fabs(psi)>4 ) return;
    for(int k = 0; k<hOrder; k++) {
      double fsin = sin(vorder*(k+1)*psi);
      double fcos = cos(vorder*(k+1)*psi);
      int indx = GetCutIndx(centbin,vtx,k);
      if(indx>=0) {
	flatX[indx]+=fcos;
	flatY[indx]+=fsin;
	++flatCnt[indx];
      }
    }
  }
  void FillOffset(double s, double c, uint m, double vtx, int centbin)
  {
    int indx = GetOffsetIndx(centbin,vtx);
    if(indx>=0) {
      xoff[indx]+=c;
      yoff[indx]+=s;
      xyoffmult[indx]+=m;
      ++xyoffcnt[indx];
    }
  }
  void FillPt(double ptval, double vtx, int centbin)
  {
  
    int indx = GetOffsetIndx(centbin,vtx);
    if(indx>=0) {
      pt[indx]+=ptval;
      pt2[indx]+=ptval*ptval;
      ++ptcnt[indx];
    }
  }

  double GetW(double pt, double vtx, int centbin)
  {
  
    int indx = GetOffsetIndx(centbin,vtx);
    if(indx>=0) {
      double ptval = GetPtDB(indx);
      double pt2val = GetPt2DB(indx);
      if(ptval>0) return pt-pt2val/ptval;
    }
    return 0.;
  }

  double GetFlatPsi(double psi, double vtx, double cent)
  {
    double correction = 0;
    for(int k = 0; k<hOrder; k++) {
      int indx = GetCutIndx(cent,vtx,k);
      if(indx>=0) correction+=(2./(double)((k+1)*vorder))*(flatXDB[indx]*sin(vorder*(k+1)*psi)-flatYDB[indx]*cos(vorder*(k+1)*psi));
    }
    psi+=correction;
    psi=bounds(psi);
    psi=bounds2(psi);
    return psi;
  }
  
  double GetOffsetPsi(double s, double c,  double vtx, double cent)
  {
    int indx = GetOffsetIndx(cent,vtx);
    double snew = s-yoffDB[indx];
    double cnew = c-xoffDB[indx];
    double psi = atan2(snew,cnew)/vorder;
    if((fabs(snew)<1e-4) && (fabs(cnew)<1e-4)) psi = 0.;
    psi=bounds(psi);
    psi=bounds2(psi);
    return psi;
  }
  
  ~HiEvtPlaneFlatten(){}
  int GetHBins(){return hbins;}
  int GetOBins(){return obins;}
  int GetNvtx(){return nvtxbins;}
  double GetVtxMin(){return minvtx;}
  double GetVtxMax(){return minvtx+nvtxbins*delvtx;}
  int GetNcent(){return hcentbins;}

  double GetX(int bin){return flatX[bin];}
  double GetY(int bin){return flatY[bin];}
  double GetXoff(int bin){return xoff[bin];}
  double GetYoff(int bin){return yoff[bin];}
  double GetXoffDB(int bin){return xoffDB[bin];}
  double GetYoffDB(int bin){return yoffDB[bin];}
  double GetXYoffcnt(int bin){return xyoffcnt[bin];}
  double GetXYoffmult(int bin){return xyoffmult[bin];}
  double GetPt(int bin){return pt[bin];}
  double GetPt2(int bin){return pt2[bin];}
  double GetPtDB(int bin){return ptDB[bin];}
  double GetPt2DB(int bin){return pt2DB[bin];}
  double GetPtcnt(int bin){return ptcnt[bin];}


  double GetCnt(int bin) {return flatCnt[bin];}
  void SetXDB(int indx, double val) {flatXDB[indx]=val;}
  void SetYDB(int indx, double val) {flatYDB[indx]=val;}
  void SetXoffDB(int indx, double val) {xoffDB[indx]=val;}
  void SetYoffDB(int indx, double val) {yoffDB[indx]=val;}
  void SetPtDB(int indx, double val) {ptDB[indx]=val;}
  void SetPt2DB(int indx, double val) {pt2DB[indx]=val;}
  Double_t bounds(Double_t ang) {
    if(ang<-pi) ang+=2.*pi;
    if(ang>pi)  ang-=2.*pi;
    return ang;
  }
  Double_t bounds2(Double_t ang) {
    double range = pi/(double) vorder;
    if(ang<-range) ang+=2*range;
    if(ang>range)  ang-=2*range;
    return ang;
  }
private:
  double flatX[MAXCUT];
  double flatY[MAXCUT];
  double flatXDB[MAXCUT];
  double flatYDB[MAXCUT];
  double flatCnt[MAXCUT];



  double xoff[MAXCUTOFF];
  double yoff[MAXCUTOFF];
  double xoffDB[MAXCUTOFF];
  double yoffDB[MAXCUTOFF];
  double xyoffcnt[MAXCUTOFF];
  uint xyoffmult[MAXCUTOFF]; 

  double pt[MAXCUTOFF];
  double pt2[MAXCUTOFF];
  double ptDB[MAXCUTOFF];
  double pt2DB[MAXCUTOFF];
  double ptcnt[MAXCUTOFF];

  int hOrder;    //flattening order
  int hcentbins; //# of centrality bins
  int centbinComp;
  int hbins; //number of bins needed for flattening
  int obins; //number of (x,y) offset bins
  int vorder; //order of flattened event plane
  double pi;

  int nvtxbins;
  double minvtx;
  double delvtx;

};



#endif
