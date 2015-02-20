
//
// $Id: EvtPlane.h,v 1.4 2009/09/08 12:33:11 edwenger Exp $
//

#ifndef DataFormats_EvtPlane_h
#define DataFormats_EvtPlane_h

#include <vector>
#include <string>
#include <math.h>

namespace reco { class EvtPlane {
  public:
    EvtPlane(int epindx=0, int level = 2, double planeA=0,double sumSin=0, double sumCos=0, double sumw = 0, double sumw2 = 0, double pe = 0, double pe2 = 0, uint mult = 0);
    virtual ~EvtPlane();
    void AddLevel(int level, double ang, double sumsin, double sumcos);
    int indx() const { return indx_ ;}
    double      angle(int level=2)   const { return (level>=0||level<=2)? angle_[level]:angle_[2]; }
    double      sumSin(int level=2)  const { return (level>=0||level<=2)? sumSin_[level]:sumSin_[2];}
    double      sumCos(int level=2)  const { return (level>=0||level<=2)? sumCos_[level]:sumCos_[2];}
    double      sumw() const { return sumw_;}
    double      sumw2() const { return sumw2_;}
    double      sumPtOrEt() const { return sumPtOrEt_;}
    double      sumPtOrEt2() const { return sumPtOrEt2_;}
    double      mult()    const { return mult_;}
    double      Qy(int level=2) const { return sumSin(level); }
    double      Qx(int level=2) const { return sumCos(level); }
    double      Q(int level=2)      const { return ((pow(Qx(level),2)+pow(Qy(level),2))>0)? sqrt(pow(Qx(level),2)+pow(Qy(level),2)): 0.;}
    double      qy(int level=2)      const { return (mult_>0)? ((level>=0||level<=2)? sumSin_[level]/sqrt((double)mult_):sumSin_[2]/sqrt((double) mult_)):0.;}
    double      qx(int level=2)      const { return (mult_>0)? ((level>=0||level<=2)? sumCos_[level]/sqrt((double)mult_):sumCos_[2]/sqrt((double) mult_)):0.;}
    double      q(int level=2)      const { return ((pow(qx(level),2)+pow(qy(level),2))>0)? sqrt(pow(qx(level),2)+pow(qy(level),2)): 0.;}
        
  private:
    int           indx_;
    double        angle_[4];
    double        sumSin_[4];
    double        sumCos_[4];
    double        sumw_;
    double        sumw2_;
    double        sumPtOrEt_;
    double        sumPtOrEt2_;
    uint          mult_;
    
  };
  
  typedef std::vector<EvtPlane> EvtPlaneCollection;
  
}

#endif 






