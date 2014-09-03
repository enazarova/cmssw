
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
    EvtPlane(double planeA=0,double sumSin=0, double sumCos=0, uint mult = 0, std::string label="");
  virtual ~EvtPlane();

  double      angle()   const { return angle_; }
  double      sumSin()  const { return sumSin_;}
  double      sumCos()  const { return sumCos_;}
  double      mult()    const { return mult_;}
  double      q() const {return ((pow(sumSin_,2)+pow(sumCos_,2))>0)? sqrt(pow(sumSin_,2)+pow(sumCos_,2)): 0.;};
  std::string label()   const { return label_; }

 

private:

  double        angle_  ;
  double        sumSin_;
  double        sumCos_;
  uint          mult_;
  std::string   label_;


};

 typedef std::vector<EvtPlane> EvtPlaneCollection;

}

#endif 






