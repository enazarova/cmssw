// //
// // $Id: EvtPlane.cc,v 1.4 2009/09/08 12:33:12 edwenger Exp $
// //
// 
// #include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
// 
// using namespace reco;
// 
// EvtPlane::EvtPlane(double planeA, std::string label)
//   : 
//    angle_(planeA),
//    label_(label)
// {
//   // default constructor
// }
// 
// 
// EvtPlane::~EvtPlane()
// {
// }
// 
// 


//
// $Id: EvtPlane.cc,v 1.4 2009/09/08 12:33:12 edwenger Exp $
//

#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

using namespace reco;
using namespace std;
EvtPlane::EvtPlane(double planeA,double sumSin, double sumCos, double sumw, uint mult, std::string label)
  : 
   angle_(planeA),
   sumSin_(sumSin),
   sumCos_(sumCos),
   sumw_(sumw),
   mult_(mult),
   label_(label)
{
  // default constructor
}


EvtPlane::~EvtPlane()
{
}


