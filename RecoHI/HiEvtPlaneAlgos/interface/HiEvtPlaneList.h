#ifndef __HiEvtPlaneList__
#define __HiEvtPlaneList__
/*
Indx       Name Order EtaMin1 EtaMax1 EtaMin2 EtaMax2 ResCalcType MomConsWeight ResCoreMate1 ResCoreMate2
   0      HFm2      2    -5.0    -3.0     0.0     0.0       3sub          no          HFp2     trackmid2
   1      HFp2      2     3.0     5.0     0.0     0.0       3sub          no          HFm2     trackmid2
   2       HF2      2    -5.0    -3.0     3.0     5.0       3sub          no       trackm2       trackp2
   3 trackmid2      2    -0.7     0.7     0.0     0.0       3sub          no          HFm2          HFp2
   4   trackm2      2    -2.0    -1.0     0.0     0.0       3sub          no          HFm2          HFp2
   5   trackp2      2     1.0     2.0     0.0     0.0       3sub          no          HFm2          HFp2
   6      HFm3      3    -5.0    -3.0     0.0     0.0       3sub          no          HFp3     trackmid3
   7      HFp3      3     3.0     5.0     0.0     0.0       3sub          no          HFm3     trackmid3
   8       HF3      3    -5.0    -3.0     3.0     5.0       3sub          no       trackm3       trackp3
   9 trackmid3      3    -0.7     0.7     0.0     0.0       3sub          no          HFm3          HFp3
  10   trackm3      3    -2.0    -1.0     0.0     0.0       3sub          no          HFm3          HFp3
  11   trackp3      3     1.0     2.0     0.0     0.0       3sub          no          HFm3          HFp3
  12      HFm4      4    -5.0    -3.0     0.0     0.0       3sub          no          HFp4     trackmid4
  13      HFp4      4     3.0     5.0     0.0     0.0       3sub          no          HFm4     trackmid4
  14       HF4      4    -5.0    -3.0     3.0     5.0       3sub          no       trackm4       trackp4
  15 trackmid4      4    -0.7     0.7     0.0     0.0       3sub          no          HFm4          HFp4
  16   trackm4      4    -2.0    -1.0     0.0     0.0       3sub          no          HFm4          HFp4
  17   trackp4      4     1.0     2.0     0.0     0.0       3sub          no          HFm4          HFp4
  18      HFm1      1    -5.0    -3.0     0.0     0.0       2sub          no          HFp1              
  19      HFp1      1     3.0     5.0     0.0     0.0       2sub          no          HFm1              
  20       HF1      1    -5.0    -3.0     3.0     5.0       3sub          no       trackm1       trackp1
  21   trackm1      1    -2.2    -1.4     0.0     0.0       3sub         yes          HFm1          HFp1
  22   trackp1      1     1.4     2.2     0.0     0.0       3sub         yes          HFm1          HFp1
  23    HFm1mc      1    -5.0    -3.0     0.0     0.0       3sub         yes        HFp1mc       trackp1
  24    HFp1mc      1     3.0     5.0     0.0     0.0       3sub         yes        HFm1mc       trackm1
*/
#include <string>
using namespace std;
namespace hi{
enum EPNamesInd {
      HFm2,	      HFp2,	       HF2,	 trackmid2,	   trackm2,	
   trackp2,	      HFm3,	      HFp3,	       HF3,	 trackmid3,	
   trackm3,	   trackp3,	      HFm4,	      HFp4,	       HF4,	
 trackmid4,	   trackm4,	   trackp4,	      HFm1,	      HFp1,	
       HF1,	   trackm1,	   trackp1,	    HFm1mc,	    HFp1mc,	   EPBLANK
};
const int RCMate1[] = {
      HFp2,	      HFm2,	   trackm2,	      HFm2,	      HFm2,	
      HFm2,	      HFp3,	      HFm3,	   trackm3,	      HFm3,	
      HFm3,	      HFm3,	      HFp4,	      HFm4,	   trackm4,	
      HFm4,	      HFm4,	      HFm4,	      HFp1,	      HFm1,	
   trackm1,	      HFm1,	      HFm1,	    HFp1mc,	    HFm1mc
};
const int RCMate2[] = {
 trackmid2,	 trackmid2,	   trackp2,	      HFp2,	      HFp2,	
      HFp2,	 trackmid3,	 trackmid3,	   trackp3,	      HFp3,	
      HFp3,	      HFp3,	 trackmid4,	 trackmid4,	   trackp4,	
      HFp4,	      HFp4,	      HFp4,	EPBLANK,	EPBLANK,	
   trackp1,	      HFp1,	      HFp1,	   trackp1,	   trackm1
};
const std::string EPNames[] = {
"HFm2",	"HFp2",	"HF2",	"trackmid2",	"trackm2",	
"trackp2",	"HFm3",	"HFp3",	"HF3",	"trackmid3",	
"trackm3",	"trackp3",	"HFm4",	"HFp4",	"HF4",	
"trackmid4",	"trackm4",	"trackp4",	"HFm1",	"HFp1",	
"HF1",	"trackm1",	"trackp1",	"HFm1mc",	"HFp1mc"
};
enum Detectors {Tracker, HF};
const int EPDet[] = {
HF,	HF,	HF,	Tracker,	Tracker,	
Tracker,	HF,	HF,	HF,	Tracker,	
Tracker,	Tracker,	HF,	HF,	HF,	
Tracker,	Tracker,	Tracker,	HF,	HF,	
HF,	Tracker,	Tracker,	HF,	HF
};
const int EPOrder[] = {
2,	2,	2,	2,	2,	
2,	3,	3,	3,	3,	
3,	3,	4,	4,	4,	
4,	4,	4,	1,	1,	
1,	1,	1,	1,	1
};
const double EPEtaMin1[] = {
-5.00,	 3.00,	-5.00,	-0.70,	-2.00,	
 1.00,	-5.00,	 3.00,	-5.00,	-0.70,	
-2.00,	 1.00,	-5.00,	 3.00,	-5.00,	
-0.70,	-2.00,	 1.00,	-5.00,	 3.00,	
-5.00,	-2.20,	 1.40,	-5.00,	 3.00
};
const double EPEtaMax1[] = {
-3.00,	 5.00,	-3.00,	 0.70,	-1.00,	
 2.00,	-3.00,	 5.00,	-3.00,	 0.70,	
-1.00,	 2.00,	-3.00,	 5.00,	-3.00,	
 0.70,	-1.00,	 2.00,	-3.00,	 5.00,	
-3.00,	-1.40,	 2.20,	-3.00,	 5.00
};
const double EPEtaMin2[] = {
 0.00,	 0.00,	 3.00,	 0.00,	 0.00,	
 0.00,	 0.00,	 0.00,	 3.00,	 0.00,	
 0.00,	 0.00,	 0.00,	 0.00,	 3.00,	
 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	
 3.00,	 0.00,	 0.00,	 0.00,	 0.00
};
const double EPEtaMax2[] = {
 0.00,	 0.00,	 5.00,	 0.00,	 0.00,	
 0.00,	 0.00,	 0.00,	 5.00,	 0.00,	
 0.00,	 0.00,	 0.00,	 0.00,	 5.00,	
 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	
 5.00,	 0.00,	 0.00,	 0.00,	 0.00
};
const std::string ResCalcType[] = {
"3sub",	"3sub",	"3sub",	"3sub",	"3sub",	
"3sub",	"3sub",	"3sub",	"3sub",	"3sub",	
"3sub",	"3sub",	"3sub",	"3sub",	"3sub",	
"3sub",	"3sub",	"3sub",	"2sub",	"2sub",	
"3sub",	"3sub",	"3sub",	"3sub",	"3sub"
};
const std::string MomConsWeight[] = {
"no",	"no",	"no",	"no",	"no",	
"no",	"no",	"no",	"no",	"no",	
"no",	"no",	"no",	"no",	"no",	
"no",	"no",	"no",	"no",	"no",	
"no",	"yes",	"yes",	"yes",	"yes"
};
static const int NumEPNames = 25;
static const int FlatOrder = 9;
static const int NumFlatCentBins = 50; //bins AFTER compression
//CentBinCompression*NumFlatCentBins = #bins from centrality module
static const int CentBinCompression = 4;
}
#endif
