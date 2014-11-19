#ifndef ___EPFlattenParams_h__
#define ___EPFlattenParams_h__

#include <TNamed.h>
#include <TFile.h>
#include <vector>
#include <map>

class EPFlat : public TObject {
 public:
  static const int MaxEPAllowed = 30;
  EPFlat(){;}
  ~EPFlat(){;}
  float x[MaxEPAllowed];
  float y[MaxEPAllowed];
  int RPNameIndx[MaxEPAllowed];
  
  ClassDef(EPFlat,1)
};

class EPFlattenParams : public TNamed {
   
 public:
   typedef std::map<int, const EPFlattenParams*> RunMap;

   EPFlattenParams(){;}
   EPFlattenParams(const char* name, const char* title, int nbins) : TNamed(name,title) {
      table_.reserve(nbins);
      for(int j = 0; j < nbins; ++j){
	 EPFlat b;
	 table_.push_back(b); 
      }
   }
      ~EPFlattenParams() {;}
      // private:
      std::vector<EPFlat> table_;
      ClassDef(EPFlattenParams,1)
};

EPFlattenParams::RunMap getEPFlatParamsFromFile(TDirectoryFile*, const char* dir, const char* tag, int firstRun = 0, int lastRun = 10);
EPFlattenParams::RunMap getEPFlatParamsFromFile(TDirectoryFile*, const char* tag, int firstRun = 0, int lastRun = 10);





#endif
