#ifndef __RPFlatParams_h__
#define __RPFlatParams_h__

#include <vector>

class RPFlatParams{
 public:
  static const int MaxEPAllowed = 30;
  struct EP {
    float x[MaxEPAllowed];
    float y[MaxEPAllowed];
    int RPNameIndx[MaxEPAllowed];
  };
  RPFlatParams(){}
  virtual ~RPFlatParams(){}
  std::vector<EP> m_table;
};

#endif

