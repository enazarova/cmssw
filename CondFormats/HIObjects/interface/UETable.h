#ifndef __UETable_h__
#define __UETable_h__

#include "CondFormats/Serialization/interface/Serializable.h"
#include <vector>

class UETable{
 public:
  UETable(){};
  std::vector<std::vector<std::vector<std::vector<std::vector<float> > > > > values;

  COND_SERIALIZABLE;
};

#endif
