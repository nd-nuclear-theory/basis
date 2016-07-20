/****************************************************************
  operator.cpp
                                 
  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // storage diagnostic
  ////////////////////////////////////////////////////////////////
  std::size_t AllocatedEntries(const MatrixVector& matrices)
  {
    std::size_t count = 0;
    for (auto iterator = matrices.begin(); iterator != matrices.end(); ++iterator)
      count += iterator->size();
    return count;
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
