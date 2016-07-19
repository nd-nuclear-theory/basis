/****************************************************************
  many_body.cpp

  Mark A. Caprio
  University of Notre Dame
****************************************************************/

#include "many_body.h"

namespace basis {

  std::tuple<int,int> TwoBodyCutoffs(basis::Rank truncation_rank, int truncation_cutoff)
  {
    // validate arguments
    assert(int(truncation_rank)<=2);
    assert(truncation_cutoff>=0);

    // extract one-body and two-body cutoffs
    int N1max = truncation_cutoff;
    int N2max;
    if (truncation_rank==basis::Rank::kOneBody)
      N2max = 2*truncation_cutoff;
    else if (truncation_rank==basis::Rank::kTwoBody)
      N2max = truncation_cutoff;
    
    // package return values
    return std::tuple<int,int>(N1max,N2max);
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
