/****************************************************************
  many_body.h

  Provides convenience definitions for many-body indexing schemes.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/19/16 (mac): Created, incorporating code from operator.h.

****************************************************************/

#ifndef MANY_BODY_H_
#define MANY_BODY_H_

#include <cassert>
#include <tuple>

namespace basis {

  ////////////////////////////////////////////////////////////////
  // many-body truncation rank
  ////////////////////////////////////////////////////////////////

  enum class Rank {kOneBody=1, kTwoBody=2, kThreeBody=3, kFourBody=4};
  // Truncation rank for many-body bases.
  //
  // To be used as a mode parameter for basis constructors, to
  // determine the many-body rank of the truncation scheme, e.g.,
  // one-body ("square") truncation, two-body ("triangular")
  // truncation.

  std::tuple<int,int> TwoBodyCutoffs(basis::Rank truncation_rank, int truncation_cutoff);
  // Calculate one-body N1max and two-body N2max cutoffs for truncation scheme.
  //
  // Example:
  //    int N1max, N2max;
  //    std::tie(N1max,N2max) = basis::TwoBodyCutoffs(truncation_rank,truncation_cutoff);

  //
  // Arguments:
  //   truncation_rank (basis::Rank) : many-body rank (kOneBody or kTwoBody)
  //   truncation_cutoff (int) : oscillator cutoff for that rank
  //
  // Returns:
  //   (std::tuple<int,int>) : one-body and two-body cutoffs (N1max,N2max)
  

  ////////////////////////////////////////////////////////////////
  // many-body normalization conversion flag
  ////////////////////////////////////////////////////////////////

  enum class NormalizationConversion {kNone, kASToNAS, kNASToAS};
  // Normalization scheme converson flag.
  //
  // To be used as value for a mode parameter to requests on-the-fly
  // conversion between antisymmetrized (AS) and normalized
  // antisymmetrized (NAS) matrix elements on input/output.



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
