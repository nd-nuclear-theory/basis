/****************************************************************
  lsjt_operator.h

  Defines functions for I/O and manipulation of relative and two-body
  interaction matrices in LSJT coupling scheme.  Written for use in
  Moshinsky transformation.

  Language: C++11
  (for tuple and map::at)
                                 
  Mark A. Caprio, University of Notre Dame.

  11/21/15 (mac): Created (lsjt_interaction.h), from code in
    moshinsky_xform_lsjt.
  6/8/16 (mac): Update for basis package and update conventions.
    - Rename to lsjt_operator.h.
    - Remove matrix-style output.
    - Remove generic template functions.

****************************************************************/

#ifndef LSJT_OPERATOR_H_
#define LSJT_OPERATOR_H_

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include "eigen3/Eigen/Core"

#include "am/am.h"
#include "basis/lsjt_scheme.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // common typedefs
  ////////////////////////////////////////////////////////////////

  typedef std::vector< Eigen::MatrixXd > MatrixVector;

  ////////////////////////////////////////////////////////////////
  // relative space
  ////////////////////////////////////////////////////////////////

  // relative space output
  void WriteRelativeOperatorLSJT(
      std::ostream& os,
      const RelativeSectorsLSJT& sectors,
      const MatrixVector& matrices
    );

  // relative space input
  void ReadRelativeOperatorLSJT(
      std::istream& is,
      const RelativeSectorsLSJT& sectors,
      MatrixVector& matrices
    );

  ////////////////////////////////////////////////////////////////
  // two-body space
  ////////////////////////////////////////////////////////////////

  void WriteTwoBodyOperatorLSJT(
      std::ostream& os,
      const TwoBodySectorsLSJT& sectors,
      const MatrixVector& matrices
    );


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
