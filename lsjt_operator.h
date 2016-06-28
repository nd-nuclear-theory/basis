/****************************************************************
  lsjt_operator.h

  Defines functions for I/O and manipulation of relative and two-body
  interaction matrices in LSJT coupling scheme.  Written for use in
  Moshinsky transformation.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  11/21/15 (mac): Created (lsjt_interaction.h), from code in
    moshinsky_xform_lsjt.
  6/8/16 (mac): Update for basis package and update conventions.
    - Rename to lsjt_operator.h.
    - Remove matrix-style output.
    - Remove generic template functions.
  6/27/16 (mac): work in progress...

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
  // relative basis
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // general relative two-body operator format
  //
  // RELATIVE version    # version = 1
  // J0 g0               # operator a.m. and parity [P0=(-)^g0]
  // T0_min T0_max       # operator isospin component (a TBO may have components T0=0,1,2)
  // Nr_max Jr_max       # basis: max Nr (relative quanta), max Jr for subspaces (<=Nr_max+1)
  //
  // Then data lines are of form:
  //
  //    T0 Nr' lr' S' Jr' T' Nr lr S Jr T matrix_element
  //
  // Here matrix_element is the JT-reduced matrix element under group
  // theory conventions (i.e., no dimension factor in the
  // Wigner-Eckart theorem):
  //
  //    < Nr' lr' S' Jr' T' || op || Nr lr S Jr T >
  //
  // For the special case of a rotational scalar (J0=0), isoscalar
  // (T0=0) operator, this is equivalently the unreduced matrix
  // element, as is often used to represent Hamiltonians in shell
  // model codes:
  //
  //    < Nr' lr' S' Jr' T'; MJr MT | op | Nr lr S Jr T; MJr MT >
  //
  // (Note that this matrix element is independent of MJr and MT for a
  // scalar, isoscalar operator.)
  //
  // Iteration order:
  //
  //   It is assumed that the RelativeSectorsLSJT was constucted with
  //   the direction=kCanonical option.
  //
  //   Then ordering is...
  //
  //     for each T0 component (T0=T0_min...T0_max)
  //       for each sector
  //         for each matrix element
  //
  //   Sectors are thus ordered lexicographically by
  //   (bra_subspace_index,ket_subspace_index) and subject to these
  //   indices (bra_subspace_index,ket_subspace_index) being in
  //   canonical order.  It is assumed that the matrix elements in the
  //   other direction can be obtained by symmetry (CAVEAT: phase
  //   factors TBD).
  //
  //   States are likewise ordered lexicographical by
  //   (bra_state_index,ket_state_index), i.e., row-major ordering.
  //   In diagonal sectors (i.e., when the bra and ket subspaces are
  //   the same), only matrix elements with indices in canonical order
  //   (upper triangle) are written to file.


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
