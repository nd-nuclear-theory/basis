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
  7/3/16 (mac): Add relative LSJT operator file I/O.

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
#include "basis/operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  //
  // general relative two-body operator format
  //
  // # RELATIVE LSJT
  // version             # version = 1
  // J0 g0               # operator a.m. and parity [P0=(-)^g0]
  // T0_min T0_max       # operator isospin components (a TBO may have components T0=0,1,2)
  // Nmax Jmax           # basis: max N (relative quanta), max J (relative a.m.) (Jmax<=Nmax+1)
  //
  // Then data lines are of form:
  //
  //    T0 N' L' S' J' T' N L S J T matrix_element
  //
  // Here matrix_element is the JT-reduced matrix element under group
  // theory conventions (i.e., no dimension factor in the
  // Wigner-Eckart theorem):
  //
  //    < N' L' S' J' T' || op || N L S J T >
  //
  // For the special case of a rotational scalar (J0=0), isoscalar
  // (T0=0) operator, this is equivalently the unreduced matrix
  // element, as is often used to represent Hamiltonians in shell
  // model codes:
  //
  //    < N' L' S' J' T'; MJ MT | op | N L S J T; MJ MT >
  //
  // (Note that this matrix element is independent of MJ and MT for a
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
  //   other direction can be obtained by symmetry.
  //
  //   CAVEAT: The user code is currently responsible for "knowing"
  //   the phase factor required to obtain matrix elements by
  //   symmetry.
  //
  //   States are likewise ordered lexicographical by
  //   (bra_state_index,ket_state_index), i.e., row-major ordering.
  //   In diagonal sectors (i.e., when the bra and ket subspaces are
  //   the same), only matrix elements with indices in canonical order
  //   (upper triangle) are read from or written to file.
  //
  // Phase convention:
  //
  //   Two phase conventions are in use for radial wave functions and
  //   thus for harmonic oscillator basis functions.  See, e.g.,
  //   footnote 8 of csbasis [PRC 86, 034312 (2012)].  The "positive
  //   at the origin" convention is used for oscillator functions in,
  //   e.g., Suhonen (3.42), Moshinsky (1.1.8) & (1.9.15), and MFDn
  //   interaction files.  The "positive at the origin" convention
  //   should therefore be used when writing matrix elements to file
  //   in the present format.  However, the "positive at infinity"
  //   convention is more natural when considering oscillator
  //   functions as members of SU(3) irreps, so conversion may be
  //   necessary.
  
  ////////////////////////////////////////////////////////////////
  // operator header
  ////////////////////////////////////////////////////////////////

  struct RelativeOperatorParametersLSJT
  // Parameters for relative operator storage.
  //
  // Data members:
  //   See comment describing header for relative operator file.
  {
    int version;
    int J0, g0;
    int T0_min, T0_max;
    int Nmax, Jmax;
  };

  ////////////////////////////////////////////////////////////////
  // operator I/O
  ////////////////////////////////////////////////////////////////

  void WriteRelativeOperatorParametersLSJT(
      std::ostream& os,
      const basis::RelativeOperatorParametersLSJT& parameters
    );
  // Write file header for relative operator in LSJT scheme.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   parameters (basis::RelativeOperatorParametersLSJT)
  //     : operator parameters

  void ReadRelativeOperatorParametersLSJT(
      std::istream& is,
      basis::RelativeOperatorParametersLSJT& parameters
    );
  // Write file header for relative operator in LSJT scheme.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   parameters (basis::RelativeOperatorParametersLSJT, output)
  //     : operator parameters

  void WriteRelativeOperatorComponentLSJT(
      std::ostream& is,
      int T0,
      const basis::RelativeSectorsLSJT& sectors,
      const basis::MatrixVector& matrices
    );
  // Write single isospin component of a relative operator in LSJT
  // scheme.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   T0 (int) : isospin for this isospin component
  //   sector (basis::RelativeSectorsLSJT) :  sectors defining operator
  //   matrices (basis::MatrixVector) :  matrices defining operator

  void ReadRelativeOperatorComponentLSJT(
      std::istream& is,
      int T0,
      const RelativeSectorsLSJT& sectors,
      MatrixVector& matrices
    );
  // Read single isospin component of a relative operator in LSJT
  // scheme.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   T0 (int) : isospin for this isospin component
  //   sector (basis::RelativeSectorsLSJT) :  sectors defining operator
  //   matrices (basis::MatrixVector, output) :  matrices defining operator

  ////////////////////////////////////////////////////////////////
  // two-body LSJT operator
  ////////////////////////////////////////////////////////////////

  void WriteTwoBodyOperatorLSJT(
      std::ostream& os,
      const TwoBodySectorsLSJT& sectors,
      const MatrixVector& matrices
    );
  // TODO


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
