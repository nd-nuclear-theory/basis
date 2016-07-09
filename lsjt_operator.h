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
  7/6/16 (mac):
    - Add symmetry phase header field and update documentation.
    - Upgrade precision on operator output.
  7/8/16 (mac): Add two-body LSJT operator file I/O.
  7/9/16 (mac): Add support for canonicalizing indices in relative
    LSJT matrix element lookup (for Hamiltonian-like operators).

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
  // Header
  //
  //   Header format:
  //
  //     # RELATIVE LSJT
  //     # ...
  //     version                               # version
  //     J0 g0 T0_min T0_max symmetry_phase    # operator tensor properties
  //     Nmax Jmax                             # relative basis truncation
  //
  //   The header may start with one or more contiguous comment lines, 
  //   which are designated by a hash character in the first column.
  //  
  //   Then, the header contains the following fields:
  //
  //     version : file format version (=1)
  //     J0 (int) : angular momentum of operator
  //     g0 (int) : parity grade of operator [P0=(-)^g0] (g0=0,1)
  //     T0_min T0_max (int) : range of operator isospin components (a two-body operator  
  //       may in general have components T0=0,1,2 from the coupling of four isospin-1/2 
  //       fermionic operators)
  //     symmetry_phase (int) : RESERVED to describe how to obtain phase for lower triangle
  //       (see "Conjugation symmetry" below) (=0)
  //     Nmax (int) : oscillator truncation of relative space (Nmax>=0)
  //     Jmax (int) : additional relative angular momentum truncation of relative space
  //       (Jmax<=Nmax+1); Jmax=Nmax+1 includes the full relative space at the given 
  //       oscillator truncation, but operators may often be truncated at lower partial
  //       waves
  //
  // Data
  //
  //   Then data lines are of the form:
  //
  //      T0 N' L' S' J' T' N L S J T JT-RME
  //
  //   Here JT-RME is the JT-reduced matrix element under group
  //   theory conventions (i.e., no dimension factor in the
  //   Wigner-Eckart theorem):
  //
  //      < N' L' S' J' T' || op || N L S J T >
  //
  //   For the special case of a rotational scalar (J0=0), isoscalar
  //   (T0=0) operator, this is equivalently the unreduced matrix
  //   element, as is often used to represent Hamiltonians in shell
  //   model codes:
  //
  //      < N' L' S' J' T'; MJ MT | op | N L S J T; MJ MT >
  //
  //   (Note that this matrix element is independent of MJ and MT for a
  //   scalar, isoscalar operator.)
  //
  // Iteration order and symmetry
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
  //   States are likewise ordered lexicographical by
  //   (bra_state_index,ket_state_index), i.e., row-major ordering.
  //   In diagonal sectors (i.e., when the bra and ket subspaces are
  //   the same), only matrix elements with indices in canonical order
  //   (upper triangle) are read from or written to file.
  //
  // Conjugation symmetry
  //
  //   Since only the upper triangle of the operator is stored, the
  //   user code is responsible for obtaining the lower triangle
  //   (conjugate) matrix elements by "symmetry".  Recall that the
  //   matrix elements are JT-reduced matrix elements, not simple
  //   matrix elements.  Therefore, conjuation will in general involve
  //   phase and dimension factors.
  //
  //   The symmetry_phase field in the header is reserved to provide
  //   information on the correct form to use for this phase factor.
  //   However, for now, only the placeholder value kHermitian(=0) is
  //   defined for symmetry_phase, and phase conventions are only
  //   well-defined for a Hamiltonian-like (J0=0, g0=0) operator.
  //
  //   symmetry_phase=kHermitian, J0=0, g0=0: For a Hamiltonian-like operator,
  //   we expect Hermiticity, i.e., symmetry of (M_J,M_T)-branched
  //   matrix elements.  Within a diagonal sector, this means that the
  //   lower triangle is obtained from the upper triangle by ordinary
  //   symmetry.  For off-diagonal sectors, the appropriate symmetry
  //   on the JT-reduced matrix elements in general includes phase and
  //   dimension factors from the isospin Clebsches:
  //
  //     <a,J,T,g || A_{T0} || a',J,T',g>
  //       = (-)^(T'-T)*Hat(T')/Hat(T)
  //         * <a',J,T',g || A_{T0} || a,J,T,g>
  //
  //   However, note that the symmetry factor only differs from unity
  //   in the isospin-changing <T=0|T=1> sectors.  These sectors
  //   conventionally only contain vanishing matrix elements for
  //   nuclear interactions.
  //
  // Radial oscillator phase convention
  //
  //   Two phase conventions are in use for radial wave functions and
  //   thus for harmonic oscillator basis functions.  See, e.g.,
  //   footnote 8 of csbasis [PRC 86, 034312 (2012)].  We adopt the
  //   "positive at the origin" convention, which should thus be used
  //   when writing matrix elements to file in the present format.
  //   The "positive at the origin" convention is followed for the
  //   oscillator functions in, e.g., Suhonen (3.42), Moshinsky
  //   (1.1.8) & (1.9.15), and MFDn interaction files.  However, the
  //   "positive at infinity" convention is more natural when
  //   considering oscillator functions as members of SU(3) irreps.
  //   So it may be necessary to convert to this convention for
  //   internal use in SU(3) calculations.
  
  ////////////////////////////////////////////////////////////////
  // operator header
  ////////////////////////////////////////////////////////////////

  enum class SymmetryPhase {kHermitian=0};

  struct RelativeOperatorParametersLSJT
  // Parameters for relative operator storage.
  //
  // Data members:
  //   See comment describing header for relative operator file.
  {
    int version;
    int J0, g0, T0_min, T0_max;
    basis::SymmetryPhase symmetry_phase;
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
  //   is (std::istream) : text-mode input stream
  //   parameters (basis::RelativeOperatorParametersLSJT, output)
  //     : operator parameters

  void WriteRelativeOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const basis::RelativeSectorsLSJT& sectors,
      const basis::MatrixVector& matrices
    );
  // Write single isospin component of a relative operator in LSJT
  // scheme.
  //
  // Side effect: The floating point precision attribute of the output
  // stream is modified.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   T0 (int) : isospin for this isospin component
  //   sector (basis::RelativeSectorsLSJT) :  sectors defining operator
  //   matrices (basis::MatrixVector) :  matrices defining operator

  void ReadRelativeOperatorComponentLSJT(
      std::istream& is,
      int T0,
      const basis::RelativeSectorsLSJT& sectors,
      basis::MatrixVector& matrices
    );
  // Read single isospin component of a relative operator in LSJT
  // scheme.
  //
  // Arguments:
  //   is (std::istream) : text-mode inpus stream
  //   T0 (int) : isospin for this isospin component
  //   sector (basis::RelativeSectorsLSJT) :  sectors defining operator
  //   matrices (basis::MatrixVector, output) :  matrices defining operator

  ////////////////////////////////////////////////////////////////
  // operator matrix element canonicalizaton
  ////////////////////////////////////////////////////////////////

  void CanonicalizeIndicesRelativeLSJT(
      const basis::RelativeSpaceLSJT& space,
      int& bra_subspace_index, int& ket_subspace_index,
      int& bra_state_index, int& ket_state_index,
      double& canonicalization_factor,
      int J0, int T0, int g0,
      basis::SymmetryPhase symmetry_phase
    );
  // Convert subspace and state indices for a matrix element to
  // canonical ("upper triangle") indices.
  //
  // This is a customized wrapper for basis::CanonicalizeIndices (see
  // operator.h), for use with RelativeLSJT operators.
  //
  // Arguments:
  //    relative_space (basis::RelativeSpaceLSJT) : space, for retrieving 
  //      subspace quantum numbers to calculate canonicalization factor
  //    bra_subspace_index, ket_subspace_index (int, input/output) :
  //      sector bra and ket subspace indices, possibly to be swapped
  //    bra_state_index, ket_state_index (int, input/output) :
  //      bra and ket state indices, possibly to be swapped if sector
  //      is diagonal sector
  //    canonicalization_factor (double, output) : phase and dimension
  //      factor arising from any swaps
  //    J0, T0, g0 (int) : operator properties
  //    symmetry_phase (basis::SymmetryPhase) : specification of
  //      matrix element conjugation properties of the operator



  ////////////////////////////////////////////////////////////////
  // two-body LSJT operator
  ////////////////////////////////////////////////////////////////

  // Note that the primary intention of the output for two-body
  // operators in LSJT scheme is for diagnostic purposes.  If we were
  // to store operators more permanently in this format, we would also
  // want to define an appropriate file header format.
  //
  // Data lines are of the form:
  //
  //   T0  N1' l1' N2' l2' L' S' J' T' g'  N1 l1 N2 l2 L S J T g  JT-RME
  //
  // Although the g label is redundant (it can be deduced from l1 and
  // l2), it is included to make the sector structure more easily
  // apparent to a human reader.
  //
  // Iteration follows the usual scheme within the basis module:
  // sectors are lexicographic by (bra,ket) subspace indices, then
  // matrix elements within a sector are lexicographic by (bra,ket)
  // state indices.
  //
  // Reminder: One should be sure to document whether one is writing
  // AS or NAS matrix elements!

  void WriteTwoBodyOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const basis::TwoBodySectorsLSJT& sectors,
      const basis::MatrixVector& matrices,
      basis::NormalizationConversion conversion_mode
    );
  // Write single isospin component of a two-body operator in LSJT
  // scheme.
  //
  // Side effect: The floating point precision attribute of the output
  // stream is modified.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   T0 (int) : isospin for this isospin component
  //   sector (basis::TwoBodySectorsLSJT) :  sectors defining operator
  //   matrices (basis::MatrixVector) :  matrices defining operator
  //   conversion (basis::NormalizationConversion) : specifies any 
  //     conversion between AS and NAS for output

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
