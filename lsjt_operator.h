/****************************************************************
  @file lsjt_operator.h

  Defines functions for I/O and manipulation of relative and two-body
  operator matrices in LSJT coupling scheme.  Written for use in
  Moshinsky transformation.

  A few functions which provide simple "bookkeeping" operations on
  operators, such as ConstructIdentityOperatorRelativeLSJT or
  GatherOperatorTwoBodyLSJTNToTwoBodyLSJT are included here.  These
  provide examples of the coding idioms required to work with
  JT-coupled operators (looping over T0 components, looking up matrix
  elements, etc.).  But functions which have more mathematical content
  to them (Racah reduction formulas, Moshinsky transformation,
  recoupling, etc.) are deferred to the dedicated Moshinsky code.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 11/21/15 (mac): Created (lsjt_interaction.h), from code in
    moshinsky_xform_lsjt.
  + 06/08/16 (mac): Update for basis package and update conventions.
    - Rename to lsjt_operator.h.
    - Remove matrix-style output.
    - Remove generic template functions.
  + 07/03/16 (mac): Add relative LSJT operator file I/O.
  + 07/06/16 (mac):
    - Add symmetry phase header field and update documentation.
    - Upgrade precision on operator output.
  + 07/08/16 (mac): Add two-body LSJT operator file I/O.
  + 07/09/16 (mac): Add support for canonicalizing indices in relative
    LSJT matrix element lookup (for Hamiltonian-like operators).
  + 07/10/16 (mac):
    - Add support for canonicalizing indices in general
      LSJT matrix element lookup (for Hamiltonian-like operators).
    - Define OperatorLabelsJT and add documentation on operators.
    - Incorporate some basic LSJT operator construction and
      manipulation functions.
  + 07/13/16 (mac): Revise code for LSJTN->LSJT gathering operation.
  + 07/20/16 (mac): Add ReadRelativeOperatorLSJT.
  + 07/22/16 (mac): Revise syntax for CanonicalizeIndicesLSJT.
  + 07/25/16 (mac): Add WriteRelativeOperatorLSJT.
  + 08/16/16 (mac): Add WriteRelativeCMOperatorComponentLSJT and
    corresponding gather function.
  + 11/04/16 (mac): Remove dependency on Eigen/Core.
  + 03/26/17 (mac): Add conjugation phase relation for
      spherical-harmonic-like operators (still kHermitian).
  + 03/28/17 (mac):
    - Add constructors for OperatorLabelsJT.
    - Add ConstructZeroOperatorRelativeLSJT.
  + 04/04/17 (mac):
    - Rename CanonicalizeIndicesLSJT to CanonicalizeIndicesJT.
    - Add overload to CanonicalizeIndicesJT which only
      canonicalizes subspaces, with no reference to states.
  + 07/01/17 (mac): Remove constraints on operator labels for
      ConstructZeroOperatorRelativeLSJT.
  + 05/05/18 (mac): Remove constraint on operator labels for
      CanonicalizeIndicesJT.
  + 09/06/18 (mac): Split out JT operator generic code to
      jt_operator from lsjt_operator.

****************************************************************/

#ifndef BASIS_LSJT_OPERATOR_H_
#define BASIS_LSJT_OPERATOR_H_

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

//#include "eigen3/Eigen/Core"

#include "am/am.h"

#include "basis/lsjt_scheme.h"
#include "basis/operator.h"
#include "basis/jt_operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  //
  // relative LSJT operator
  //
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // general relative LSJT two-body operator file format
  //
  ////////////////////////////////////////////////////////////////
  //
  // Header
  //
  //   Header format:
  //
  //     # RELATIVE LSJT
  //     # ...
  //     version                                 # version
  //     J0 g0 T0_min T0_max symmetry_phase_mode # operator tensor properties
  //     Nmax Jmax                               # relative basis truncation
  //
  //   The header may start with one or more contiguous comment lines,
  //   which are designated by a hash character in the first column.
  //
  //   Then, the header contains the following fields:
  //
  //     version : file format version (1 = this version)
  //
  //     J0 (int) : angular momentum of operator
  //
  //     g0 (int) : parity grade of operator [P0=(-)^g0] (g0=0,1)
  //
  //     T0_min T0_max (int) : range of operator isospin components (a
  //       two-body operator may in general have components T0=0,1,2
  //       from the coupling of four isospin-1/2 fermionic operators)
  //
  //     symmetry_phase_mode (int) : RESERVED to describe how to
  //       obtain phase for lower triangle (see "Conjugation symmetry"
  //       below); currently only symmetry_phase_mode=0 is defined
  //
  //       Code note: An enum type basis::SymmetryPhaseMode is defined
  //       for this field, with kHermitian=0.
  //
  //     Nmax (int) : oscillator truncation of relative space
  //     (Nmax>=0)
  //
  //     Jmax (int) : additional relative angular momentum truncation
  //       of relative space (Jmax<=Nmax+1); Jmax=Nmax+1 includes the
  //       full relative space at the given oscillator truncation, but
  //       operators may often be truncated at lower partial waves
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
  //   The symmetry_phase_mode field in the header is reserved to
  //   provide information on the correct form to use for this phase
  //   factor.  However, for now, only the placeholder value
  //   0 (=kHermitian) is defined for symmetry_phase_mode, and phase
  //   conventions are only well-defined for a Hamiltonian-like (J0=0,
  //   g0=0) operator.
  //
  //   symmetry_phase_mode=0 (=kHermitian), J0=0, g0=0: For a
  //   Hamiltonian-like operator, we expect Hermiticity, i.e.,
  //   symmetry of (M_J,M_T)-branched matrix elements.  Within a
  //   diagonal sector, this means that the lower triangle is obtained
  //   from the upper triangle by ordinary symmetry.  For off-diagonal
  //   sectors, the appropriate symmetry on the JT-reduced matrix
  //   elements in general includes phase and dimension factors from
  //   the isospin Clebsches:
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
  //   For an operator which transforms under conjugation as a
  //   spherical harmonic, if we look purely at the angular momentum
  //   dependence, we in general obtain the conjugated RME (under
  //   group theory conventions for the RME)
  //
  //          <J||A_{J0}||J'> = (-)^(J'-J)*Hat(J')/Hat(J)*<J'||A_{J0}||J>
  //
  //   This relation applies to both the M1 and E2 operators under
  //   Condon-Shortley phase conventions (Suhonen Ch. 6) for these
  //   operators.  It is derived from the W-E theorem, symmetry of CG
  //   coefficient, and conjugation properties of these operators.
  //
  //   Putting this together, for a JT basis, a spherical
  //   harmonic-like operator conjugates as
  //
  //     <a,J,T,g || A_{J0,T0} || a',J,T',g>
  //       = (-)^(J'-J)*Hat(J')/Hat(J)
  //         * (-)^(T'-T)*Hat(T')/Hat(T)
  //         * <a',J,T',g || A_{J0,T0} || a,J,T,g>
  //
  //   Note that Hermitian conjugation behavior is recovered in the
  //   special case of J0=0.
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
  // operator labeling information
  ////////////////////////////////////////////////////////////////

  struct RelativeOperatorParametersLSJT
    : OperatorLabelsJT
  // Parameters for relative operator storage.
  //
  // Contains operator tensorial properties (inherited from
  // OperatorLabelsJT), plus relative basis truncation parameters.
  {

    RelativeOperatorParametersLSJT()
      // default constructor
      : Nmax(0), Jmax(0)
      {}

    RelativeOperatorParametersLSJT(
        const basis::OperatorLabelsJT& operator_labels, int Nmax_, int Jmax_)
      // Construct using given operator labels, plus given version and basis parameters.
      : OperatorLabelsJT(operator_labels), Nmax(Nmax_), Jmax(Jmax_)
    {}

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
  //   os (std::ostream): text-mode output stream
  //   parameters (basis::RelativeOperatorParametersLSJT)
  //     : operator parameters

  void ReadRelativeOperatorParametersLSJT(
      std::istream& is,
      basis::RelativeOperatorParametersLSJT& parameters
    );
  // Read file header for relative operator in LSJT scheme.
  //
  // Arguments:
  //   is (std::istream): text-mode input stream
  //   parameters (basis::RelativeOperatorParametersLSJT, output)
  //     : operator parameters

  void WriteRelativeOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const basis::RelativeSectorsLSJT& sectors,
      const basis::OperatorBlocks<double>& matrices
    );
  // Write single isospin component of a relative operator in LSJT
  // scheme.
  //
  // Side effect: The floating point precision attribute of the output
  // stream is modified.
  //
  // Arguments:
  //   os (std::ostream): text-mode output stream
  //   T0 (int): isospin for this isospin component
  //   sector (basis::RelativeSectorsLSJT):  sectors defining operator
  //   matrices (basis::OperatorBlocks<double>):  matrices defining operator

  void ReadRelativeOperatorComponentLSJT(
      std::istream& is,
      int T0,
      const basis::RelativeSectorsLSJT& sectors,
      basis::OperatorBlocks<double>& matrices
    );
  // Read single isospin component of a relative operator in LSJT
  // scheme.
  //
  // Arguments:
  //   is (std::istream): text-mode inpus stream
  //   T0 (int): isospin for this isospin component
  //   sector (basis::RelativeSectorsLSJT):  sectors defining operator
  //   matrices (basis::OperatorBlocks<double>, output):  matrices defining operator

  void ReadRelativeOperatorLSJT(
      const std::string& relative_filename,
      basis::RelativeSpaceLSJT& relative_space,
      basis::OperatorLabelsJT& operator_labels,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      bool verbose
    );
  // Set up and read relative operator.
  //
  // Arguments:
  //   parameters (Parameters): includes tensorial properties of operator
  //      choice of operator to use
  //   relative_space (..., output): target space, based on parameters in file
  //   operator_labels (basis::OperatorLabelsJT, output): operator labels, from file
  //   relative_component_sectors (..., output): target sectors
  //   relative_component_matrices (..., output): target matrices
  //   verbose (bool): whether or not to include diagnostic output

  void WriteRelativeOperatorLSJT(
      const std::string& relative_filename,
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::OperatorLabelsJT& operator_labels,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      bool verbose
    );
  // Set up and read relative operator.
  //
  // Arguments:
  //   parameters (Parameters): includes tensorial properties of operator
  //      choice of operator to use
  //   relative_space (...): target space
  //   operator_labels (basis::OperatorLabelsJT): operator labels
  //   relative_component_sectors (..., output): source sectors
  //   relative_component_matrices (..., output): source matrices
  //   verbose (bool): whether or not to include diagnostic output


  ////////////////////////////////////////////////////////////////
  // relative LSJT operator construction
  ////////////////////////////////////////////////////////////////

  void ConstructZeroOperatorRelativeLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices
    );
  // Construct zero operator in relative LSJT basis.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (input): tensorial properties of operator
  //   relative_space (input): target space
  //   relative_component_sectors (output): target sectors
  //   relative_component_matrices (output): target matrices

  void ConstructIdentityOperatorRelativeLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices
    );
  // Construct identity operator in relative LSJT basis.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (input): tensorial properties of operator
  //   relative_space (input): target space
  //   relative_component_sectors (output): target sectors
  //   relative_component_matrices (output): target matrices

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  //
  // relative-cm (two-body) LSJT operator
  //
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // general relative-cm two-body operator file format
  //
  ////////////////////////////////////////////////////////////////
  //
  // Header
  //
  //   Header format:
  //
  //     # RELATIVE-CM LSJT
  //     # ...
  //     version                                 # version
  //     J0 g0 T0_min T0_max symmetry_phase_mode # operator tensor properties
  //     Nmax                                    # relative-cm basis truncation
  //
  //   The header may start with one or more contiguous comment lines,
  //   which are designated by a hash character in the first column.
  //
  //   Then, the header contains the following fields:
  //
  //     version : file format version (1 = this version)
  //
  //     J0 (int) : angular momentum of operator
  //
  //     g0 (int) : parity grade of operator [P0=(-)^g0] (g0=0,1)
  //
  //     T0_min T0_max (int) : range of operator isospin components (a
  //       two-body operator may in general have components T0=0,1,2
  //       from the coupling of four isospin-1/2 fermionic operators)
  //
  //     symmetry_phase_mode (int) : RESERVED to describe how to
  //       obtain phase for lower triangle (see "Conjugation symmetry"
  //       below); currently only symmetry_phase_mode=0 is defined
  //
  //       Code note: An enum type basis::SymmetryPhaseMode is defined
  //       for this field, with kHermitian=0.
  //
  //     Nmax (int) : oscillator truncation of two-body space (Nmax>=0)
  //
  // Data
  //
  //   Data lines are of the form:
  //
  //     T0  Nr' lr' Nc' lc' L' S' J' T' g'  Nr lr Nc lc L S J T g  JT-RME
  //
  //   Although the g label is redundant (it can be deduced from lr and
  //   lc), it is included to make the sector structure more easily
  //   apparent to a human reader.
  //
  //   Here JT-RME is the JT-reduced matrix element under group
  //   theory conventions (i.e., no dimension factor in the
  //   Wigner-Eckart theorem):
  //
  //      < Nr' lr' Nc' lc' L' S' J' T' || op || Nr lr Nc lc L S J T >
  //
  //   Iteration follows the usual scheme within the basis module: sectors are
  //   lexicographic by (bra,ket) subspace indices, then matrix elements within
  //   a sector are lexicographic by (bra,ket) state indices.
  //
  //   It is assumed that the RelativeCMSectorsLSJT was constucted with the
  //   direction=kCanonical option.
  //
  //   Note: The concept of AS or NAS matrix elements does not apply in
  //   relative-cm states, so no normalization conversion is needed.
  //
  // For more information...
  //
  //   See the documentation of the *relative* LSJT operator file format for
  //   further general discussion: "Iteration order and symmetry", "Conjugation
  //   symmetry", and "Radial oscillator phase convention".  These discussions
  //   are generic to the LSJT scheme and apply the same to relative or
  //   relative-cm operators.

  ////////////////////////////////////////////////////////////////
  // relative-cm LSJT operator -- gather N blocks
  ////////////////////////////////////////////////////////////////

  void GatherOperatorRelativeCMLSJTNToRelativeCMLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
      const std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& relative_cm_lsjtn_component_matrices,
      const basis::RelativeCMSpaceLSJT& relative_cm_lsjt_space,
      std::array<basis::RelativeCMSectorsLSJT,3>& relative_cm_lsjt_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_cm_lsjt_component_matrices
    );
  // Assemble relative-cm representation of operator in LSJT basis,
  // from relative-cm representation in LSJTN basis, i.e., gathering
  // the matrix elements from different N blocks.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Symmetry: The lower triangle of diagonal sectors is
  // zero-initialized, but not populated.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT): tensorial properties of operator
  //   relative_cm_lsjtn_space (...): source space
  //   relative_cm_lsjtn_component_sectors (...): source sectors
  //   relative_cm_lsjtn_component_matrices (...): source matrices
  //   relative_cm_lsjt_space (...): target space
  //   relative_cm_lsjt_component_sectors (..., output): target sectors
  //   relative_cm_lsjt_component_matrices (..., output): target matrices

  ////////////////////////////////////////////////////////////////
  // relative-cm LSJT operator output
  ////////////////////////////////////////////////////////////////

  // Note that the primary intention of the output for relative-cm
  // operators is for diagnostic purposes.
  //
  // Data lines are of the form:
  //
  //   T0  Nr' lr' Nc' lc' L' S' J' T' g'  Nr lr Nc lc L S J T g  JT-RME
  //
  // Although the g label is redundant (it can be deduced from lr and
  // lc), it is included to make the sector structure more easily
  // apparent to a human reader.
  //
  // Iteration follows the usual scheme within the basis module:
  // sectors are lexicographic by (bra,ket) subspace indices, then
  // matrix elements within a sector are lexicographic by (bra,ket)
  // state indices.
  //
  // Note: The concept of AS or NAS matrix elements does not apply in
  // relative-cm states, so no normalization conversion is needed.

  void WriteRelativeCMOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const basis::RelativeCMSectorsLSJT& sectors,
      const basis::OperatorBlocks<double>& matrices
    );
  // Write single isospin component of a relative-cm operator in LSJT
  // scheme.
  //
  // Side effect: The floating point precision attribute of the output
  // stream is modified.
  //
  // Arguments:
  //   os (std::ostream): text-mode output stream
  //   T0 (int): isospin for this isospin component
  //   sectors (basis::RelativeCMSectorsLSJT): sectors defining operator
  //   matrices (basis::OperatorBlocks<double>): matrices defining operator

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  //
  // two-body (lab-frame) LSJT operator
  //
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // two-body LSJT operator -- gather N blocks
  ////////////////////////////////////////////////////////////////

  void GatherOperatorTwoBodyLSJTNToTwoBodyLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceLSJTN& two_body_lsjtn_space,
      const std::array<basis::TwoBodySectorsLSJTN,3>& two_body_lsjtn_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& two_body_lsjtn_component_matrices,
      const basis::TwoBodySpaceLSJT& two_body_lsjt_space,
      std::array<basis::TwoBodySectorsLSJT,3>& two_body_lsjt_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& two_body_lsjt_component_matrices
    );
  // Assemble two-body representation of operator in LSJT basis, from
  // two-body representation in LSJTN basis, i.e., gathering the
  // matrix elements from different N blocks.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Symmetry: The lower triangle of diagonal sectors is
  // zero-initialized, but not populated.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT): tensorial properties of operator
  //   two_body_lsjtn_space (...): source space
  //   two_body_lsjtn_component_sectors (...): source sectors
  //   two_body_lsjtn_component_matrices (...): source matrices
  //   two_body_lsjt_space (...): target space
  //   two_body_lsjt_component_sectors (..., output): target sectors
  //   two_body_lsjt_component_matrices (..., output): target matrices

  ////////////////////////////////////////////////////////////////
  // two-body LSJT operator output
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
      const basis::OperatorBlocks<double>& matrices,
      basis::NormalizationConversion conversion_mode
    );
  // Write single isospin component of a two-body operator in LSJT
  // scheme.
  //
  // Side effect: The floating point precision attribute of the output
  // stream is modified.
  //
  // Arguments:
  //   os (std::ostream): text-mode output stream
  //   T0 (int): isospin for this isospin component
  //   sectors (basis::TwoBodySectorsLSJT): sectors defining operator
  //   matrices (basis::OperatorBlocks<double>): matrices defining operator
  //   conversion (basis::NormalizationConversion): specifies any
  //     conversion between AS and NAS for output


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
