/// @file
/****************************************************************
  lsjt_operator.h

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
  7/10/16 (mac):
    - Add support for canonicalizing indices in general
      LSJT matrix element lookup (for Hamiltonian-like operators).
    - Define OperatorLabelsJT and add documentation on operators.
    - Incorporate some basic LSJT operator construction and
      manipulation functions.
  7/13/16 (mac): Revise code for LSJTN->LSJT gathering operation.
  7/20/16 (mac): Add ReadRelativeOperatorLSJT.
  7/22/16 (mac): Revise syntax for CanonicalizeIndicesLSJT.
  7/25/16 (mac): Add WriteRelativeOperatorLSJT.
  8/16/16 (mac): Add WriteRelativeCMOperatorComponentLSJT and
    corresponding gather function.
  11/4/16 (mac): Remove dependency on Eigen/Core.

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
  //     version : file format version (=1)
  //     J0 (int) : angular momentum of operator
  //     g0 (int) : parity grade of operator [P0=(-)^g0] (g0=0,1)
  //     T0_min T0_max (int) : range of operator isospin components (a two-body operator
  //       may in general have components T0=0,1,2 from the coupling of four isospin-1/2
  //       fermionic operators)
  //     symmetry_phase_mode (int) : RESERVED to describe how to obtain phase for lower triangle
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
  //   The symmetry_phase_mode field in the header is reserved to provide
  //   information on the correct form to use for this phase factor.
  //   However, for now, only the placeholder value kHermitian(=0) is
  //   defined for symmetry_phase_mode, and phase conventions are only
  //   well-defined for a Hamiltonian-like (J0=0, g0=0) operator.
  //
  //   symmetry_phase_mode=kHermitian, J0=0, g0=0: For a Hamiltonian-like operator,
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
  //
  // internal representation of an operator in JT scheme
  //
  ////////////////////////////////////////////////////////////////
  //
  // The data structures to be used for the internal representation of
  // an operator in a JT coupled basis consist of:
  //
  //   operator_labels (OperatorLabelsJT) : general information needed
  //     to work with operator (see comment on relative LSJT format
  //     above for a description of these labels):
  //
  //         J0, g0, T0_min, T0_max, symmetry_phase_mode
  //
  //     For a RelativeLSJT operator, one can use the daughter type
  //     RelativeOperatorParametersLSJT, which contains additional
  //     fields describing the Nmax and Jmax truncation.  This
  //     daughter type is primarily for internal use by the relative
  //     operator file I/O functions (and may otherwise be ignored).
  //
  //   component_sectors (std::array<tSectorsType,3>) : the sector
  //     enumerations for each isospin component (T0=0,1,2) of the
  //     operator -- Here, tSectorsType represents the indexing
  //     Sectors type for the particular basis being used, e.g.,
  //     RelativeCMSectorsLSJTN.  In general, the sector enumerations
  //     should contain "upper trianglar" sectors only (see notes on
  //     relative LSJT file format).
  //
  //   component_matrices (std::array<basis::MatrixVector,3>) : the
  //     matrix representations for each isospin component (T0=0,1,2)
  //     of the operator -- In general, for diagonal sectors, only the
  //     matrix elements in the "upper triangle" are guaranteed to be
  //     populated.  However, individual functions are free to further
  //     specify that they require input matrices or give output
  //     matrices in which the lower triangle is filled as well.
  //
  // Not all three isospin components (T0=0,1,2) need actually be
  // *used* for a given operator, depending on the parameter values
  // T0_min and T0_max.  However, we specify that the arrays have the
  // fixed size of 3 to keep the access scheme (i.e.,
  // component_sectors[T0] and component_matrices[T0]) consistent.
  // Any unused isospin components will just be default initialized
  // and can be ignored.
  //
  // A function which works with this operator *might* also need
  // access to the full space definition (e.g., RelativeCMSSpaceLSJTN)
  // to look up subspaces directly by their labels, etc., if the
  // required information is not all readily available through
  // component_sectors.

  ////////////////////////////////////////////////////////////////
  // operator labeling information
  ////////////////////////////////////////////////////////////////

  enum class SymmetryPhaseMode {kHermitian=0};

  struct OperatorLabelsJT
  // Labels giving tensorial properties for generic JT operator.
  //
  // See comments on relative LSJT file format and internal
  // representation of an operator in JT scheme for description.
  {
    int J0, g0, T0_min, T0_max;
    basis::SymmetryPhaseMode symmetry_phase_mode;
  };

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

  void ReadRelativeOperatorLSJT(
      const std::string& relative_filename,
      basis::RelativeSpaceLSJT& relative_space,
      basis::OperatorLabelsJT& operator_labels,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      bool verbose
    );
  // Set up and read relative operator.
  //
  // Arguments:
  //   parameters (Parameters) : includes tensorial properties of operator
  //      choice of operator to use
  //   relative_space (..., output) : target space, based on parameters in file
  //   operator_labels (basis::OperatorLabelsJT, output) : operator labels, from file
  //   relative_component_sectors (..., output) : target sectors
  //   relative_component_matrices (..., output) : target matrices
  //   verbose (bool) : whether or not to include diagnostic output

  void WriteRelativeOperatorLSJT(
      const std::string& relative_filename,
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::OperatorLabelsJT& operator_labels,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      const std::array<basis::MatrixVector,3>& relative_component_matrices,
      bool verbose
    );
  // Set up and read relative operator.
  //
  // Arguments:
  //   parameters (Parameters) : includes tensorial properties of operator
  //      choice of operator to use
  //   relative_space (...) : target space
  //   operator_labels (basis::OperatorLabelsJT) : operator labels
  //   relative_component_sectors (..., output) : source sectors
  //   relative_component_matrices (..., output) : source matrices
  //   verbose (bool) : whether or not to include diagnostic output

  ////////////////////////////////////////////////////////////////
  // clearing operator data
  ////////////////////////////////////////////////////////////////

  template <typename tJTSectors>
    void ClearOperatorJT(
        std::array<tJTSectors,3> component_sectors,
        std::array<basis::MatrixVector,3> component_matrices
      )
    // Delete all sector and matrix data for all isospin components of
    // an operator in **JT scheme.
    {
      for (int T0=0; T0<=2; ++T0)
        {
          component_sectors[T0] = tJTSectors();
          component_matrices[T0].resize(0);
        }
    }


  ////////////////////////////////////////////////////////////////
  // operator matrix element canonicalizaton
  ////////////////////////////////////////////////////////////////

  template <typename tLSJTSpace>
    void CanonicalizeIndicesLSJT(
        const tLSJTSpace& space,
        int& bra_subspace_index, int& ket_subspace_index,
        int& bra_state_index, int& ket_state_index,
        double& canonicalization_factor,
        int J0, int T0, int g0,
        basis::SymmetryPhaseMode symmetry_phase_mode
      )
    // Convert subspace and state indices for a matrix element to
    // canonical ("upper triangle") indices.
    //
    // This is a customized wrapper for basis::CanonicalizeIndices (see
    // operator.h), for use with RelativeLSJT operators.
    //
    // Template parameters:
    //   tLSJTSpace : type for the LSJT space from which the subspaces are
    //     drawn (subspaces must provide the T(), J(), and other needed
    //      accessor)
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
    //    J0, T0, g0 (int) : operator tensorial properties
    //    symmetry_phase_mode (basis::SymmetryPhaseMode) : operator
    //      conjugation symmetry
    {

      // canonicalize indices
      bool swapped_subspaces, swapped_states;
      basis::CanonicalizeIndices(
          bra_subspace_index, ket_subspace_index,
          swapped_subspaces,
          bra_state_index, ket_state_index,
          swapped_states
        );

      // calculate canonicalization factor
      //
      // Beware that the indices now describe the "new" bra and ket
      // *after* any swap, so one must take care in matching up bra and
      // ket labels to those in any formula describing the symmetry.

      // check that case is covered
      //
      // Phase definitions are currently only provided for
      // Hamiltonian-like operators.
      assert(
          (symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian)
          && (J0==0) && (g0==0)
        );

      // case: Hamiltonian-like operator
      //
      // Recall symmetry relation:
      //
      //     <a,J,T,g || A_{T0} || a',J,T',g>
      //       = (-)^(T'-T)*Hat(T')/Hat(T)
      //         * <a',J,T',g || A_{T0} || a,J,T,g>

      canonicalization_factor = 1.;
      if (swapped_subspaces)
        {

          // retrieve sector labels (*after* swap, i.e., canonical m.e. on RHS)
          const typename tLSJTSpace::SubspaceType& bra_subspace = space.GetSubspace(
              bra_subspace_index
            );
          const typename tLSJTSpace::SubspaceType& ket_subspace = space.GetSubspace(
              ket_subspace_index
            );
          int Tp = bra_subspace.T();
          int T = ket_subspace.T();

          canonicalization_factor *= ParitySign(Tp-T)*Hat(Tp)/Hat(T);
        }
    }


  template <typename tLSJTSpace>
    std::tuple<int,int,int,int,double> CanonicalizeIndicesLSJT(
        const tLSJTSpace& space,
        int J0, int T0, int g0,
        basis::SymmetryPhaseMode symmetry_phase_mode,
        int subspace_index_bra, int subspace_index_ket,
        int state_index_bra, int state_index_ket
      )
  // Convert subspace and state indices for a matrix element to
  // canonical ("upper triangle") indices.
  //
  // This is a customized wrapper for basis::CanonicalizeIndices (see
  // operator.h), for use with RelativeLSJT operators.
  //
  // Template parameters:
  //   tLSJTSpace : type for the LSJT space from which the subspaces are
  //     drawn (subspaces must provide the T(), J(), and other needed
  //      accessor)
  //
  // Arguments:
  //   relative_space (basis::RelativeSpaceLSJT) : space, for retrieving
  //     subspace quantum numbers to calculate canonicalization factor
  //    J0, T0, g0 (int) : operator tensorial properties
  //    symmetry_phase_mode (basis::SymmetryPhaseMode) : operator
  //      conjugation symmetry
  //   bra_subspace_index, ket_subspace_index (int) :
  //     naive sector bra and ket subspace indices, possibly to be swapped
  //   bra_state_index, ket_state_index (int) :
  //     naive bra and ket state indices, possibly to be swapped if sector
  //     is diagonal sector
  //
  // Returns:
  //   (std::tuple<int,int,int,int,double>) : canonicalized indices and swap
  //      flag as:
  //
  //        subspace_index_bra,subspace_index_ket,
  //        state_index_bra,state_index_ket,
  //        canonicalization_factor
  {

    // Note: We use the local copies of the indices (on the stack)
    // as working variables.  Slightly unnerving.

    // Note: The operator labels cannot be bundled as an
    // OperatorLabelsJT, since we need a fixed T0 value, not a range
    // of T0 values.

    // canonicalize indices
    bool swapped_subspaces;
    std::tie(
        subspace_index_bra,subspace_index_ket,
        state_index_bra,state_index_ket,
        swapped_subspaces
      )
      = basis::CanonicalizeIndices(
          subspace_index_bra, subspace_index_ket,
          state_index_bra, state_index_ket
        );

    // calculate canonicalization factor
    //
    // Beware that the indices now describe the "new" bra and ket
    // *after* any swap, so one must take care in matching up bra and
    // ket labels to those in any formula describing the symmetry.

    // check that case is covered
    //
    // Phase definitions are currently only provided for
    // Hamiltonian-like operators.
    assert(
        (symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian)
        && (J0==0) && (g0==0)
      );

    // case: Hamiltonian-like operator
    //
    // Recall symmetry relation:
    //
    //     <a,J,T,g || A_{T0} || a',J,T',g>
    //       = (-)^(T'-T)*Hat(T')/Hat(T)
    //         * <a',J,T',g || A_{T0} || a,J,T,g>

    double canonicalization_factor = 1.;
    if (swapped_subspaces)
      {

        // retrieve sector labels (*after* swap, i.e., canonical m.e. on RHS)
        const typename tLSJTSpace::SubspaceType& subspace_bra = space.GetSubspace(
            subspace_index_bra
          );
        const typename tLSJTSpace::SubspaceType& subspace_ket = space.GetSubspace(
            subspace_index_ket
          );
        int Tp = subspace_bra.T();
        int T = subspace_ket.T();

        canonicalization_factor *= ParitySign(Tp-T)*Hat(Tp)/Hat(T);
      }

    // bundle return values
    return std::tuple<int,int,int,int,double>(
        subspace_index_bra,subspace_index_ket,
        state_index_bra,state_index_ket,
        canonicalization_factor
      );

  }

  ////////////////////////////////////////////////////////////////
  // relative LSJT operator construction
  ////////////////////////////////////////////////////////////////

  void ConstructIdentityOperatorRelativeLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices
    );
  // Construct identity operator in relative LSJT basis.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   relative_space (...) : target space
  //   relative_component_sectors (..., output) : target sectors
  //   relative_component_matrices (..., output) : target matrices

  ////////////////////////////////////////////////////////////////
  // relative-cm LSJT operator -- gather N blocks
  ////////////////////////////////////////////////////////////////

  void GatherOperatorRelativeCMLSJTNToRelativeCMLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
      const std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
      const std::array<basis::MatrixVector,3>& relative_cm_lsjtn_component_matrices,
      const basis::RelativeCMSpaceLSJT& relative_cm_lsjt_space,
      std::array<basis::RelativeCMSectorsLSJT,3>& relative_cm_lsjt_component_sectors,
      std::array<basis::MatrixVector,3>& relative_cm_lsjt_component_matrices
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
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   relative_cm_lsjtn_space (...) : source space
  //   relative_cm_lsjtn_component_sectors (...) : source sectors
  //   relative_cm_lsjtn_component_matrices (...) : source matrices
  //   relative_cm_lsjt_space (...) : target space
  //   relative_cm_lsjt_component_sectors (..., output) : target sectors
  //   relative_cm_lsjt_component_matrices (..., output) : target matrices

  ////////////////////////////////////////////////////////////////
  // two-body LSJT operator -- gather N blocks
  ////////////////////////////////////////////////////////////////

  void GatherOperatorTwoBodyLSJTNToTwoBodyLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceLSJTN& two_body_lsjtn_space,
      const std::array<basis::TwoBodySectorsLSJTN,3>& two_body_lsjtn_component_sectors,
      const std::array<basis::MatrixVector,3>& two_body_lsjtn_component_matrices,
      const basis::TwoBodySpaceLSJT& two_body_lsjt_space,
      std::array<basis::TwoBodySectorsLSJT,3>& two_body_lsjt_component_sectors,
      std::array<basis::MatrixVector,3>& two_body_lsjt_component_matrices
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
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   two_body_lsjtn_space (...) : source space
  //   two_body_lsjtn_component_sectors (...) : source sectors
  //   two_body_lsjtn_component_matrices (...) : source matrices
  //   two_body_lsjt_space (...) : target space
  //   two_body_lsjt_component_sectors (..., output) : target sectors
  //   two_body_lsjt_component_matrices (..., output) : target matrices

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
  // Although the g label is redundant (it can be deduced from l1 and
  // l2), it is included to make the sector structure more easily
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
      const basis::MatrixVector& matrices
    );
  // Write single isospin component of a relative-cm operator in LSJT
  // scheme.
  //
  // Side effect: The floating point precision attribute of the output
  // stream is modified.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   T0 (int) : isospin for this isospin component
  //   sectors (basis::RelativeCMSectorsLSJT) : sectors defining operator
  //   matrices (basis::MatrixVector) : matrices defining operator

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
  //   sectors (basis::TwoBodySectorsLSJT) : sectors defining operator
  //   matrices (basis::MatrixVector) : matrices defining operator
  //   conversion (basis::NormalizationConversion) : specifies any
  //     conversion between AS and NAS for output


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
