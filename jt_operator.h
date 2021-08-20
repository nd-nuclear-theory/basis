/****************************************************************
  @file jt_operator.h

  Defines generic JT operator labels and templatized functions for manipulation
  of JT operators.  These are generic to LSJT and JJTT schemes.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 11/21/15 (mac): Created (lsjt_interaction.h), from code in
    moshinsky_xform_lsjt.
  + 06/08/16 (mac): Update for basis package and update conventions.
    - Rename to lsjt_operator.h.
    - Remove matrix-style output.
    - Remove generic template functions.
  + 07/10/16 (mac):
    - Define OperatorLabelsJT and add documentation on operators.
  + 11/04/16 (mac): Remove dependency on Eigen/Core.
  + 03/26/17 (mac): Add conjugation phase relation for
      spherical-harmonic-like operators (still kHermitian).
  + 03/28/17 (mac):
    - Add constructors for OperatorLabelsJT.
  + 04/04/17 (mac):
    - Rename CanonicalizeIndicesLSJT to CanonicalizeIndicesJT.
    - Add overload to CanonicalizeIndicesJT which only
      canonicalizes subspaces, with no reference to states.
  + 05/05/18 (mac): Remove constraint on operator labels for
      CanonicalizeIndicesJT.
  + 09/06/18 (mac): Create jt_operator, by splitting out JT
      operator generic code from lsjt_operator.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.

****************************************************************/

#ifndef BASIS_JT_OPERATOR_H_
#define BASIS_JT_OPERATOR_H_

#include <cstddef>
#include <cassert>
#include <cstdlib>
#include <array>
#include <tuple>

#include "operator.h"

#include "am/am.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  //
  // generic constructs for operators in a JT scheme
  //
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // internal representation of an operator in JT scheme
  //
  ////////////////////////////////////////////////////////////////
  //
  // The data structures to be used for the internal representation of
  // an operator in a JT coupled basis consist of:
  //
  //   operator_labels (OperatorLabelsJT): general information needed
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
  //   component_sectors (std::array<tSectorsType,3>): the sector
  //     enumerations for each isospin component (T0=0,1,2) of the
  //     operator -- Here, tSectorsType represents the indexing
  //     Sectors type for the particular basis being used, e.g.,
  //     RelativeCMSectorsLSJTN.  In general, the sector enumerations
  //     should contain "upper trianglar" sectors only (see notes on
  //     relative LSJT file format).
  //
  //   component_matrices (std::array<basis::OperatorBlocks<double>,3>): the
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

  OperatorLabelsJT()
  // default constructor
  : J0(0), g0(0), T0_min(0), T0_max(0), symmetry_phase_mode(basis::SymmetryPhaseMode::kHermitian)
    {}

  OperatorLabelsJT(int J0_, int g0_, int T0_min_, int T0_max_, basis::SymmetryPhaseMode symmetry_phase_mode_)
  // explicit constructor
  : J0(J0_), g0(g0_), T0_min(T0_min_), T0_max(T0_max_), symmetry_phase_mode(symmetry_phase_mode_)
    {}

    int J0, g0, T0_min, T0_max;
    basis::SymmetryPhaseMode symmetry_phase_mode;
  };

  ////////////////////////////////////////////////////////////////
  // clearing operator data
  ////////////////////////////////////////////////////////////////

  template <typename tJTSectors>
    void ClearOperatorJT(
        std::array<tJTSectors,3> component_sectors,
        std::array<basis::OperatorBlocks<double>,3> component_matrices
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

  template <typename tJTSpace>
    std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,bool,double> CanonicalizeIndicesJT(
        const tJTSpace& space,
        int J0, int T0, int g0,
        basis::SymmetryPhaseMode symmetry_phase_mode,
        std::size_t subspace_index_bra, std::size_t subspace_index_ket,
        std::size_t state_index_bra, std::size_t state_index_ket
      )
  // Convert subspace and state indices for a matrix element to
  // canonical ("upper triangle") indices.
  //
  // This is a customized wrapper for basis::CanonicalizeIndices (see
  // operator.h), for use with RelativeJT operators.
  //
  // Template parameters:
  //   tJTSpace: type for the JT space from which the subspaces are
  //     drawn (subspaces must provide the T(), J(), and other needed
  //      accessor)
  //
  // Arguments:
  //   relative_space (basis::RelativeSpaceJT): space, for retrieving
  //     subspace quantum numbers to calculate canonicalization factor
  //   J0, T0, g0 (int): operator tensorial properties
  //   symmetry_phase_mode (basis::SymmetryPhaseMode): operator
  //      conjugation symmetry
  //   bra_subspace_index, ket_subspace_index (std::size_t):
  //     naive sector bra and ket subspace indices, possibly to be swapped
  //   bra_state_index, ket_state_index (std::size_t):
  //     naive bra and ket state indices, possibly to be swapped if sector
  //     is diagonal sector
  //
  // Returns:
  //   canonicalized indices and swap flag and phase as:
  //
  //        subspace_index_bra,subspace_index_ket,
  //        state_index_bra,state_index_ket,
  //        swapped_subspaces,
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
    // Hamiltonian-like or more generally spherical-harmonic-like
    // (kHermitian) operators.
    assert(symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian);

    double canonicalization_factor = 1.;
    if (swapped_subspaces)
      {

        // retrieve sector labels (*after* swap, i.e., canonical m.e. on RHS)
        const typename tJTSpace::SubspaceType& subspace_bra = space.GetSubspace(
            subspace_index_bra
          );
        const typename tJTSpace::SubspaceType& subspace_ket = space.GetSubspace(
            subspace_index_ket
          );
        int Tp = subspace_bra.T();
        int T = subspace_ket.T();
        int Jp = subspace_bra.J();
        int J = subspace_ket.J();

        canonicalization_factor *= ParitySign(Tp-T)*Hat(Tp)/Hat(T);
        canonicalization_factor *= ParitySign(Jp-J)*Hat(Jp)/Hat(J);
      }

    // bundle return values
    return std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,bool,double>(
        subspace_index_bra,subspace_index_ket,
        state_index_bra,state_index_ket,
        swapped_subspaces,
        canonicalization_factor
      );
  }

  template <typename tJTSpace>
    std::tuple<std::size_t,std::size_t,bool,double> CanonicalizeIndicesJT(
        const tJTSpace& space,
        int J0, int T0, int g0,
        basis::SymmetryPhaseMode symmetry_phase_mode,
        std::size_t subspace_index_bra, std::size_t subspace_index_ket
      )
  // Convert subspace indices for a matrix element to
  // canonical ("upper triangle") indices.
  //
  // This is an overloaded wrapper to the full version (which takes
  // subspace and state indices) and discards the state indices.
  //
  // Template parameters:
  //   tJTSpace: type for the JT space from which the subspaces are
  //     drawn (subspaces must provide the T(), J(), and other needed
  //      accessor)
  //
  // Arguments:
  //   relative_space (basis::RelativeSpaceJT): space, for retrieving
  //     subspace quantum numbers to calculate canonicalization factor
  //   J0, T0, g0 (int): operator tensorial properties
  //   symmetry_phase_mode (basis::SymmetryPhaseMode): operator
  //      conjugation symmetry
  //   bra_subspace_index, ket_subspace_index (std::size_t):
  //     naive sector bra and ket subspace indices, possibly to be swapped
  //
  // Returns:
  //   canonicalized indices and swap flag and phase as:
  //
  //        subspace_index_bra,subspace_index_ket,
  //        swapped_subspaces,
  //        canonicalization_factor
  {

    // invoke canonicalization with dummy state indices
    bool swapped_subspaces;
    double canonicalization_factor;
    std::tie(
        subspace_index_bra,subspace_index_ket,
        std::ignore,std::ignore,
        swapped_subspaces,
        canonicalization_factor
      )
      = CanonicalizeIndicesJT(
          space,
          J0,T0,g0,symmetry_phase_mode,
          subspace_index_bra,subspace_index_ket,
          0,0
        );

    // bundle return values
    return std::tuple<std::size_t,std::size_t,bool,double>(
        subspace_index_bra,subspace_index_ket,
        swapped_subspaces,
        canonicalization_factor
      );

  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
