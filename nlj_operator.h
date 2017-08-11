/************************************************************//**
  @file nlj_operator.h

  Defines functions for manipulation of operators for nlj orbitals.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 11/1/16 (mac): Created.
  + 11/4/16 (mac): Add informative error messages on failed
    MatrixElementLJPN lookup.
  + 11/21/16 (mac): Break out matrix element indexing lookup
    from value lookup.

****************************************************************/

#ifndef BASIS_NLJ_OPERATOR_H_
#define BASIS_NLJ_OPERATOR_H_

#include "basis/nlj_orbital.h"
#include "basis/operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // matrix element lookup
  ////////////////////////////////////////////////////////////////

  std::tuple<int,int,int>
  MatrixElementIndicesLJPN(
      const basis::OrbitalSpaceLJPN& bra_orbital_space,
      const basis::OrbitalSpaceLJPN& ket_orbital_space,
      const basis::OrbitalSectorsLJPN& sectors,
      const basis::FullOrbitalLabels& bra_labels,
      const basis::FullOrbitalLabels& ket_labels
    );
  // Look up radial matrix element indices by (species,n,l,j) labels.
  //
  // Returns basis::kNone for failed index lookups.
  //
  // Arguments:
  //   bra_orbital_space (basis::OrbitalSpaceLJPN): LJPN bra space
  //   ket_orbital_space (basis::OrbitalSpaceLJPN): LJPN ket space
  //   sectors (basis::OrbitalSectorsLJPN): LJPN sectors for operator
  //   matrices (basis::OperatorBlocks<double>): matrices corresponding to these sectors
  //   bra_labels (basis::FullOrbitalLabels): bra state labels
  //   ket_labels (basis::FullOrbitalLabels): ket state labels
  //
  // Returns:
  //   (std::tuple<int,int,int>): the matrix element indices
  //     (sector_index,bra_state_index,ket_state_index)

  double MatrixElementLJPN(
      const basis::OrbitalSpaceLJPN& bra_orbital_space,
      const basis::OrbitalSpaceLJPN& ket_orbital_space,
      const basis::OrbitalSectorsLJPN& sectors,
      const basis::OperatorBlocks<double>& matrices,
      const basis::OrbitalStatePN& bra, const basis::OrbitalStatePN& ket
    );
  // Look up radial matrix element by (species,n,l,j) labels.
  //
  // Fails with assertion if matrix element is not found, i.e., not
  // supported by the radial operator's indexing.
  //
  // Arguments:
  //   bra_orbital_space (basis::OrbitalSpaceLJPN): LJPN bra space
  //   ket_orbital_space (basis::OrbitalSpaceLJPN): LJPN ket space
  //   sectors (basis::OrbitalSectorsLJPN): LJPN sectors for operator
  //   matrices (basis::OperatorBlocks<double>): matrices corresponding to these sectors
  //   bra (basis::OrbitalStatePN): bra state (in PN space)
  //   ket (basis::OrbitalStatePN): ket state (in PN space)
  //
  // Returns:
  //   (double): the matrix element

  ////////////////////////////////////////////////////////////////
  /// @}
  ////////////////////////////////////////////////////////////////
}  // namespace basis

#endif  // NLJ_ORBITAL_H_
