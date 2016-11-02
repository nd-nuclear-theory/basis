/// @file
/****************************************************************
  nlj_operator.h

  Defines functions for manipulation of operators for nlj orbitals.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 11/1/16 (mac): Created.

****************************************************************/

#ifndef NLJ_OPERATOR_H_
#define NLJ_OPERATOR_H_

#include "basis/nlj_orbital.h"
#include "basis/operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // matrix element lookup
  ////////////////////////////////////////////////////////////////

  double MatrixElementLJPN(
      const basis::OrbitalSpaceLJPN& bra_orbital_space,
      const basis::OrbitalSpaceLJPN& ket_orbital_space,
      const basis::OrbitalSectorsLJPN& sectors,
      const basis::MatrixVector& matrices,
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
  //   matrices (basis::MatrixVector): matrices corresponding to these sectors
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
