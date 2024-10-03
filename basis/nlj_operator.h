/************************************************************//**
  @file nlj_operator.h

  Defines functions for manipulation of operators for nlj orbitals.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 11/1/16 (mac): Created.
  + 11/4/16 (mac): Add informative error messages on failed MatrixElementLJPN lookup.
  + 11/21/16 (mac): Break out matrix element indexing lookup from value lookup.
  + 09/22/17 (pjf): Look up state indices honestly in MatrixElementIndicesLJPN
  + 02/20/17 (pjf): Add OneBodyOperatorLJPN.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.
****************************************************************/

#ifndef BASIS_NLJ_OPERATOR_H_
#define BASIS_NLJ_OPERATOR_H_

#include <cstddef>
#include <string>
#include <tuple>

#include "nlj_orbital.h"
#include "operator.h"

namespace basis {

  enum class OneBodyOperatorType : char {
    kRadial = 'R',
    kSpherical = 'S'
  };

  ////////////////////////////////////////////////////////////////
  // matrix element lookup
  ////////////////////////////////////////////////////////////////

  std::tuple<std::size_t,std::size_t,std::size_t>
  MatrixElementIndicesLJPN(
      const OrbitalSpaceLJPN& bra_orbital_space,
      const OrbitalSpaceLJPN& ket_orbital_space,
      const OrbitalSectorsLJPN& sectors,
      const FullOrbitalLabels& bra_labels,
      const FullOrbitalLabels& ket_labels
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
  //   (std::tuple<std::size_t,std::size_t,std::size_t>): the matrix element indices
  //     (sector_index,bra_state_index,ket_state_index)

  double MatrixElementLJPN(
      const OrbitalSpaceLJPN& bra_orbital_space,
      const OrbitalSpaceLJPN& ket_orbital_space,
      const OrbitalSectorsLJPN& sectors,
      const OperatorBlocks<double>& matrices,
      const OrbitalStatePN& bra, const OrbitalStatePN& ket
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
  // one-body operator indexing and storage
  ////////////////////////////////////////////////////////////////

  struct OneBodyOperatorLJPN
  {
    OneBodyOperatorType operator_type;
    OrbitalSpaceLJPN bra_orbital_space;
    OrbitalSpaceLJPN ket_orbital_space;
    OrbitalSectorsLJPN sectors;
    OperatorBlocks<double> matrices;
    std::string name;

    double get_matrix_element(const OrbitalStatePN& bra, const OrbitalStatePN& ket) const
    {
      return MatrixElementLJPN(bra_orbital_space, ket_orbital_space, sectors, matrices, bra, ket);
    }
  };

  ////////////////////////////////////////////////////////////////
  /// @}
  ////////////////////////////////////////////////////////////////
}  // namespace basis

#endif  // BASIS_NLJ_OPERATOR_H_
