/****************************************************************
  nlj_orbital.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <cstddef>
#include "nlj_operator.h"

namespace basis {


  std::tuple<std::size_t,std::size_t,std::size_t>
  MatrixElementIndicesLJPN(
      const basis::OrbitalSpaceLJPN& bra_orbital_space,
      const basis::OrbitalSpaceLJPN& ket_orbital_space,
      const basis::OrbitalSectorsLJPN& sectors,
      const basis::FullOrbitalLabels& bra_labels,
      const basis::FullOrbitalLabels& ket_labels
    )
  {
    // extract state labels
    basis::OrbitalSpeciesPN bra_orbital_species, ket_orbital_species;
    int bra_n, ket_n, bra_l, ket_l;
    HalfInt bra_j, ket_j;
    std::tie(bra_orbital_species,bra_n,bra_l,bra_j) = bra_labels;
    std::tie(ket_orbital_species,ket_n,ket_l,ket_j) = ket_labels;

    // look up LJPN sector
    std::size_t bra_subspace_index = bra_orbital_space.LookUpSubspaceIndex(
        typename basis::OrbitalSubspaceLJPN::SubspaceLabelsType(bra_orbital_species,bra_l,bra_j)
      );
    std::size_t ket_subspace_index = ket_orbital_space.LookUpSubspaceIndex(
        typename basis::OrbitalSubspaceLJPN::SubspaceLabelsType(ket_orbital_species,ket_l,ket_j)
      );
    std::size_t sector_index = sectors.LookUpSectorIndex(bra_subspace_index,ket_subspace_index);

    // short-circuit on missing sector
    if (sector_index==basis::kNone)
      {
        std::size_t bra_state_index = basis::kNone;
        std::size_t ket_state_index = basis::kNone;
        return std::make_tuple(sector_index,bra_state_index,ket_state_index);
      }

    // set up sector alias
    const typename basis::OrbitalSectorsLJPN::SectorType& sector = sectors.GetSector(sector_index);

    // retrieve state indices

    std::size_t bra_state_index, ket_state_index;
    assert((bra_n>=0)&&(ket_n>=0));  // pjf says: don't be dumb
    bra_state_index
      = sector.bra_subspace().LookUpStateIndex(basis::OrbitalStateLJPNLabels(bra_n));
    ket_state_index
      = sector.ket_subspace().LookUpStateIndex(basis::OrbitalStateLJPNLabels(ket_n));

    return std::make_tuple(sector_index,bra_state_index,ket_state_index);

  }

  double MatrixElementLJPN(
      const basis::OrbitalSpaceLJPN& bra_orbital_space,
      const basis::OrbitalSpaceLJPN& ket_orbital_space,
      const basis::OrbitalSectorsLJPN& sectors,
      const basis::OperatorBlocks<double>& matrices,
      const basis::OrbitalStatePN& bra, const basis::OrbitalStatePN& ket
    )
  {
    // look up LJPN sector
    std::size_t sector_index, bra_state_index, ket_state_index;
    std::tie(sector_index,bra_state_index,ket_state_index)
      = MatrixElementIndicesLJPN(
          bra_orbital_space,ket_orbital_space,sectors,
          bra.full_labels(),ket.full_labels()
        );

    // validate returned indices
    //   with helpful messages for tracking down indexing lookup failures

    // trap failed sector lookup
    if (sector_index==basis::kNone)
      {
        std::cerr
          << "ERROR: basis::MatrixElementLJPN: missing sector while looking up radial matrix element "
          << bra.LabelStr() << " " << ket.LabelStr()
          << std::endl;
        std::exit(EXIT_FAILURE);
      }

    // set up sector alias
    const typename basis::OrbitalSectorsLJPN::SectorType& sector = sectors.GetSector(sector_index);

    // trapping subspace overruns
    if ((bra_state_index==basis::kNone)||(ket_state_index==basis::kNone))
      {
        std::cerr
          << "ERROR: basis::MatrixElementLJPN: radial quantum number not found in orbital subspace "
          << bra.LabelStr() << " " << ket.LabelStr()
          << std::endl;
      }

    double matrix_element = matrices[sector_index](bra_state_index,ket_state_index);

    return matrix_element;

  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
