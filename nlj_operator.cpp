/****************************************************************
  nlj_orbital.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "nlj_operator.h"

namespace basis {


  double MatrixElementLJPN(
      const basis::OrbitalSpaceLJPN& bra_orbital_space,
      const basis::OrbitalSpaceLJPN& ket_orbital_space,
      const basis::OrbitalSectorsLJPN& sectors,
      const basis::MatrixVector& matrices,
      const basis::OrbitalStatePN& bra, const basis::OrbitalStatePN& ket
    )
  {
    // extract state labels
    basis::OrbitalSpeciesPN bra_orbital_species = bra.orbital_species();
    basis::OrbitalSpeciesPN ket_orbital_species = ket.orbital_species();
    int bra_n, ket_n, bra_l, ket_l;
    HalfInt bra_j, ket_j;
    std::tie(bra_n,bra_l,bra_j) = bra.labels();
    std::tie(ket_n,ket_l,ket_j) = ket.labels();
    
    // look up LJPN sector
    int bra_subspace_index = bra_orbital_space.LookUpSubspaceIndex(
        typename basis::OrbitalSubspaceLJPN::SubspaceLabelsType(bra_orbital_species,bra_l,bra_j)
      );
    int ket_subspace_index = ket_orbital_space.LookUpSubspaceIndex(
        typename basis::OrbitalSubspaceLJPN::SubspaceLabelsType(ket_orbital_species,ket_l,ket_j)
      );
    int sector_index = sectors.LookUpSectorIndex(bra_subspace_index,ket_subspace_index);
    assert(sector_index!=basis::kNone);  // trap failed lookup
    const typename basis::OrbitalSectorsLJPN::SectorType& sector = sectors.GetSector(sector_index);

    // retrieve LJPN matrix element
    //
    // We rely on the assumption that the label n is equivalent to the
    // index in an LJPN subspace.
    assert(bra_n<sector.bra_subspace().size());
    assert(ket_n<sector.ket_subspace().size());
    double matrix_element = matrices[sector_index](bra_n,ket_n);

    return matrix_element;

  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
