/****************************************************************
  jjjpn_operator.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "jjjpn_operator.h"

#include <cmath>
#include <cassert>
#include <cstddef>
#include <array>
#include <iomanip>
#include <ostream>

#include "basis.h"
#include "jjjpn_scheme.h"
#include "many_body.h"
#include "operator.h"
#include "proton_neutron.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body jjJpn state lookup (with canonicalization) 
  ////////////////////////////////////////////////////////////////

  std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,double> CanonicalizeTwoBodyOrbitalIndicesJJJPN(
      const OrbitalSpacePN& space,
      int J,
      std::size_t subspace_index1, std::size_t subspace_index2,
      std::size_t state_index1, std::size_t state_index2
    )
  {

    // Note: We use the local copies of the indices (on the stack)
    // as working variables.  Slightly unnerving.

    // swap if necessary
    bool swapped = !(
        std::make_tuple(subspace_index1, state_index1)
        <= std::make_tuple(subspace_index2, state_index2)
      );
    if (swapped)
      {
        std::swap(subspace_index1, subspace_index2);
        std::swap(state_index1, state_index2);  // so index stays with sector
      }

    // calculate canonicalization factor
    double canonicalization_factor = +1.;
    if (swapped)
      {
        const HalfInt j1 = space.GetSubspace(subspace_index1).GetState(state_index1).j();
        const HalfInt j2 = space.GetSubspace(subspace_index2).GetState(state_index2).j();
        canonicalization_factor = ParitySign(J-j1-j2);
      }

    // bundle return values
    return std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,double>(
        subspace_index1, subspace_index2,
        state_index1, state_index2,
        canonicalization_factor
      );
  }

  ////////////////////////////////////////////////////////////////
  // two-body jjJpn RME lookup (with canonicalization) 
  ////////////////////////////////////////////////////////////////
  
  std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,double> CanonicalizeIndicesJJJPN(
      const TwoBodySpaceJJJPN& space,
      int J0, int g0,
      std::size_t subspace_index_bra, std::size_t subspace_index_ket,
      std::size_t state_index_bra, std::size_t state_index_ket
    )
  {

    // Note: We use the local copies of the indices (on the stack)
    // as working variables.  Slightly unnerving.

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

    double canonicalization_factor = 1.;
    if (swapped_subspaces)
      {

        // retrieve sector labels (*after* swap, i.e., canonical m.e. on RHS)
        const TwoBodySubspaceJJJPN& subspace_bra = space.GetSubspace(
            subspace_index_bra
          );
        const TwoBodySubspaceJJJPN& subspace_ket = space.GetSubspace(
            subspace_index_ket
          );
        int J_bra = subspace_bra.J();
        int J_ket = subspace_ket.J();

        canonicalization_factor *= ParitySign(J_bra-J_ket)*Hat(J_bra)/Hat(J_ket);
      }

    // bundle return values
    return std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,double>(
        subspace_index_bra, subspace_index_ket,
        state_index_bra, state_index_ket,
        canonicalization_factor
      );
  }

  ////////////////////////////////////////////////////////////////
  // two-body jjJpn operator output
  ////////////////////////////////////////////////////////////////

  void WriteTwoBodyOperatorJJJPN(
      std::ostream& os,
      const basis::TwoBodySectorsJJJPN& sectors,
      const basis::OperatorBlocks<double>& matrices,
      basis::NormalizationConversion conversion_mode,
      int indexing_base
    )
  {

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
        const typename TwoBodySectorsJJJPN::SectorType& sector = sectors.GetSector(sector_index);
        const typename TwoBodySectorsJJJPN::SubspaceType& bra_subspace = sector.bra_subspace();
        const typename TwoBodySectorsJJJPN::SubspaceType& ket_subspace = sector.ket_subspace();

        // verify that sector is canonical
        //
        // This is a check that the caller's sector construction
        // followed the specification that only "upper triangle"
        // sectors are stored.
        assert(sector.IsUpperTriangle());

        // iterate over matrix elements
        for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
          for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
            {

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // retrieve states
              const basis::TwoBodyStateJJJPN bra(bra_subspace,bra_index);
              const basis::TwoBodyStateJJJPN ket(ket_subspace,ket_index);

              // determine matrix element normalization factor
              double conversion_factor = 1.;
              if (conversion_mode == basis::NormalizationConversion::kASToNAS)
                {
                  if (bra.two_body_species()!=basis::TwoBodySpeciesPN::kPN)
                    if (bra.index1()==bra.index2())
                      conversion_factor *= (1/std::sqrt(2.));
                  if (ket.two_body_species()!=basis::TwoBodySpeciesPN::kPN)
                    if (ket.index1()==ket.index2())
                      conversion_factor *= (1/std::sqrt(2.));
                }
              else if (conversion_mode == basis::NormalizationConversion::kNASToAS)
                {
                  if (bra.two_body_species()!=basis::TwoBodySpeciesPN::kPN)
                    if (bra.index1()==bra.index2())
                      conversion_factor *= (std::sqrt(2.));
                  if (ket.two_body_species()!=basis::TwoBodySpeciesPN::kPN)
                    if (ket.index1()==ket.index2())
                      conversion_factor *= (std::sqrt(2.));
                }

              // extract matrix element
              const double matrix_element = conversion_factor*matrices[sector_index](bra_index,ket_index);

              // generate output line
              const int width = 3;
              const int precision = 8;
              os << std::setprecision(precision);
              os
                << " " << std::setw(width) << bra.index1() + indexing_base
                << " " << std::setw(width) << bra.index2() + indexing_base
                << " " << std::setw(width) << bra.J()
                << " " << std::setw(width) << bra.g()
                << " " << std::setw(width) << kTwoBodySpeciesPNCodeTz[int(bra.two_body_species())]
                << " " << "    "
                << " " << std::setw(width) << ket.index1() + indexing_base
                << " " << std::setw(width) << ket.index2() + indexing_base
                << " " << std::setw(width) << ket.J()
                << " " << std::setw(width) << ket.g()
                << " " << std::setw(width) << kTwoBodySpeciesPNCodeTz[int(ket.two_body_species())]
                << " " << "    "
                << " " << std::showpoint << std::scientific << std::setw(width+precision+5)
                << matrix_element
                << std::endl;

            }

      }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
