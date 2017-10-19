/****************************************************************
  jjjpn_operator.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "jjjpn_operator.h"

#include <iomanip>
#include <iostream>

namespace basis {

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
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
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
        for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
          for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
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
                  if ((bra.index1()==bra.index2())&&(bra.two_body_species()!=basis::TwoBodySpeciesPN::kPN))
                    conversion_factor *= (1/sqrt(2.));
                  if ((ket.index1()==ket.index2())&&(ket.two_body_species()!=basis::TwoBodySpeciesPN::kPN))
                    conversion_factor *= (1/sqrt(2.));
                }
              else if (conversion_mode == basis::NormalizationConversion::kNASToAS)
                {
                  if ((bra.index1()==bra.index2())&&(bra.two_body_species()!=basis::TwoBodySpeciesPN::kPN))
                    conversion_factor *= (sqrt(2.));
                  if ((ket.index1()==ket.index2())&&(ket.two_body_species()!=basis::TwoBodySpeciesPN::kPN))
                    conversion_factor *= (sqrt(2.));
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
