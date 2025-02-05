/****************************************************************

  jjjt_operator.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <ostream>
#include <tuple>
#include <vector>

#include "basis/jjjttz_operator.h"

#include "am/halfint.h"
#include <Eigen/Core>
#include "basis/basis.h"
#include "basis/many_body.h"
#include "basis/jjjttz_scheme.h"
#include "basis/jt_operator.h"
#include "basis/operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body JJJT operator -- gather N blocks
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // two-body JJJTTz operator output
  ////////////////////////////////////////////////////////////////

  void WriteTwoBodyOperatorJJJTTz(
      std::ostream& os,
      const TwoBodySectorsJJJTTz& sectors,
      const OperatorBlocks<double>& matrices,
      basis::NormalizationConversion conversion_mode
    )
  {

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
        const typename TwoBodySectorsJJJTTz::SectorType& sector = sectors.GetSector(sector_index);
        const typename TwoBodySectorsJJJTTz::SubspaceType& bra_subspace = sector.bra_subspace();
        const typename TwoBodySectorsJJJTTz::SubspaceType& ket_subspace = sector.ket_subspace();

        // verify that sector is canonical
        //
        // This is a check that the caller's sector construction
        // followed the specification that only "upper triangle"
        // sectors are stored.
        assert(sector.bra_subspace_index()<=sector.ket_subspace_index());

        // iterate over matrix elements
        for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
          for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
            {

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // define states
              const basis::TwoBodyStateJJJTTz bra(bra_subspace,bra_index);
              const basis::TwoBodyStateJJJTTz ket(ket_subspace,ket_index);

              // determine matrix element normalization factor
              double conversion_factor = 1.;
              if (conversion_mode == basis::NormalizationConversion::kASToNAS)
                {
                  if ((bra.N1()==bra.N2())&&(bra.l1()==bra.l2())&&(bra.j1()==bra.j2()))
                    conversion_factor *= (1/sqrt(2.));
                  if ((ket.N1()==ket.N2())&&(ket.l1()==ket.l2())&&(ket.j1()==ket.j2()))
                    conversion_factor *= (1/sqrt(2.));
                }
              else if (conversion_mode == basis::NormalizationConversion::kNASToAS)
                {
                  if ((bra.N1()==bra.N2())&&(bra.l1()==bra.l2())&&(bra.j1()==bra.j2()))
                    conversion_factor *= sqrt(2.);
                  if ((ket.N1()==ket.N2())&&(ket.l1()==ket.l2())&&(ket.j1()==ket.j2()))
                    conversion_factor *= sqrt(2.);
                }

              // extract matrix element
              const double matrix_element = conversion_factor*matrices[sector_index](bra_index,ket_index);

              // generate output line
              const int width = 3;
              const int precision = 8;  // for approximately single precision output
              os << std::setprecision(precision);
              os
                << " " << std::setw(width) << bra.N1()
                << " " << std::setw(width) << bra.l1()
                << " " << std::showpoint << std::fixed << std::setprecision(1) << std::setw(4) << float(bra.j1())
                << " " << std::setw(width) << bra.N2()
                << " " << std::setw(width) << bra.l2()
                << " " << std::showpoint << std::fixed << std::setprecision(1) << std::setw(4) << float(bra.j2())
                << " " << std::setw(width) << bra.J()
                << " " << std::setw(width) << bra.T()
                << " " << std::setw(width) << bra.Tz()
                << " " << std::setw(width) << bra.g()
                << " " << "    "
                << " " << std::setw(width) << ket.N1()
                << " " << std::setw(width) << ket.l1()
                << " " << std::showpoint << std::fixed << std::setprecision(1) << std::setw(4) << float(ket.j1())
                << " " << std::setw(width) << ket.N2()
                << " " << std::setw(width) << ket.l2()
                << " " << std::showpoint << std::fixed << std::setprecision(1) << std::setw(4) << float(ket.j2())
                << " " << std::setw(width) << ket.J()
                << " " << std::setw(width) << ket.T()
                << " " << std::setw(width) << bra.Tz()
                << " " << std::setw(width) << ket.g()
                << " " << "    "
                << " " << std::showpoint << std::scientific << std::setprecision(precision) << matrix_element
                << std::endl;

            }

      }
  }

  ////////////////////////////////////////////////////////////////
  // two-body JJJTTz operator: find index from labels and save matrix element
  ////////////////////////////////////////////////////////////////

  void SetTwoBodyOperatorMatrixElementJJJTTz( // requires inputs to be canonical
      const basis::TwoBodySpaceJJJTTz& space,
      const TwoBodySubspaceJJJTTzLabels& subspace_labels_bra,
      const TwoBodySubspaceJJJTTzLabels& subspace_labels_ket,
      const TwoBodyStateJJJTTzLabels& state_labels_bra,
      const TwoBodyStateJJJTTzLabels& state_labels_ket,
      const TwoBodySectorsJJJTTz& sectors,
      OperatorBlocks<double>& matrices,
      double matrix_element
    )
  {
    // bool swap_braket=false;
    std::size_t bra_subspace_index=space.LookUpSubspaceIndex(subspace_labels_bra);
    std::size_t ket_subspace_index=space.LookUpSubspaceIndex(subspace_labels_ket);
    if ((bra_subspace_index==kNone)||(ket_subspace_index==kNone))
      // gracefully handle lookup of forbidden subspace
      //
      // E.g., in ob-2 truncation, J=5, T=1, sector does not exist, since
      // (d5/2)^2 cannot couple to J+T even, but me2j I/O loop naively includes
      // this matrix element.
      if (matrix_element!=0) {
        std::cerr << "ERROR: Attempting to set nonzero value for matrix element in nonexistent sector" << std::endl;
        exit(EXIT_FAILURE);
        // return;
      } else {
        // std::cout << "0" << std::endl;
        return;
      }
    std::size_t sector_index=sectors.LookUpSectorIndex(bra_subspace_index,ket_subspace_index);
    // if (sector_index==kNone) {
    //   sector_index=sectors.LookUpSectorIndex(ket_subspace_index,bra_subspace_index);
    //   swap_braket=true;
    // }
    basis::TwoBodySectorsJJJTTz::SectorType sector=sectors.GetSector(sector_index);
    
    // look up matrix element
    const basis::TwoBodySubspaceJJJTTz& bra_subspace=sector.bra_subspace();
    const basis::TwoBodySubspaceJJJTTz& ket_subspace=sector.ket_subspace();
    std::size_t bra_state_index=bra_subspace.LookUpStateIndex(state_labels_bra);
    std::size_t ket_state_index=ket_subspace.LookUpStateIndex(state_labels_ket);
    // assert((bra_state_index!=kNone)&&(ket_state_index!=kNone));
    if ((bra_state_index==kNone)||(ket_state_index==kNone)) {
      // gracefully handle lookup of forbidden state
      if (matrix_element!=0) {
        // std::cout << "matrix_element!=0" << std::endl;
        // std::cout << "space debug string " << std::endl << space.DebugStr() << std::endl;
        // std::cout << "sectors debug string " << std::endl << sectors.DebugStr() << std::endl;
        // std::cout << "sector index " << sector_index << std::endl;
        // std::cout << "bra_subspace " << bra_subspace.LabelStr() << std::endl;
        // std::cout << "ket_subspace " << ket_subspace.LabelStr() << std::endl;
        // std::cout << "bra_subspace " << std::endl << bra_subspace.DebugStr() << std::endl;
        // std::cout << "ket_subspace " << std::endl << ket_subspace.DebugStr() << std::endl;
        // std::cout << "state_labels_bra " << std::get<0>(state_labels_bra) << " " << std::get<1>(state_labels_bra) << " " << std::get<2>(state_labels_bra) << " " << std::get<3>(state_labels_bra) << std::endl;
        // std::cout << "state_labels_ket " << std::get<0>(state_labels_ket) << " " << std::get<1>(state_labels_ket) << " " << std::get<2>(state_labels_ket) << " " << std::get<3>(state_labels_ket) << std::endl;
        // std::cout << "bra state index " << bra_state_index << std::endl;
        // std::cout << "ket state index " << ket_state_index << std::endl;
        // std::cout << "matrix element " << matrix_element << std::endl;
        std::cerr << "ERROR: Can't find indexes for a non-zero matrix element." << std::endl;
        exit(EXIT_FAILURE);
        // return;
      } else {
        // std::cout << "0" << std::endl;
        return;
      }
    }
    matrices[sector_index](bra_state_index,ket_state_index)=matrix_element;
  }

  double GetTwoBodyOperatorMatrixElementJJJTTz(
      const basis::TwoBodySpaceJJJTTz& space,
      const TwoBodySubspaceJJJTTzLabels& subspace_labels_bra,
      const TwoBodySubspaceJJJTTzLabels& subspace_labels_ket,
      const TwoBodyStateJJJTTzLabels& state_labels_bra,
      const TwoBodyStateJJJTTzLabels& state_labels_ket,
      const TwoBodySectorsJJJTTz& sectors,
      const OperatorBlocks<double>& matrices
    )
  {

    // look up sector
    std::size_t bra_subspace_index=space.LookUpSubspaceIndex(subspace_labels_bra);
    std::size_t ket_subspace_index=space.LookUpSubspaceIndex(subspace_labels_ket);
    if ((bra_subspace_index==kNone)||(ket_subspace_index==kNone))
      // gracefully handle lookup of forbidden subspace
      //
      // E.g., in ob-2 truncation, J=5, T=1, sector does not exist, since
      // (d5/2)^2 cannot couple to J+T even, but me2j I/O loop naively includes
      // this matrix element.
      return 0;
    std::size_t sector_index=sectors.LookUpSectorIndex(bra_subspace_index,ket_subspace_index);
    assert(sector_index!=kNone);

    // look up matrix element
    basis::TwoBodySectorsJJJTTz::SectorType sector=sectors.GetSector(sector_index);
    const basis::TwoBodySubspaceJJJTTz& bra_subspace=sector.bra_subspace();
    const basis::TwoBodySubspaceJJJTTz& ket_subspace=sector.ket_subspace();
    // std::cout << "bra_subspace " << bra_subspace.LabelStr() << std::endl;
    // std::cout << "ket_subspace " << ket_subspace.LabelStr() << std::endl;
    // std::cout << "bra_subspace " << std::endl << bra_subspace.DebugStr() << std::endl;
    // std::cout << "ket_subspace " << std::endl << ket_subspace.DebugStr() << std::endl;
    std::size_t bra_state_index=bra_subspace.LookUpStateIndex(state_labels_bra);
    std::size_t ket_state_index=ket_subspace.LookUpStateIndex(state_labels_ket);
    // std::cout << "bra state index " << bra_state_index << std::endl;
    // std::cout << "ket state index " << ket_state_index << std::endl;
    if ((bra_state_index==kNone)||(ket_state_index==kNone)) {
      // gracefully handle lookup of forbidden state
      return 0;
    }
    return matrices[sector_index](bra_state_index,ket_state_index);
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
