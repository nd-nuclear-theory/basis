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

#include "jjjt_operator.h"

#include "am/halfint.h"
#include <Eigen/Core>
#include "basis.h"
#include "many_body.h"
#include "jjjt_scheme.h"
#include "jt_operator.h"
#include "operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body JJJT operator -- gather N blocks
  ////////////////////////////////////////////////////////////////

  inline
  void RecastLabelsTwoBodyJJJTToTwoBodyJJJTN(
      const TwoBodySubspaceJJJTLabels& two_body_jjjt_subspace_labels,
      const TwoBodyStateJJJTLabels& two_body_jjjt_state_labels,
      TwoBodySubspaceJJJTNLabels& two_body_jjjtn_subspace_labels,
      TwoBodyStateJJJTNLabels& two_body_jjjtn_state_labels
    )
  // Recast labels for (subspace,state) from TwoBodyJJJT scheme to
  // TwoBodyJJJTN scheme.
  //
  // This direction of conversion (from target labels to source
  // labels) is as needed for looking up the source matrix element.
  //
  // This switches N from being the most-significant state label (that
  // is, least-rapidly-varying in the lexicographic ordering scheme),
  // to being the least-significant subspace label (that is,
  // most-rapidly-varying):
  //
  //   (L,S,J,T,g) : ([N],N1,l1,N2,l2) -> (L,S,J,T,g,N) : (N1,l1,N2,l2)
  //
  // As a state label, N is actually an implicit label, not stored,
  // but effectively the first label for ordering purposes.  Its value
  // may be recovered as N1+N2.
  //
  // Arguments:
  //   two_body_jjjt_subspace_labels (...) : source subspace labels
  //   two_body_jjjt_state_labels (...) : source state labels
  //   two_body_jjjtn_subspace_labels (...,output) : target subspace labels
  //   two_body_jjjtn_state_labels (...,output) : target state labels
  {
    // extract labels
    int J, T, g;
    std::tie(J,T,g) = two_body_jjjt_subspace_labels;
    int N1, l1, N2, l2;
    HalfInt j1, j2;
    std::tie(N1,j1,N2,j2) = two_body_jjjt_state_labels;

    // repackage labels
    int N = N1+N2;
    two_body_jjjtn_subspace_labels = TwoBodySubspaceJJJTNLabels(J,T,g,N);
    two_body_jjjtn_state_labels = two_body_jjjt_state_labels;
  }

  void GatherOperatorTwoBodyJJJTNToTwoBodyJJJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceJJJTN& two_body_jjjtn_space,
      const std::array<basis::TwoBodySectorsJJJTN,3>& two_body_jjjtn_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjtn_component_matrices,
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices
    )
  {
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      two_body_jjjt_component_sectors[T0]
        = basis::TwoBodySectorsJJJT(two_body_jjjt_space,operator_labels.J0,T0,operator_labels.g0);

      // populate matrices
      two_body_jjjt_component_matrices[T0].resize(two_body_jjjt_component_sectors[T0].size());
      for (std::size_t sector_index=0; sector_index<two_body_jjjt_component_sectors[T0].size(); ++sector_index)
        {
          // retrieve target sector
          const basis::TwoBodySectorsJJJT::SectorType& two_body_jjjt_sector
            = two_body_jjjt_component_sectors[T0].GetSector(sector_index);

          // initialize matrix
          Eigen::MatrixXd& two_body_jjjt_matrix = two_body_jjjt_component_matrices[T0][sector_index];
          two_body_jjjt_matrix = Eigen::MatrixXd::Zero(
              two_body_jjjt_sector.bra_subspace().size(),
              two_body_jjjt_sector.ket_subspace().size()
            );

          // populate matrix elements
          for (std::size_t bra_index = 0; bra_index < two_body_jjjt_sector.bra_subspace().size(); ++bra_index)
            for (std::size_t ket_index = 0; ket_index < two_body_jjjt_sector.ket_subspace().size(); ++ket_index)
              // for each target matrix element
              {

                // ensure canonical matrix element if diagonal sector
                if (two_body_jjjt_sector.IsDiagonal())
                  if (!(bra_index<=ket_index))
                    continue;

                // retrieve target states
                basis::TwoBodyStateJJJT two_body_jjjt_bra(two_body_jjjt_sector.bra_subspace(),bra_index);
                basis::TwoBodyStateJJJT two_body_jjjt_ket(two_body_jjjt_sector.ket_subspace(),ket_index);

                // extract source bra labels
                TwoBodySubspaceJJJTLabels two_body_jjjt_subspace_labels_bra
                  = two_body_jjjt_bra.subspace().labels();
                TwoBodyStateJJJTLabels two_body_jjjt_state_labels_bra
                  = two_body_jjjt_bra.labels();
                TwoBodySubspaceJJJTNLabels two_body_jjjtn_subspace_labels_bra;
                TwoBodyStateJJJTNLabels two_body_jjjtn_state_labels_bra;
                RecastLabelsTwoBodyJJJTToTwoBodyJJJTN(
                    two_body_jjjt_subspace_labels_bra,
                    two_body_jjjt_state_labels_bra,
                    two_body_jjjtn_subspace_labels_bra,
                    two_body_jjjtn_state_labels_bra
                  );

                // extract source bra indices
                std::size_t two_body_jjjtn_subspace_index_bra
                  = two_body_jjjtn_space.LookUpSubspaceIndex(
                      two_body_jjjtn_subspace_labels_bra
                    );
                std::size_t two_body_jjjtn_state_index_bra
                  = two_body_jjjtn_space.GetSubspace(two_body_jjjtn_subspace_index_bra).LookUpStateIndex(
                      two_body_jjjtn_state_labels_bra
                    );

                // extract source ket labels
                TwoBodySubspaceJJJTLabels two_body_jjjt_subspace_labels_ket
                  = two_body_jjjt_ket.subspace().labels();
                TwoBodyStateJJJTLabels two_body_jjjt_state_labels_ket
                  = two_body_jjjt_ket.labels();
                TwoBodySubspaceJJJTNLabels two_body_jjjtn_subspace_labels_ket;
                TwoBodyStateJJJTNLabels two_body_jjjtn_state_labels_ket;
                RecastLabelsTwoBodyJJJTToTwoBodyJJJTN(
                    two_body_jjjt_subspace_labels_ket,
                    two_body_jjjt_state_labels_ket,
                    two_body_jjjtn_subspace_labels_ket,
                    two_body_jjjtn_state_labels_ket
                  );

                // extract source ket indices
                std::size_t two_body_jjjtn_subspace_index_ket
                  = two_body_jjjtn_space.LookUpSubspaceIndex(
                      two_body_jjjtn_subspace_labels_ket
                    );
                std::size_t two_body_jjjtn_state_index_ket
                  = two_body_jjjtn_space.GetSubspace(two_body_jjjtn_subspace_index_ket).LookUpStateIndex(
                      two_body_jjjtn_state_labels_ket
                    );

                // look up matrix element
                std::size_t two_body_jjjtn_sector_index
                  = two_body_jjjtn_component_sectors[T0].LookUpSectorIndex(
                      two_body_jjjtn_subspace_index_bra,
                      two_body_jjjtn_subspace_index_ket
                    );

                const Eigen::MatrixXd& two_body_jjjtn_matrix
                  = two_body_jjjtn_component_matrices[T0][two_body_jjjtn_sector_index];
                double two_body_jjjtn_matrix_element = two_body_jjjtn_matrix(
                    two_body_jjjtn_state_index_bra,two_body_jjjtn_state_index_ket
                  );

                two_body_jjjt_matrix(bra_index,ket_index) = two_body_jjjtn_matrix_element;

              }
        }
    }

  }

  ////////////////////////////////////////////////////////////////
  // two-body JJJT operator output
  ////////////////////////////////////////////////////////////////

  void WriteTwoBodyOperatorComponentJJJT(
      std::ostream& os,
      int T0,
      const TwoBodySectorsJJJT& sectors,
      const OperatorBlocks<double>& matrices,
      basis::NormalizationConversion conversion_mode
    )
  {

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
        const typename TwoBodySectorsJJJT::SectorType& sector = sectors.GetSector(sector_index);
        const typename TwoBodySectorsJJJT::SubspaceType& bra_subspace = sector.bra_subspace();
        const typename TwoBodySectorsJJJT::SubspaceType& ket_subspace = sector.ket_subspace();

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
              const basis::TwoBodyStateJJJT bra(bra_subspace,bra_index);
              const basis::TwoBodyStateJJJT ket(ket_subspace,ket_index);

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
                << " " << std::setw(width) << T0
                << " " << "  "
                << " " << std::setw(width) << bra.N1()
                << " " << std::setw(width) << bra.l1()
                << " " << std::showpoint << std::fixed << std::setprecision(1) << std::setw(4) << float(bra.j1())
                << " " << std::setw(width) << bra.N2()
                << " " << std::setw(width) << bra.l2()
                << " " << std::showpoint << std::fixed << std::setprecision(1) << std::setw(4) << float(bra.j2())
                << " " << std::setw(width) << bra.J()
                << " " << std::setw(width) << bra.T()
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
                << " " << std::setw(width) << ket.g()
                << " " << "    "
                << " " << std::showpoint << std::scientific << std::setprecision(precision) << matrix_element
                << std::endl;

            }

      }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
