/****************************************************************
  lsjt_operator_test.cpp

  Mark A. Caprio, University of Notre Dame.
  See interaction_lsjt.h for documentation.

****************************************************************/

#include <iomanip>
#include <sstream>

#include "lsjt_operator.h"

////////////////////////////////////////////////////////////////
// templatized generic definitions
////////////////////////////////////////////////////////////////

// TODO: move these generic definitions elsewhere!  not part of lsjt

namespace basis {

  // output
  //
  // with generic index labels only... maybe someday to generalize using
  // a label member function on the subspace

  template <typename tSectorsType>
  void WriteOperator(
      std::ostream& os,
      const tSectorsType& sectors,
      const MatrixVector& matrices
    )
  {

    // output formatting
    const int lw = 3;

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename tSectorsType::SectorType& sector = sectors.GetSector(sector_index);
	const typename tSectorsType::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename tSectorsType::SubspaceType& ket_subspace = sector.ket_subspace();

        // write sector indices
	const int bra_subspace_index = sector.bra_subspace_index();
	const int ket_subspace_index = sector.ket_subspace_index();
	os 
	  << "sector"
          << " " << std::setw(lw) << sector_index
	  << " " << std::setw(lw) << bra_subspace_index
	  << " " << std::setw(lw) << ket_subspace_index
	  << std::endl;

	// iterate over matrix elements
	for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
	  for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
	    {
	      const double matrix_element = matrices[sector_index](bra_index,ket_index);
	      os 
		<< "  " 
                << std::setw(lw) << bra_index
		<< " " << std::setw(lw) << ket_index
		<< " " << std::showpoint << std::scientific << matrix_element
		<< std::endl;
	    
	    }

      };
  }



  // zero operator

  template <typename tSectorsType>
  void SetOperatorToZero(
      const tSectorsType& sectors,
      MatrixVector& matrices
    )
  {

    // clear vector of matrices
    matrices.clear();

    // iterator over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename tSectorsType::SectorType& sector = sectors.GetSector(sector_index);
	const typename tSectorsType::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename tSectorsType::SubspaceType& ket_subspace = sector.ket_subspace();

        // construct zero matrix
	Eigen::MatrixXd sector_matrix;
	sector_matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());
	matrices.push_back(sector_matrix);
      }
  }

} // namespace


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // relative state tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {

      std::cout << "Relative space" << std::endl;
      int Nr_max = 2;
      int Jr_max = 3;
      basis::RelativeSpaceLSJT space(Nr_max,Jr_max);
      int J0 = 0;  // also try J0=2 for quadrupole operator
      int T0 = 0;
      int g0 = 0;
      basis::RelativeSectorsLSJT sectors(space,J0,T0,g0);
      basis::MatrixVector matrices;

      std::cout << "Zeros" << std::endl;
      basis::SetOperatorToZero(sectors,matrices);
      basis::WriteRelativeOperatorLSJT(std::cout,sectors,matrices);

      std::cout << "Identity" << std::endl;
      // TO REPLACE: basis::SetOperatorToNaiveIdentityReduced(sectors,matrices);
      basis::WriteRelativeOperatorLSJT(std::cout,sectors,matrices);
      // std::cout << "...printed as matrices..." << std::endl;
      // WriteRelativeOperatorMatrices(std::cout,sectors,matrices,6,3);

      // writeout/readback test
      std::cout << "Readback" << std::endl;

      std::ostringstream os;
      basis::WriteRelativeOperatorLSJT(os,sectors,matrices);
      // std::cout << os.str() << std::endl; // debugging: inspect stream contents

      std::istringstream is(os.str());
      basis::MatrixVector matrices2;
      basis::ReadRelativeOperatorLSJT(is,sectors,matrices2);
      basis::WriteRelativeOperatorLSJT(std::cout,sectors,matrices2);

    }

  ////////////////////////////////////////////////////////////////
  // two-body state tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {

      std::cout << "Two-body space" << std::endl;
      int Nmax = 2;
      basis::TwoBodySpaceLSJT space(Nmax);
      int J0 = 0;  // also try J0=2 for quadrupole operator
      int T0 = 0;
      int g0 = 0;
      basis::TwoBodySectorsLSJT sectors(space,J0,T0,g0);
      basis::MatrixVector matrices;

      std::cout << "Identity" << std::endl;
      // TO REPLACE: basis::SetOperatorToNaiveIdentityReduced(sectors,matrices);
      basis::WriteTwoBodyOperatorLSJT(std::cout,sectors,matrices);
      //std::cout << "...printed as matrices..." << std::endl;
      //basis::WriteTwoBodyOperatorMatricesLSJT(std::cout,sectors,matrices,9,6);

    }


  // termination
  return 0;
}
