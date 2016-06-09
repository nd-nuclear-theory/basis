/****************************************************************

  lsjt_operator.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/


#include "lsjt_operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // relative space
  ////////////////////////////////////////////////////////////////

  void WriteRelativeOperatorLSJT(
      std::ostream& os,
      const RelativeSectorsLSJT& sectors,
      const MatrixVector& matrices
    )
  {

    // output formatting
    const int lw = 3;

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
	const typename RelativeSectorsLSJT::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename RelativeSectorsLSJT::SubspaceType& ket_subspace = sector.ket_subspace();

	// iterate over matrix elements
	for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
	  for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
	    {

              // define states
	      const basis::RelativeStateLSJT bra(bra_subspace,bra_index);
	      const basis::RelativeStateLSJT ket(ket_subspace,ket_index);
	    
              // generate output line
	      const double matrix_element = matrices[sector_index](bra_index,ket_index);
	      os 
		<< " " << std::setw(lw) << bra.Nr()
		<< " " << std::setw(lw) << bra.L() 
		<< " " << std::setw(lw) << bra.S() 
		<< " " << std::setw(lw) << bra.J() 
		<< " " << std::setw(lw) << bra.T() 
		<< " " << std::setw(lw) << bra.g()
		<< " " << "    "
		<< " " << std::setw(lw) << ket.Nr()
		<< " " << std::setw(lw) << ket.L() 
		<< " " << std::setw(lw) << ket.S() 
		<< " " << std::setw(lw) << ket.J() 
		<< " " << std::setw(lw) << ket.T() 
		<< " " << std::setw(lw) << ket.g()
		<< " " << "    "
		<< " " << std::showpoint << std::scientific << matrix_element
		<< std::endl;
	    
	    }

      };
  }

  void ReadRelativeOperatorLSJT(
      std::istream& is,
      const RelativeSectorsLSJT& sectors,
      MatrixVector& matrices
    )
  {

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
	const typename RelativeSectorsLSJT::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename RelativeSectorsLSJT::SubspaceType& ket_subspace = sector.ket_subspace();

        // construct zero matrix for sector
	Eigen::MatrixXd sector_matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());

	// retrieve matrix elements
	for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
	  for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
	    {

              // define states
	      const basis::RelativeStateLSJT bra(bra_subspace,bra_index);
	      const basis::RelativeStateLSJT ket(ket_subspace,ket_index);
	    
	      // read numbers from input line
	      int
 		bra_Nr,
		bra_L, 
		bra_S, 
		bra_J, 
		bra_T, 
		bra_g,
		ket_Nr,
		ket_L, 
		ket_S, 
		ket_J, 
		ket_T, 
		ket_g;
	      double matrix_element;
	      is 
		>> bra_Nr
		>> bra_L
		>> bra_S
		>> bra_J
		>> bra_T
		>> bra_g
		>> ket_Nr
		>> ket_L
		>> ket_S
		>> ket_J
		>> ket_T
		>> ket_g
		>> matrix_element;

	      // validate labels
	      bool expected_labels = true
		&& (bra_Nr==bra.Nr())
		&& (bra_L==bra.L()) 
		&& (bra_S==bra.S()) 
		&& (bra_J==bra.J()) 
		&& (bra_T==bra.T()) 
		&& (bra_g==bra.g())
		&& (ket_Nr==ket.Nr())
		&& (ket_L==ket.L()) 
		&& (ket_S==ket.S()) 
		&& (ket_J==ket.J()) 
		&& (ket_T==ket.T()) 
		&& (ket_g==ket.g());
	      assert(expected_labels);

              // save matrix element
	      sector_matrix(bra_index,ket_index) = matrix_element;
	    }

        // store matrix for sector
	matrices.push_back(sector_matrix);

      }
  }

  ////////////////////////////////////////////////////////////////
  // two-body space
  ////////////////////////////////////////////////////////////////


  void WriteTwoBodyOperatorLSJT(
      std::ostream& os,
      const TwoBodySectorsLSJT& sectors,
      const MatrixVector& matrices
    )
  {

    // output formatting
    const int lw = 3;

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename TwoBodySectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
	const typename TwoBodySectorsLSJT::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename TwoBodySectorsLSJT::SubspaceType& ket_subspace = sector.ket_subspace();

	// iterate over matrix elements
	for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
	  for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
	    {

              // define states
	      const basis::TwoBodyStateLSJT bra(bra_subspace,bra_index);
	      const basis::TwoBodyStateLSJT ket(ket_subspace,ket_index);
	    
	      const double matrix_element = matrices[sector_index](bra_index,ket_index);
	      os 
		<< " " << std::setw(lw) << bra.N1()
		<< " " << std::setw(lw) << bra.l1()
		<< " " << std::setw(lw) << bra.N2()
		<< " " << std::setw(lw) << bra.l2()
		<< " " << std::setw(lw) << bra.L() 
		<< " " << std::setw(lw) << bra.S() 
		<< " " << std::setw(lw) << bra.J() 
		<< " " << std::setw(lw) << bra.T() 
		<< " " << std::setw(lw) << bra.g()
		<< " " << "    "
		<< " " << std::setw(lw) << ket.N1()
		<< " " << std::setw(lw) << ket.l1()
		<< " " << std::setw(lw) << ket.N2()
		<< " " << std::setw(lw) << ket.l2()
		<< " " << std::setw(lw) << ket.L() 
		<< " " << std::setw(lw) << ket.S() 
		<< " " << std::setw(lw) << ket.J() 
		<< " " << std::setw(lw) << ket.T() 
		<< " " << std::setw(lw) << ket.g()
		<< " " << "    "
		<< " " << std::showpoint << std::scientific << matrix_element
		<< std::endl;
	    
	    }

      };
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
