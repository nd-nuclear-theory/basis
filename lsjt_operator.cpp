/****************************************************************

  lsjt_operator.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/

#include "lsjt_operator.h"

namespace basis {

  void WriteRelativeOperatorParametersLSJT(
      std::ostream& os,
      const basis::RelativeOperatorParametersLSJT& parameters
    )
  {
    assert(parameters.version==1);
    os 
      << "# RELATIVE LSJT" << std::endl
      << " " << parameters.version << std::endl
      << " " << parameters.J0 << " " << parameters.g0 << std::endl
      << " " << parameters.T0_min << " " << parameters.T0_max << std::endl
      << " " << parameters.Nmax << " " << parameters.Jmax << std::endl;
  }

  void ReadRelativeOperatorParametersLSJT(
      std::istream& is,
      basis::RelativeOperatorParametersLSJT& parameters
    )
  {

    std::string line;

    // line 1: version -- and gobble any comment lines
    while (std::getline(is,line), line[0]=='#') {};
    std::stringstream(line) >> parameters.version;
    assert(parameters.version==1);
    
    // line 2: operator tensor properties
    std::getline(is,line);
    std::stringstream(line) >> parameters.J0 >> parameters.g0;

    // line 3: isospin components
    std::getline(is,line);
    std::stringstream(line) >> parameters.T0_min >> parameters.T0_max;

    // line 4: truncation
    std::getline(is,line);
    std::stringstream(line) >> parameters.Nmax >> parameters.Jmax;
  }

  void WriteRelativeOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const RelativeSectorsLSJT& sectors,
      const MatrixVector& matrices
    )
  {

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

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // define states
	      const basis::RelativeStateLSJT bra(bra_subspace,bra_index);
	      const basis::RelativeStateLSJT ket(ket_subspace,ket_index);
	    
              // generate output line
              const int width = 3;
	      const double matrix_element = matrices[sector_index](bra_index,ket_index);
	      os 
		<< " " << std::setw(width) << T0
		<< " " << "  "
		<< " " << std::setw(width) << bra.N()
		<< " " << std::setw(width) << bra.L() 
		<< " " << std::setw(width) << bra.S() 
		<< " " << std::setw(width) << bra.J() 
		<< " " << std::setw(width) << bra.T() 
		<< " " << "  "
		<< " " << std::setw(width) << ket.N()
		<< " " << std::setw(width) << ket.L() 
		<< " " << std::setw(width) << ket.S() 
		<< " " << std::setw(width) << ket.J() 
		<< " " << std::setw(width) << ket.T() 
		<< " " << "  "
		<< " " << std::showpoint << std::scientific << matrix_element
		<< std::endl;
	    
	    }

      };
  }

  void ReadRelativeOperatorComponentLSJT(
      std::istream& is,
      int T0,
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

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // define states
	      const basis::RelativeStateLSJT bra(bra_subspace,bra_index);
	      const basis::RelativeStateLSJT ket(ket_subspace,ket_index);
	    
	      // read numbers from input line
	      int
                input_T0,
 		bra_N,
		bra_L, 
		bra_S, 
		bra_J, 
		bra_T, 
		ket_N,
		ket_L, 
		ket_S, 
		ket_J, 
		ket_T;
	      double matrix_element;
	      is 
                >> input_T0
		>> bra_N
		>> bra_L
		>> bra_S
		>> bra_J
		>> bra_T
		>> ket_N
		>> ket_L
		>> ket_S
		>> ket_J
		>> ket_T
		>> matrix_element;

	      // validate labels
	      bool expected_labels = true
		&& (input_T0==T0)
		&& (bra_N==bra.N())
		&& (bra_L==bra.L())
		&& (bra_S==bra.S())
		&& (bra_J==bra.J())
		&& (bra_T==bra.T())
		&& (ket_N==ket.N())
		&& (ket_L==ket.L())
		&& (ket_S==ket.S())
		&& (ket_J==ket.J())
		&& (ket_T==ket.T());
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
	    
              const int width = 3;
              const double matrix_element = matrices[sector_index](bra_index,ket_index);
	      os 
		<< " " << std::setw(width) << bra.N1()
		<< " " << std::setw(width) << bra.l1()
		<< " " << std::setw(width) << bra.N2()
		<< " " << std::setw(width) << bra.l2()
		<< " " << std::setw(width) << bra.L() 
		<< " " << std::setw(width) << bra.S() 
		<< " " << std::setw(width) << bra.J() 
		<< " " << std::setw(width) << bra.T() 
		<< " " << std::setw(width) << bra.g()
		<< " " << "    "
		<< " " << std::setw(width) << ket.N1()
		<< " " << std::setw(width) << ket.l1()
		<< " " << std::setw(width) << ket.N2()
		<< " " << std::setw(width) << ket.l2()
		<< " " << std::setw(width) << ket.L() 
		<< " " << std::setw(width) << ket.S() 
		<< " " << std::setw(width) << ket.J() 
		<< " " << std::setw(width) << ket.T() 
		<< " " << std::setw(width) << ket.g()
		<< " " << "    "
		<< " " << std::showpoint << std::scientific << matrix_element
		<< std::endl;
	    
	    }

      };
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace