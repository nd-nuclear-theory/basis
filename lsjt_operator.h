/****************************************************************
  lsjt_operator.h

  Defines functions for I/O and manipulation of relative and two-body
  interaction matrices in LSJT coupling scheme.  Written for use in
  Moshinsky transformation.

  Language: C++11
  (for tuple and map::at)
                                 
  Mark A. Caprio, University of Notre Dame.

  11/21/15 (mac): Created (lsjt_interaction.h), from code in
    moshinsky_xform_lsjt.
  6/8/16 (mac): Update for basis package and update conventions.
    - Rename to lsjt_operator.h.
    - Remove matrix-style output.

****************************************************************/

#ifndef LSJT_OPERATOR_H_
#define LSJT_OPERATOR_H_

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include "eigen3/Eigen/Core"

#include "am/am.h"
#include "basis/lsjt_scheme.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // common typedefs
  ////////////////////////////////////////////////////////////////

  typedef std::vector< Eigen::MatrixXd > MatrixVector;


  ////////////////////////////////////////////////////////////////
  // templatized generic definitions
  ////////////////////////////////////////////////////////////////

  // TODO: move these generic definitions elsewhere!  not part of lsjt

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

  // identity operator operator

  template <typename tSectorsType>
    void SetOperatorToNaiveIdentityReduced(
        const tSectorsType& sectors,
        MatrixVector& matrices
      )
  // Provides spherical tensor reduced identity.
  //
  // CAVEAT: Assumes that the matrix elements are RMEs in a
  // *normalized* basis.  In particular, this is not applicable if basis is two-body
  // AS, as opposed to NAS, basis.
  //
  // tSpaceType::SubspaceType must have a J accessor.
  //
  // Assertion: Bra and ket J are same for all sectors.
  {

    // clear vector of matrices
    matrices.clear();

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename tSectorsType::SectorType& sector = sectors.GetSector(sector_index);
	const typename tSectorsType::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename tSectorsType::SubspaceType& ket_subspace = sector.ket_subspace();

        // verify sector is diagonal in J
	int J = ket_subspace.J();
	int J_check = bra_subspace.J();
	assert(J == J_check);

        // generate matrix for sector
	Eigen::MatrixXd sector_matrix;
	if (sector.bra_subspace_index() == sector.ket_subspace_index())
	  {
	    int J = ket_subspace.J();
	    sector_matrix = Hat(J) * Eigen::MatrixXd::Identity(bra_subspace.size(),ket_subspace.size());
	  }
	else
	  {
	    sector_matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());
	  }
	matrices.push_back(sector_matrix);
      }
  }

  ////////////////////////////////////////////////////////////////
  // relative space
  ////////////////////////////////////////////////////////////////

  // relative space output
  void WriteRelativeOperatorLSJT(
      std::ostream& os,
      const RelativeSectorsLSJT& sectors,
      const MatrixVector& matrices
    );

  // relative space input
  void ReadRelativeOperatorLSJT(
      std::istream& is,
      const RelativeSectorsLSJT& sectors,
      MatrixVector& matrices
    );

  ////////////////////////////////////////////////////////////////
  // two-body space
  ////////////////////////////////////////////////////////////////

  void WriteTwoBodyOperatorLSJT(
      std::ostream& os,
      const TwoBodySectorsLSJT& sectors,
      const MatrixVector& matrices
    );


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
