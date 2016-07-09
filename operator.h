/****************************************************************
  operator.h

  Provides generic definitions for storage of operator matrix elements
  making use of the basis package's indexing scheme and Eigen matrices.

  Libraries: Eigen 3

  Note: Eigen should be installed and the compiler include path such
  that the subdirectory chain "eigen3/Eigen/" can be found in a
  directory in the search path.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/2/16 (mac): Created from code in lsjt_operator.h.

****************************************************************/

#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <vector>

#include "eigen3/Eigen/Core"

#include "basis/indexing.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // matrix storage convenience typedef
  ////////////////////////////////////////////////////////////////

  typedef std::vector<Eigen::MatrixXd> MatrixVector;

  ////////////////////////////////////////////////////////////////
  // conversion flag enum
  ////////////////////////////////////////////////////////////////

  enum class NormalizationConversion {kNone, kASToNAS, kNASToAS};
  // Used as value for a mode parameter to requests on-the-fly
  // conversion between antisymmetrized (AS) and normalized
  // antisymmetrized (NAS) matrix elements on input/output.

  ////////////////////////////////////////////////////////////////
  // zero operator
  ////////////////////////////////////////////////////////////////

  template <typename tSectorsType>
  void SetOperatorToZero(
      const tSectorsType& sectors,
      MatrixVector& matrices
    )
    // Set operator to zero operator.
    //
    // The matrices for all sectors are set to zero.
    //
    // Arguments:
    //   sectors (BaseSectors daughter type) : the set of sectors on
    //     which the operator is defined
    //   matrices (MatrixVector, output) : matrices to hold operator

  {

    // clear vector of matrices
    matrices.clear();

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename tSectorsType::SectorType& sector = sectors.GetSector(sector_index);

        // extract sector subspaces
	const typename tSectorsType::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename tSectorsType::SubspaceType& ket_subspace = sector.ket_subspace();

        // generate matrix for sector
	Eigen::MatrixXd sector_matrix;
        sector_matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());

        // store matrix
	matrices.push_back(sector_matrix);

      }
  }


  ////////////////////////////////////////////////////////////////
  // identity operator
  ////////////////////////////////////////////////////////////////

  template <typename tSectorsType>
  void SetOperatorToIdentity(
      const tSectorsType& sectors,
      MatrixVector& matrices
    )
    // Set operator to "naive" identity operator.
    //
    // The matrices for off-diagonal sectors are set to zero
    // matrices, and the matrices for diagonal sectors are set to
    // identity matrices.
    //
    // This is a "naive" identity operator, in that it may or may
    // not actually represent the identity operator, depending upon
    // normalization conventions.  For instance, for matrix elements
    // taken in an orthonormal basis this *will* represent the
    // identity operator, and likewise for reduced matrix elements
    // under group theory (Rose) normalization conventions in an
    // orthonormal basis.  But dimension factors will be missing if
    // the operator is represented by its reduced matrix elements
    // under angular-momentum (Edmonds) normalization conventions.
    // And normalization factors will be missing if the basis consists
    // of, e.g., two-body antisymmetrized (AS), as opposed to
    // normalized antisymmetrized (NAS), states.
    //
    // Arguments:
    //   sectors (BaseSectors daughter type) : the set of sectors on
    //     which the operator is defined
    //   matrices (MatrixVector, output) : matrices to hold operator

  {

    // clear vector of matrices
    matrices.clear();

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename tSectorsType::SectorType& sector = sectors.GetSector(sector_index);

        // extract sector subspaces
	const typename tSectorsType::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename tSectorsType::SubspaceType& ket_subspace = sector.ket_subspace();

        // generate matrix for sector
	Eigen::MatrixXd sector_matrix;
	if (sector.IsDiagonal())
	  {
	    sector_matrix = Eigen::MatrixXd::Identity(bra_subspace.size(),ket_subspace.size());
	  }
	else
	  {
	    sector_matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());
	  }

        // store matrix
	matrices.push_back(sector_matrix);

      }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
