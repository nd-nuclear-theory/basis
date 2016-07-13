/****************************************************************
  operator.h

  Provides generic definitions for storage of operator matrix elements
  making use of the basis package's indexing scheme and Eigen matrices.

  Libraries: Eigen 3

  Note: Eigen should be installed, and the compiler include path
  should be set ("-I") to include the appropriate *parent* directory
  such that the subdirectory chain "eigen3/Eigen/" stems from this
  directory.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/2/16 (mac): Created from code in lsjt_operator.h.
  7/9/16 (mac): Add support for canonicalizing indices in matrix
    element lookup.
  7/12/16 (mac): Add diagonal constant operator.

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
    // RELATIONS: This could be absorbed as a special case of
    // SetOperatorToDiagonalConstant.
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
    // RELATIONS: This could be absorbed as a special case of
    // SetOperatorToDiagonalConstant.
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
  // constant operator
  ////////////////////////////////////////////////////////////////

  template <typename tSectorsType>
  void SetOperatorToDiagonalConstant(
      const tSectorsType& sectors,
      MatrixVector& matrices,
      double c
    )
    // Set operator to "naive" diagonal constant operator.
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
    //   c (double) : the diagonal constant

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
	    sector_matrix = c*Eigen::MatrixXd::Identity(bra_subspace.size(),ket_subspace.size());
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
  // matrix element lookup
  ////////////////////////////////////////////////////////////////

  inline
  void CanonicalizeIndices(
      int& bra_subspace_index, int& ket_subspace_index,
      bool& swapped_subspaces, 
      int& bra_state_index, int& ket_state_index,
      bool& swapped_states
    )
  // Convert subspace and state indices for a matrix element to
  // canonical ("upper triangle") indices.
  //
  // Note: For specific indexing schemes, it might be more convenient
  // to define a customized "canonicalizaton" function, which, rather
  // than simply flagging the swaps through boolean variable, actually
  // calculates any necessary phase an dimension factors arising from
  // the swaps.
  //
  // Arguments:
  //    bra_subspace_index, ket_subspace_index (int, input/output) :
  //      sector bra and ket subspace indices, possibly to be swapped
  //    swapped_subspaces (bool, output) : whether or not subspace
  //      indices were swapped (off-diagonal sectors)
  //    bra_state_index, ket_state_index (int, input/output) :
  //      bra and ket state indices, possibly to be swapped if sector
  //      is diagonal sector
  //    swapped_states (bool, output) : whether or not state
  //      indices were swapped (diagonal sectors)
  {
    // process subspace indices (off-diagonal sectors)
    swapped_subspaces = !(bra_subspace_index <= ket_subspace_index);
    if (swapped_subspaces)
      {
        std::swap(bra_subspace_index,ket_subspace_index);
        std::swap(bra_state_index,ket_state_index);  // so index stays with sector
      }

    // process state indices (diagonal sectors)
    swapped_states = (bra_subspace_index == ket_subspace_index) 
      & !(bra_state_index <= ket_state_index);
    if (swapped_states)
      std::swap(bra_state_index,ket_state_index);
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
