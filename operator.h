/**************************************************************************//**
  @file operator.h

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

  + 7/2/16 (mac): Created from code in lsjt_operator.h.
  + 7/9/16 (mac): Add support for canonicalizing indices in matrix
    element lookup.
  + 7/12/16 (mac): Add diagonal constant operator.
  + 7/19/16 (mac):
    - Extract AS/NAS conversion enum to many_body.h.
    - Add diagnostic function AllocatedEntries.
  + 7/22/16 (mac): Revise syntax for CanonicalizeIndices.
  + 7/25/16 (mac): Add diagnostic function UpperTriangularEntries.
  + 11/1/16 (mac): Reduce dependency from Eigen/Core to Eigen/Dense.

****************************************************************/

#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <vector>

#include "eigen3/Eigen/Dense"

#include "basis/basis.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // matrix storage convenience typedef
  ////////////////////////////////////////////////////////////////

  typedef std::vector<Eigen::MatrixXd> MatrixVector;

  ////////////////////////////////////////////////////////////////
  // storage diagnostics
  ////////////////////////////////////////////////////////////////
  std::size_t AllocatedEntries(const MatrixVector& matrices);
    // Count entries in a vector of matrices.
    //
    // This may be larger than the nominal number of matrix elements
    // defining the operator, if only upper-triangle matrix elements
    // are relevant for the operator, since full square/rectangular
    // matrices are always allocated.
    //
    // Arguments:
    //   matrices (basis::MatrixVector) : vector of matrices
    //
    // Returns:
    //   (std::size_t) : number of allocated matrix entries


  template <typename tSectors>
    int UpperTriangularEntries(
        const tSectors& sectors
      )
    // Count entries in the upper triangular portion of a set of sectors.
    //
    // Lower triangular sectors are ignored.  In diagonal sectors,
    // only upper-triangular entries are counted.
    //
    // Arguments:
    //   sectors (tSectors) : container for sectors
    //
    // Returns:
    //   (std::size_t) : number of upper triangular matrix entries
    {
      int total_entries = 0;
      for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
        {
          // make reference to sector for convenience
          const typename tSectors::SectorType& sector
            = sectors.GetSector(sector_index);

          // count sector entries
          int sector_entries = 0;
          if (sector.IsDiagonal())
            // diagonal sector
            {
              int dimension = sector.ket_subspace().size();
              sector_entries = dimension*(dimension+1)/2;
            }
          else if (sector.IsUpperTriangle())
            // upper triangle sector (but not diagonal)
            {
              int bra_dimension = sector.bra_subspace().size();
              int ket_dimension = sector.ket_subspace().size();
              sector_entries = bra_dimension*ket_dimension;
            }

          total_entries += sector_entries;
        }

      return total_entries;
    }

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
    std::tuple<int,int,int,int,bool>
    CanonicalizeIndices(
        int subspace_index_bra, int subspace_index_ket,
        int state_index_bra, int state_index_ket
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
    //   subspace_index_bra, subspace_index_ket (input):
    //     naive sector bra and ket subspace indices, possibly to be swapped
    //   state_index_bra, state_index_ket (input):
    //     naive bra and ket state indices, possibly to be swapped if sector
    //     is diagonal
    //
    // Returns:
    //   canonicalized indices and swap flag as:
    //
    //        subspace_index_bra,subspace_index_ket,
    //        state_index_bra,state_index_ket,
    //        swapped_subspaces
    {
      // Note: We use the local copies of the indices (on the stack)
      // as working variables.  Slightly unnerving.

      // process subspace indices (off-diagonal sectors)
      bool swapped_subspaces = !(subspace_index_bra <= subspace_index_ket);
      if (swapped_subspaces)
        {
          std::swap(subspace_index_bra,subspace_index_ket);
          std::swap(state_index_bra,state_index_ket);  // so index stays with sector
        }

      // process state indices (diagonal sectors)
      bool swapped_states = (subspace_index_bra == subspace_index_ket)
        & !(state_index_bra <= state_index_ket);
      if (swapped_states)
        std::swap(state_index_bra,state_index_ket);

      // bundle return values
      return std::tuple<int,int,int,int,bool>(
          subspace_index_bra,subspace_index_ket,
          state_index_bra,state_index_ket,
          swapped_subspaces
        );
    }



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
