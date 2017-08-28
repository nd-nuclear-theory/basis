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
  + 4/11/17 (mac):
    - Replace MatrixVector (hard coded double-precision) with
      templatized OperatorBlocks (generic precision).
    - Retemplatize functions to work with generic precision
      operator blocks instead of MatrixXd.
  + 08/11/17 (pjf): Emit warnings if deprecated member functions are used.

****************************************************************/

#ifndef BASIS_OPERATOR_H_
#define BASIS_OPERATOR_H_

#include <vector>

#include "eigen3/Eigen/Dense"

#include "basis/basis.h"

// emit warnings on deprecated
#include "mcutils/deprecated.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // matrix storage convenience typedef
  ////////////////////////////////////////////////////////////////

  // typedefs for operator blocks
  //
  // EX: basis::OperatorBlock<double>  // single matrix
  // EX: basis::OperatorBlocks<double>  // matrices for all sectors

  // Note: Typedefs cannot be templatized, but C++11 introduced
  // parametrized type aliases, using the "using" keyword.
  //
  // template <typename tFloat>
  //   typedef Eigen::Matrix<tFloat,Eigen::Dynamic,Eigen::Dynamic> OperatorBlocks;  // WRONG!

  // Note: When templatizing a function which takes an OperatorBlocks
  // argument, the following idiom should be successful
  //
  //  template <typename tFloat>
  //    MyFunction(basis::OperatorBlocks<tFloat>& matrices)

  template <typename tFloat>
    using OperatorBlock = Eigen::Matrix<tFloat,Eigen::Dynamic,Eigen::Dynamic>;

  template <typename tFloat>
    using OperatorBlocks = std::vector<OperatorBlock<tFloat>>;

  // legacy typedef -- DEPRECATED
  DEPRECATED("use basis::OperatorBlocks<double> instead")
  typedef basis::OperatorBlocks<double> MatrixVector;


  ////////////////////////////////////////////////////////////////
  // storage diagnostics
  ////////////////////////////////////////////////////////////////
  template <typename tFloat>
    std::size_t AllocatedEntries(const basis::OperatorBlocks<tFloat>& matrices)
    // Count entries in a vector of matrices.
    //
    // This may be larger than the nominal number of matrix elements
    // defining the operator, if only upper-triangle matrix elements
    // are relevant for the operator, since full square/rectangular
    // matrices are always allocated.
    //
    // Arguments:
    //   matrices (input): vector of matrices
    //
    // Returns:
    //   (std::size_t): number of allocated matrix entries
    {
      std::size_t count = 0;
      for (auto iterator = matrices.begin(); iterator != matrices.end(); ++iterator)
        count += iterator->size();
      return count;
    }

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
    //   sectors (input): container for sectors
    //
    // Returns:
    //   number of upper triangular matrix entries
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

  template <typename tSectorsType, typename tFloat>
  void SetOperatorToZero(
      const tSectorsType& sectors,
      basis::OperatorBlocks<tFloat>& matrices
    )
    // Set operator blocks to zero.
    //
    // Arguments:
    //   sectors (input): the set of sectors on which the operator is
    //     defined
    //   matrices (output): matrices to hold blocks
  {

    // clear vector of matrices
    matrices.clear();
    matrices.resize(sectors.size());

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename tSectorsType::SectorType& sector = sectors.GetSector(sector_index);

        // extract sector subspaces
	const typename tSectorsType::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename tSectorsType::SubspaceType& ket_subspace = sector.ket_subspace();

        // generate matrix for sector
        matrices[sector_index] = basis::OperatorBlock<tFloat>::Zero(bra_subspace.size(),ket_subspace.size());

      }
  }

  ////////////////////////////////////////////////////////////////
  // identity operator
  ////////////////////////////////////////////////////////////////

  template <typename tSectorsType, typename tFloat>
  void SetOperatorToIdentity(
      const tSectorsType& sectors,
      basis::OperatorBlocks<tFloat>& matrices
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
    //   sectors (input): the set of sectors on
    //     which the operator is defined
    //   matrices (output): matrices to hold operator

  {

    // clear vector of matrices
    matrices.clear();
    matrices.resize(sectors.size());

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename tSectorsType::SectorType& sector = sectors.GetSector(sector_index);

        // extract sector subspaces
	const typename tSectorsType::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename tSectorsType::SubspaceType& ket_subspace = sector.ket_subspace();

        // generate matrix for sector
	if (sector.IsDiagonal())
	  {
            matrices[sector_index] = basis::OperatorBlock<tFloat>::Identity(bra_subspace.size(),ket_subspace.size());
	  }
	else
	  {
            matrices[sector_index] = basis::OperatorBlock<tFloat>::Zero(bra_subspace.size(),ket_subspace.size());
	  }

      }
  }

  ////////////////////////////////////////////////////////////////
  // constant operator
  ////////////////////////////////////////////////////////////////

  template <typename tSectorsType, typename tFloat>
  void SetOperatorToDiagonalConstant(
      const tSectorsType& sectors,
      basis::OperatorBlocks<tFloat>& matrices,
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
    //   sectors (input): the set of sectors on
    //     which the operator is defined
    //   matrices (output): matrices to hold operator
    //   c (input): the diagonal constant

  {

    // clear vector of matrices
    matrices.clear();
    matrices.resize(sectors.size());

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename tSectorsType::SectorType& sector = sectors.GetSector(sector_index);

        // extract sector subspaces
	const typename tSectorsType::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename tSectorsType::SubspaceType& ket_subspace = sector.ket_subspace();

        // generate matrix for sector
	if (sector.IsDiagonal())
	  {
            matrices[sector_index] = c*basis::OperatorBlock<tFloat>::Identity(bra_subspace.size(),ket_subspace.size());
	  }
	else
	  {
            matrices[sector_index] = basis::OperatorBlock<tFloat>::Zero(bra_subspace.size(),ket_subspace.size());
	  }

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
