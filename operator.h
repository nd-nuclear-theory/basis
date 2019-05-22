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
  + 09/07/18 (pjf): Fix incomplete templatization of SetOperatorToDiagonalConstant.
  + 12/19/18 (pjf): Add operator arithmetic manipulation functions.
  + 02/01/19 (pjf): Update comment on OperatorLinearCombination.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.

****************************************************************/

#ifndef BASIS_OPERATOR_H_
#define BASIS_OPERATOR_H_

#include <cstddef>
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
    std::size_t UpperTriangularEntries(
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
      std::size_t total_entries = 0;
      for (std::size_t sector_index=0; sector_index<sectors.size(); ++sector_index)
        {
          // make reference to sector for convenience
          const typename tSectors::SectorType& sector
            = sectors.GetSector(sector_index);

          // count sector entries
          std::size_t sector_entries = 0;
          if (sector.IsDiagonal())
            // diagonal sector
            {
              std::size_t dimension = sector.ket_subspace().size();
              sector_entries = dimension*(dimension+1)/2;
            }
          else if (sector.IsUpperTriangle())
            // upper triangle sector (but not diagonal)
            {
              std::size_t bra_dimension = sector.bra_subspace().size();
              std::size_t ket_dimension = sector.ket_subspace().size();
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
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
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
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
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
      tFloat c
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
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
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
  // operator linear combination
  ////////////////////////////////////////////////////////////////

  // new features in C++14 allow us to manipulate the Eigen expressions
  // to get better performance; if decltype(auto) is not available, use
  // naive
  #if __cpp_decltype_auto
    template <typename tFloat>
    decltype(auto) OperatorBlockLinearCombination(
        const std::size_t& sector_index,
        const tFloat& a,
        const OperatorBlocks<tFloat>& matrices
      )
      // Base case for linear combination of block. Simply return
      // constant times matrix.
      //
      // WARNING: NOT INTENDED TO BE CALLED BY USER CODE.
      // See OperatorLinearCombination for further documentation.
    {
      return a*matrices[sector_index];
    }

    template <typename tFloat, typename... Args>
    decltype(auto) OperatorBlockLinearCombination(
        const std::size_t& sector_index,
        const tFloat& a,
        const OperatorBlocks<tFloat>& matrices,
        Args... args
      )
      // Peel off one layer of inputs from parameter pack and
      // add matrix to sum, recurse.
      //
      // WARNING: NOT INTENDED TO BE CALLED BY USER CODE.
      // See OperatorLinearCombination for further documentation.
    {
      return a*matrices[sector_index] + OperatorBlockLinearCombination(sector_index, args...);
    }
  #else
    template <typename tFloat>
    OperatorBlock<tFloat> OperatorBlockLinearCombination(
      const std::size_t& sector_index,
      const tFloat& a,
      const OperatorBlocks<tFloat>& matrices
      )
    // Base case for linear combination of block. Simply return
    // constant times matrix.
    //
    // WARNING: NOT INTENDED TO BE CALLED BY USER CODE.
    // See OperatorLinearCombination for further documentation.
    {
      return (a*matrices[sector_index]).eval();
    }

    template <typename tFloat, typename... Args>
    OperatorBlock<tFloat> OperatorBlockLinearCombination(
        const std::size_t& sector_index,
        const tFloat& a,
        const OperatorBlocks<tFloat>& matrices,
        Args... args
      )
    // Peel off one layer of inputs from parameter pack and
    // add matrix to sum, recurse.
    //
    // WARNING: NOT INTENDED TO BE CALLED BY USER CODE.
    // See OperatorLinearCombination for further documentation.
    {
      return (a*matrices[sector_index]
                + OperatorBlockLinearCombination(sector_index, args...)).eval();
    }
  #endif

  template <typename tSectorsType, typename tFloat, typename... Args>
  void OperatorLinearCombination(
      const tSectorsType& sectors,
      OperatorBlocks<tFloat>& output_matrices,
      Args... args
    )
    // Set output operator as linear combination of input operators.
    //
    // The output matrices are set to $\sum_i a_i O_i$. Any number of
    // input coefficients and operators matrices may be passed.
    //
    // This variadic template function uses parameter packs to allow
    // any arbitrary number of terms in the linear combination. On
    // modern compilers (with decltype(auto) defined) this allows
    // us to delay evaluation of Eigen's template language and optimize
    // matrix operations.
    //
    // Arguments:
    //   sectors (input): the set of sectors on
    //     which the operators are defined
    //   output_matrices (output): matrices to hold operator
    //   Args (input, parameter pack): pairs defined by
    //     a (input): coefficient of operator term
    //     matrices (input): matrices for operator term
  {
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
    {
      output_matrices[sector_index] = OperatorBlockLinearCombination(sector_index, args...);
    }
  }

  ////////////////////////////////////////////////////////////////
  // operator scalar multiply
  ////////////////////////////////////////////////////////////////

  template <typename tSectorsType, typename tFloat>
  void ScalarMultiplyOperator(
      const tSectorsType& sectors,
      basis::OperatorBlocks<tFloat>& matrices,
      const tFloat& c
    )
  // Multiply operator by a constant.
  //
  // All matrices are multiplied by a constant scalar value.
  //
  // Arguments:
  //   sectors (input): the set of sectors on
  //     which the operator is defined
  //   matrices (input/output): matrices to be multiplied
  //   c (input): the diagonal constant
  {
    OperatorLinearCombination(sectors, matrices, c, matrices);
  }

  ////////////////////////////////////////////////////////////////
  // operator addition
  ////////////////////////////////////////////////////////////////

  template <typename tSectorsType, typename tFloat>
  void AddOperators(
      const tSectorsType& sectors,
      const basis::OperatorBlocks<tFloat>& matrices_a,
      const basis::OperatorBlocks<tFloat>& matrices_b,
      basis::OperatorBlocks<tFloat>& output_matrices
    )
  // Add two operators and place into new output matrices.
  //
  // Loop over sectors, adding matrices from operators A and B and placing
  // in output matrix. This assumes that both operators share common sectors.
  //
  // Arguments:
  //   sectors (input): the set of sectors on
  //     which the operator is defined
  //   matrices_a (input): matrices for first operator
  //   matrices_b (input): matrices for second operator
  //   output_matrices (output): matrices for sum operator
  {
    OperatorLinearCombination(sectors, output_matrices, 1, matrices_a, 1, matrices_b);
  }

  ////////////////////////////////////////////////////////////////
  // operator addition
  ////////////////////////////////////////////////////////////////

  template <typename tSectorsType, typename tFloat>
  void SubtractOperators(
      const tSectorsType& sectors,
      const basis::OperatorBlocks<tFloat>& matrices_a,
      const basis::OperatorBlocks<tFloat>& matrices_b,
      basis::OperatorBlocks<tFloat>& output_matrices
    )
  // Subtract operator B from operator A and place into new output matrices.
  //
  // Arguments:
  //   sectors (input): the set of sectors on
  //     which the operator is defined
  //   matrices_a (input): matrices for first operator
  //   matrices_b (input): matrices for second operator
  //   output_matrices (output): matrices for sum operator
  {
    OperatorLinearCombination(sectors, output_matrices, 1, matrices_a, -1, matrices_b);
  }

  ////////////////////////////////////////////////////////////////
  // matrix element lookup
  ////////////////////////////////////////////////////////////////

  inline
    std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,bool>
    CanonicalizeIndices(
        std::size_t subspace_index_bra, std::size_t subspace_index_ket,
        std::size_t state_index_bra, std::size_t state_index_ket
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
      return std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,bool>(
          subspace_index_bra,subspace_index_ket,
          state_index_bra,state_index_ket,
          swapped_subspaces
        );
    }



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
