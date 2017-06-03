/************************************************************//**
  hypersector.h

  Enumeration of sectors for multiple operators.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 4/4/17 (mac): Created, building on code from basis.h.
  + 4/11/17 (mac): Add basic hyperoperator support --
    OperatorHypersectors and SetHyperoperatorToZero.
  + 4/22/17 (aem): Fix error in hypersector constructor.
****************************************************************/

#ifndef BASIS_HYPERSECTOR_H_
#define BASIS_HYPERSECTOR_H_

#include "basis/basis.h"
#include "basis/operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // hypersector indexing
  ////////////////////////////////////////////////////////////////

  // Note: A "hypersector" represents a triplet consisting of a pair
  // of "state" subspaces (which thus define a "sector" or block in
  // the matrix representation of an operator on the space) and an
  // "operator" subspace.  The hypersector may also be labeled with an
  // optional multiplicity index.

  template <typename tSubspaceType, typename tOperatorSubspaceType>
    class BaseHypersector
    // Store indexing and subspace reference information for a
    // hypersector.
    //
    // Allows for multiplicity index on sector, for symmetry sectors
    // of groups with outer multiplicities.
    {

      ////////////////////////////////////////////////////////////////
      // typedefs
      ////////////////////////////////////////////////////////////////

      public:
      typedef tSubspaceType SubspaceType;
      typedef tOperatorSubspaceType OperatorSubspaceType;
      typedef std::tuple<int,int,int,int> KeyType;

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseHypersector(
          int bra_subspace_index, int ket_subspace_index, int operator_subspace_index,
          const SubspaceType& bra_subspace, const SubspaceType& ket_subspace, const OperatorSubspaceType& operator_subspace,
          int multiplicity_index=1
        )
        : bra_subspace_index_(bra_subspace_index), ket_subspace_index_(ket_subspace_index),
        operator_subspace_index_(operator_subspace_index), multiplicity_index_(multiplicity_index)
      {
        bra_subspace_ptr_ = &bra_subspace;
        ket_subspace_ptr_ = &ket_subspace;
        operator_subspace_ptr_ = &operator_subspace;
      }

      ////////////////////////////////////////////////////////////////
      // accessors
      ////////////////////////////////////////////////////////////////

      inline KeyType Key() const
      // Return tuple key identifying sector for sorting/lookup
      // purposes.
      {
        return KeyType(bra_subspace_index(),ket_subspace_index(),operator_subspace_index(),multiplicity_index());
      }

      int bra_subspace_index() const {return bra_subspace_index_;}
      int ket_subspace_index() const {return ket_subspace_index_;}
      int operator_subspace_index() const {return operator_subspace_index_;}
      // Return integer index of bra/ket subspace.

      const SubspaceType& bra_subspace() const {return *bra_subspace_ptr_;}
      const SubspaceType& ket_subspace() const {return *ket_subspace_ptr_;}
      const OperatorSubspaceType& operator_subspace() const {return *operator_subspace_ptr_;}
      // Return reference to bra/ket/operator subspace.

      int multiplicity_index() const {return multiplicity_index_;}
      // Return multiplicity index of this sector.

      inline bool IsDiagonal() const
      // Test if sector is diagonal (i.e., within a single subspace).
      {
        return (bra_subspace_index()==ket_subspace_index());
      }

      inline bool IsUpperTriangle() const
      // Test if sector is in upper triangle (includes diagonal).
      {
        return (bra_subspace_index()<=ket_subspace_index());
      }

      private:
      int bra_subspace_index_, ket_subspace_index_, operator_subspace_index_;
      const SubspaceType* bra_subspace_ptr_;
      const SubspaceType* ket_subspace_ptr_;
      const OperatorSubspaceType* operator_subspace_ptr_;
      int multiplicity_index_;
    };

  // BaseHypersectors -- container to hold a set of hypersectors with
  // reverse lookup by hypersector labels
  //
  // Template arguments:
  //   tSpaceType (typename): type for bra/ket space
  //   tOperatorSpaceType (typename): type for operator space

  template <typename tSpaceType, typename tOperatorSpaceType>
    class BaseHypersectors
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      typedef tSpaceType SpaceType;
      typedef tOperatorSpaceType OperatorSpaceType;
      typedef typename tSpaceType::SubspaceType SubspaceType;
      typedef typename tOperatorSpaceType::SubspaceType OperatorSubspaceType;
      typedef BaseHypersector<SubspaceType,OperatorSubspaceType> HypersectorType;

      ////////////////////////////////////////////////////////////////
      // sector lookup and retrieval
      ////////////////////////////////////////////////////////////////

      const HypersectorType& GetHypersector(int hypersector_index) const
      // Given sector index, return reference to sector itself.
      {
        return hypersectors_[hypersector_index];
      };

      bool ContainsHypersector(
          int bra_subspace_index, int ket_subspace_index, int operator_subspace_index,
          int multiplicity_index=1
        ) const
      // Given the labels for a sector, returns whether or not the sector
      // is found within the the sector set.
      {
        typename HypersectorType::KeyType
          key(bra_subspace_index,ket_subspace_index,operator_subspace_index,multiplicity_index);
        return lookup_.count(key);
      };

      int LookUpHypersectorIndex(
          int bra_subspace_index, int ket_subspace_index, int operator_subspace_index,
          int multiplicity_index=1
        ) const
      // Given the labels for a sector, look up its index within the
      // sector set.
      //
      // If no such labels are found, basis::kNone is returned.
      {
        typename HypersectorType::KeyType
          key(bra_subspace_index,ket_subspace_index,operator_subspace_index,multiplicity_index);
        auto pos = lookup_.find(key);
        if (pos==lookup_.end())
          return kNone;
        else
          return pos->second;
      };

      ////////////////////////////////////////////////////////////////
      // size retrieval
      ////////////////////////////////////////////////////////////////

      int size() const
      // Return number of hypersectors within hypersector set.
      {
        return hypersectors_.size();
      };

      ////////////////////////////////////////////////////////////////
      // diagnostic strings
      ////////////////////////////////////////////////////////////////

      std::string DebugStr() const;
      // Generate string dump of contents, for debugging purposes.
      //
      // Requires subspace to have a LabelStr() member function.

      protected:

      ////////////////////////////////////////////////////////////////
      // hypersector push (for initial construction)
      ////////////////////////////////////////////////////////////////

      void PushHypersector(const HypersectorType& hypersector)
      // Create indexing information (in both directions, index <->
      // labels) for a sector.
      {
        lookup_[hypersector.Key()] = hypersectors_.size(); // index for lookup
        hypersectors_.push_back(hypersector);  // save hypersector
      };

      ////////////////////////////////////////////////////////////////
      // internal storage
      ////////////////////////////////////////////////////////////////

      // Hypersectors (accessible by index)
      std::vector<HypersectorType> hypersectors_;

      // sector index lookup by subspace indices
#ifdef BASIS_HASH
      std::unordered_map<typename HypersectorType::KeyType,int,boost::hash<typename HypersectorType::KeyType>> lookup_;
#else
      std::map<typename HypersectorType::KeyType,int> lookup_;
#endif

    };

  template <typename tSpaceType, typename tOperatorSpaceType>
    std::string BaseHypersectors<tSpaceType,tOperatorSpaceType>::DebugStr() const
    {
      std::ostringstream os;
      for (int hypersector_index=0; hypersector_index<size(); ++hypersector_index)
        {
          const HypersectorType& hypersector = GetHypersector(hypersector_index);
          
          os << "  hypersector " << hypersector_index
             << "  bra index " << hypersector.bra_subspace_index()
             << " labels " << hypersector.bra_subspace().LabelStr()
             << " dim " << hypersector.bra_subspace().size()
             << "  ket index " << hypersector.ket_subspace_index()
             << " labels " << hypersector.ket_subspace().LabelStr()
             << " dim " << hypersector.ket_subspace().size()
             << "  operator index " << hypersector.operator_subspace_index()
             << " labels " << hypersector.operator_subspace().LabelStr()
             << " dim " << hypersector.operator_subspace().size()
             << "  multiplicity index " << hypersector.multiplicity_index()
             << std::endl;
        }
      return os.str();
    }

  ////////////////////////////////////////////////////////////////
  // operator storage by hyperblock
  ////////////////////////////////////////////////////////////////


  // typedef for operator hyperblocks
  //
  // EX:
  //   basis::OperatorHyperblocks<double> matrices;  // matrices for all operator planes of all hypersectors
  //   ...
  //   matrix_element = matrices[hypersector_index][operator_index](bra_index,ket_index);  // operator_index is within given hypersector's operator subspace
  
  template <typename tFloat>
    using OperatorHyperblocks = std::vector<std::vector<OperatorBlock<tFloat>>>;

  ////////////////////////////////////////////////////////////////
  // zero operator
  ////////////////////////////////////////////////////////////////

  template <typename tHypersectorsType, typename tFloat>
  void SetHyperoperatorToZero(
      const tHypersectorsType& hypersectors,
      basis::OperatorHyperblocks<tFloat>& matrices
    )
    // Set operator hyperblocks to zero.
    //
    // Arguments:
    //   hypersectors (input): the set of hypersectors on which the
    //     operators are defined
    //   matrices (output): matrices to hold hyperblocks
  {

    // clear vector of matrices
    matrices.clear();
    matrices.resize(hypersectors.size());

    // iterate over hypersectors
    for (int hypersector_index = 0; hypersector_index < hypersectors.size(); ++hypersector_index)
      {

        // extract hypersector
	const typename tHypersectorsType::HypersectorType& hypersector = hypersectors.GetHypersector(hypersector_index);

        // extract hypersector subspaces
	const typename tHypersectorsType::SubspaceType& bra_subspace = hypersector.bra_subspace();
	const typename tHypersectorsType::SubspaceType& ket_subspace = hypersector.ket_subspace();
	const typename tHypersectorsType::OperatorSubspaceType& operator_subspace = hypersector.operator_subspace();

        // generate matrices for hypersector (by operator)
        matrices[hypersector_index].resize(operator_subspace.size());
        for (int operator_index = 0; operator_index < operator_subspace.size(); ++operator_index)
          matrices[hypersector_index][operator_index]
            = basis::OperatorBlock<tFloat>::Zero(bra_subspace.size(),ket_subspace.size());
      }
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
