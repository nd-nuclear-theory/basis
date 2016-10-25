/// @file
/****************************************************************
  indexing.h

  Defines template base classes for indexing quantum mechanical states
  arranged into subspaces (defined by good symmetry quantum numbers).

  General scheme:

    The foundational class is the "subspace" type, representing a
    subspace with good quantum numbers (e.g., LSJT).  This class sets
    up and stores indexing of state quantum numbers within that
    subspace.

    Access to this indexing is via instantiating a "state".  A state
    has to refer to a "subspace" in order to know its quantum numbers.
    Its definition therefore depends upon the "subspace" type.  The
    actual information content of a state consists of a pointer to the
    subspace and an index representing the state within that subspace.

    Going in the other direction in terms of the hierarchy of
    mathematical structures, a "space" is a collection "subspaces".
    Its definition therefore *also* depends upon the "subspace" type.

    Finally, in order to keep track of the matrix elements of an
    operator (of good JP quantum numbers), it is necessary to
    enumerate the allowed sectors, i.e., the pairs of "subspaces"
    within the "space" which are connected under the angular momentum
    and parity selection rules.  These sectors are indexed by a
    "sector enumeration".  Its definition depends upon the "space"
    type, and therefore indirectly on the "subspace" type.

    So notice that the hiererchy of dependencies of template
    definitions

        state <- subspace -> space -> sectors

    does not quite parallel the mathematical hiererchy

        state -> subspace -> space -> sectors.

  Possible extensions:

    Arguably, each subspace and space should have defined a
    TruncationLabelsType, with a corresponding data member and
    accessors.  The constructor would then record the truncation
    scheme for posterity and later reference.

    For standardized debugging output, add a String conversion to
    interface of each type (state, subspace, ...), then use this in
    standardized Print functions.

  Compilation directives:

    INDEXING_HASH: In lookup tables, replace std::map with
      std::unordered_map, using boost hash extensions.

  Library dependences:

    boost -- only if INDEXING_HASH enabled

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  11/21/15 (mac): Created (indexing_base.h), abstracted from code in indexing_lsjt.
  3/5/16 (mac): Update header comment.
  3/9/16 (mac): Relax protection of BaseSpace indexing containers
    from private to protected to support sp3rlib.
  3/11/16 (mac): Add subspace size() method and space
    TotalDimension() and ContainsSubspace() methods.
  5/3/16 (mac): Add boost-style hashing option ad hoc on subspace labels.
  6/8/16 (mac): Extract to new basis module:
    - Rename file from indexing_base.h to indexing.h.
    - Change namespace from shell to basis.
    - Systematize switching between map and unordered_map via INDEXING_HASH
      directive.
    - Remove deprecated subspace method Dimension().
  6/16/16 (mac): Fix missing typename qualifier.
  6/24/16 (mac): Add direction enum class.
  6/30/16 (mac):
    - Flesh out documentation of member functions.
    - Add utility member functions ContainsState, ContainsSector,
      IsDiagonal.
  7/5/16 (mac): Add pass-through typedef StateLabelsType to BaseState.
  7/13/16 (mac): Clean up accessor naming: Add labels() accessor to
    subspace and state, deprecating subspace GetSubspaceLabels() and
    state GetStateLabels().  Add subspace() accessor to state, deprecating
    Subspace().
  7/25/16 (mac): Add utility member function IsUpperTriangle.
  10/25/16 (mac): Implement return of flag value for failed lookup.

****************************************************************/

#ifndef INDEXING_H_
#define INDEXING_H_

#include <cassert>
#include <tuple>
#include <vector>

#ifdef INDEXING_HASH
#include <unordered_map>
#include "boost/functional/hash.hpp"
#else
#include <map>
#endif

namespace basis {

  static const int kNone = -1;
  // Flag value for missing target in index lookups.

  ////////////////////////////////////////////////////////////////
  // generic subspace
  ////////////////////////////////////////////////////////////////

  // BaseSubspace -- holds indexing of states within a symmetry
  // subspace
  //
  // The derived class is expected to set up a constructor and
  // friendlier accessors for the individual labels.
  //
  // Template arguments:
  //   tSubspaceLabelsType : tuple for subspace labels, e.g., std::tuple<int,int,int,int,int>
  //   tStateLabelsType : tuple for state labels, e.g., std::tuple<int>
  //
  // Note: Even if only a single integer label is needed, we must use
  // tuple<int> (as opposed to plain int) to make the two forms of the
  // state constructor syntactically distinct.

  template <typename tSubspaceLabelsType, typename tStateLabelsType>
    class BaseSubspace
  {

    public:

    ////////////////////////////////////////////////////////////////
    //  common type definitions
    ////////////////////////////////////////////////////////////////

    typedef tSubspaceLabelsType SubspaceLabelsType;
    typedef tStateLabelsType StateLabelsType;

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    // default constructor
    //   Implicitly invoked by derived class.
    BaseSubspace() : dimension_(0) {}

    ////////////////////////////////////////////////////////////////
    // retrieval
    ////////////////////////////////////////////////////////////////

    const SubspaceLabelsType& labels() const
    // Return the labels of the subspace itself.
    {
      return labels_;
    }

    const SubspaceLabelsType& GetSubspaceLabels() const
    // Return the labels of the subspace itself.
    //
    // DEPRECATED in favor of labels().
    //
    // To fix: grep GetSubspaceLabels -R --include="*.cpp"
    {
      return labels_;
    }

    const StateLabelsType& GetStateLabels(int index) const
    // Retrieve the labels of a state within the subspace, given its
    // index within the subspace.
    //
    // Note: It is not normally expected that this member function
    // should be used directly.  The idea is instead to work with the
    // *state* type (derived from BaseState) associated with this
    // subspace.  One would rather instantiate a *state*, identified
    // by this index, then query the state for various specific
    // labels, as they are needed, using accessors provided by the
    // state.
    {
      return state_table_[index];
    }

    bool ContainsState(const StateLabelsType& state_labels) const
    // Given the labels for a state, returns whether or not the state
    // is found within the subspace.
    {
      return lookup_.count(state_labels);
    };

    int LookUpStateIndex(const StateLabelsType& state_labels) const
    // Given the labels for a state, look up its index within the
    // subspace.
    //
    // If no such labels are found, basis::kNone is returned.
    {
      // PREVIOUSLY: trap failed lookup with assert for easier debugging
      // assert(ContainsState(state_labels));
      // return lookup_.at(state_labels);

      auto pos = lookup_.find(state_labels);
      if (pos==lookup_.end())
        return kNone;
      else
        return pos->second;
    };

    ////////////////////////////////////////////////////////////////
    // size retrieval
    ////////////////////////////////////////////////////////////////

    int size() const
    // Return the size of the subspace.
    {
      return dimension_;
    }

    protected:

    ////////////////////////////////////////////////////////////////
    // state label push (for initial construction)
    ////////////////////////////////////////////////////////////////

    void PushStateLabels(const StateLabelsType& state_labels)
    // Create indexing information (in both directions, index <->
    // labels) for a state.
    {
      lookup_[state_labels] = dimension_; // index for lookup
      state_table_.push_back(state_labels);  // save state
      dimension_++;
    };

    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

    // subspace properties
    SubspaceLabelsType labels_;
    int dimension_;

    // state labels (accessible by index)
    std::vector<StateLabelsType> state_table_;

    // state index lookup by labels
#ifdef INDEXING_HASH
    std::unordered_map<StateLabelsType,int,boost::hash<StateLabelsType>> lookup_;
#else
    std::map<StateLabelsType,int> lookup_;
#endif

  };

  ////////////////////////////////////////////////////////////////
  // generic state realized within subspace
  ////////////////////////////////////////////////////////////////

  // BaseState -- realization of a state withinin a given subspace
  //
  // The derived class is expected to set up a constructor and
  // friendlier accessors for the individual labels.
  //
  //
  // The space (and the indexing it provides) is *not* copied into the
  // state but rather stored by pointer reference.  It should
  // therefore exist for the lifetime of the state object.
  //
  // Template arguments:
  //   tSubspaceType : subspace type in which this state lives

  template <typename tSubspaceType>
    class BaseState
    {

      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      typedef tSubspaceType SubspaceType;
      typedef typename SubspaceType::StateLabelsType StateLabelsType;

      ////////////////////////////////////////////////////////////////
      // general constructors
      ////////////////////////////////////////////////////////////////

      // default constructor -- disabled
      BaseState();

      // copy constructor -- synthesized

      // constructors

      BaseState(const SubspaceType& subspace, int index)
        // Construct state, given index index within subspace.
        : subspace_ptr_(&subspace), index_(index)
      {
        assert(ValidIndex());
      }

      BaseState(const SubspaceType& subspace, const StateLabelsType& state_labels)
        // Construct state, by reverse lookup on labels within subspace.
        {

          // Debugging: Delegation to BaseState(subspace,index)
          // fails. Argument index is saved in index_ as far as the
          // subordinate BaseState(...,index) call is concerned, but
          // after return to present calling constructor, index_
          // contains garbage.  (Why???)
          //
          // BaseState(subspace,index);

          subspace_ptr_ = &subspace;
          index_ = subspace.LookUpStateIndex(state_labels);
          assert(index_!=basis::kNone);
        }

      ////////////////////////////////////////////////////////////////
      // retrieval
      ////////////////////////////////////////////////////////////////

      const SubspaceType& subspace() const
      // Return reference to subspace in which this state lies.
      {return *subspace_ptr_;}

      const SubspaceType& Subspace() const
      // Return reference to subspace in which this state lies.
      //
      // DEPRECATED in favor of subspace().
      {return *subspace_ptr_;}

      const StateLabelsType& labels() const
      // Return labels of this state.
      //
      // Note: It is not normally expected that this member function
      // should be used directly (see related comment for
      // BaseSubspace::GetStateLabels).  Rather, derived types will
      // provide accessors for convenient access to individual labels.
      {
        return Subspace().GetStateLabels(index());
      }

      const StateLabelsType& GetStateLabels() const
      // Return labels of this state.
      //
      // DEPRECATED in favor of labels().
      //
      // Note: It is not normally expected that this member function
      // should be used directly (see related comment for
      // BaseSubspace::GetStateLabels).  Rather, derived types will
      // provide accessors for convenient access to individual labels.
      {
        return Subspace().GetStateLabels(index());
      }

      int index() const
      // Retrieve integer index of state within subspace.
      {return index_;}

      ////////////////////////////////////////////////////////////////
      // generic iteration support -- DISABLED
      ////////////////////////////////////////////////////////////////

      // currently DISABLED, pending decision about whether or not
      // this is a good thing
      //
      // Example:
      //
      //   for (RelativeStateLSJT state(space); state.ValidIndex(); ++state)
      //     {
      //           std::cout << state.index() << " " << state.N() << std::endl;
      //     };
      //
      // This usage requires ValidIndex() to be made public.

      // BaseState(const SubspaceType& subspace)
      // // Construct state, defaulting to 0th state in space.
      // // Meant for use with iterator-style iteration over states.
      //  {
      //        space_ptr = subspace;
      //        index_ = 0;
      //  }

      // BaseState& operator ++();
      // // Provide prefix increment operator.
      // //
      // // Meant for use with iterator-style iteration over states.
      // {
      //         ++index_;
      //         return *this;
      // }

      ////////////////////////////////////////////////////////////////
      // validation
      ////////////////////////////////////////////////////////////////

      private:

      int ValidIndex() const
      // Verify whether or not state indexing lies within allowed
      // dimension.
      //
      // For use on construction (or possibly with iteration).
      {
        return index() < Subspace().size();
      }

      private:

      ////////////////////////////////////////////////////////////////
      // private storage
      ////////////////////////////////////////////////////////////////

      const SubspaceType* subspace_ptr_;  // subspace in which state lies
      int index_;   // 0-based index within space

    };



  ////////////////////////////////////////////////////////////////
  // generic space
  ////////////////////////////////////////////////////////////////

  // BaseSpace -- container to hold subspaces, with reverse lookup by
  // subspace labels
  //
  // Template arguments:
  //   tSubspaceType (typename) : type for subspace

  template <typename tSubspaceType>
    class BaseSpace
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      typedef tSubspaceType SubspaceType;

      ////////////////////////////////////////////////////////////////
      // subspace lookup and retrieval
      ////////////////////////////////////////////////////////////////

      bool ContainsSubspace(const typename SubspaceType::SubspaceLabelsType& subspace_labels) const
      // Given the labels for a subspace, returns whether or not the
      // subspace is found within the space.
      {
        return lookup_.count(subspace_labels);
      }

      int LookUpSubspaceIndex(const typename SubspaceType::SubspaceLabelsType& subspace_labels) const
      // Given the labels for a subspace, look up its index within the
      // space.
      //
      // If no such labels are found, basis::kNone is returned.
      {

        // PREVIOUSLY: trap failed lookup with assert for easier debugging
        // assert(ContainsSubspace(subspace_labels));
        // return lookup_.at(subspace_labels);

        auto pos = lookup_.find(subspace_labels);
        if (pos==lookup_.end())
          return kNone;
        else
          return pos->second;
      };

      const SubspaceType& LookUpSubspace(const typename SubspaceType::SubspaceLabelsType& subspace_labels) const
      // Given the labels for a subspace, retrieve a reference to the
      // subspace.
      //
      // If no such labels are found, an exception will result
      // (enforced by LookUpSubspaceIndex).
      {

        int subspace_index = LookUpSubspaceIndex(subspace_labels);
        assert(subspace_index!=kNone);
        return subspaces_[subspace_index];
      };

      const SubspaceType& GetSubspace(int i) const
      // Given the index for a subspace, return a reference to the
      // subspace.
      {
        return subspaces_[i];
      };

      ////////////////////////////////////////////////////////////////
      // size retrieval
      ////////////////////////////////////////////////////////////////

      int size() const
      // Return the number of subspaces within the space.
      {
        return subspaces_.size();
      };

      int TotalDimension() const
      // Return the total dimension of all subspaces within the space.
      {
        int dimension = 0;
        for (int i=0; i<size(); ++i)
          dimension += GetSubspace(i).size();
        return dimension;
      }

      protected:

      ////////////////////////////////////////////////////////////////
      // subspace push (for initial construction)
      ////////////////////////////////////////////////////////////////

      void PushSubspace(const SubspaceType& subspace)
      // Create indexing information (in both directions, index <->
      // labels) for a subspace.
      {
        lookup_[subspace.labels()] = subspaces_.size(); // index for lookup
        subspaces_.push_back(subspace);  // save space
      };

      ////////////////////////////////////////////////////////////////
      // internal storage
      ////////////////////////////////////////////////////////////////

      // subspaces (accessible by index)
      std::vector<SubspaceType> subspaces_;

      // subspace index lookup by labels
#ifdef INDEXING_HASH
      std::unordered_map<typename SubspaceType::SubspaceLabelsType,int,boost::hash<typename SubspaceType::SubspaceLabelsType>> lookup_;
#else
      std::map<typename SubspaceType::SubspaceLabelsType,int> lookup_;
#endif
    };

  ////////////////////////////////////////////////////////////////
  // sector indexing
  ////////////////////////////////////////////////////////////////

  // Note: A "sector" is a pair of subspaces (which thus define a
  // "sector" or block in the matrix representation of an operator on
  // the space).
  //
  // The sector may also be labeled with a multiplicity index.  The
  // concept of a multiplicity index on a sector applies when the
  // symmetry group has outer multiplicites, so reduced matrix
  // elements are labeled not only by the bra and ket but also by a
  // multiplicity index.


  // BaseSector -- provide storage of information for a single sector

  template <typename tSubspaceType>
    class BaseSector
    // Store indexing and subspace reference information for a sector.
    //
    // Allows for multiplicity index on sector, for symmetry sectors
    // of groups with outer multiplicities.
    {

      ////////////////////////////////////////////////////////////////
      // typedefs
      ////////////////////////////////////////////////////////////////

      public:
      typedef tSubspaceType SubspaceType;
      typedef std::tuple<int,int,int> KeyType;

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSector(
          int bra_subspace_index, int ket_subspace_index,
          const SubspaceType& bra_subspace, const SubspaceType& ket_subspace,
          int multiplicity_index=1
        )
        : bra_subspace_index_(bra_subspace_index), ket_subspace_index_(ket_subspace_index), multiplicity_index_(multiplicity_index)
      {
        bra_subspace_ptr_ = &bra_subspace;
        ket_subspace_ptr_ = &ket_subspace;
      }

      ////////////////////////////////////////////////////////////////
      // accessors
      ////////////////////////////////////////////////////////////////

      inline KeyType Key() const
      // Return tuple key identifying sector for sorting/lookup
      // purposes.
      {
        return KeyType(bra_subspace_index(),ket_subspace_index(),multiplicity_index());
      }

      int bra_subspace_index() const {return bra_subspace_index_;}
      int ket_subspace_index() const {return ket_subspace_index_;}
      // Return integer index of bra/ket subspace.

      const SubspaceType& bra_subspace() const {return *bra_subspace_ptr_;}
      const SubspaceType& ket_subspace() const {return *ket_subspace_ptr_;}
      // Return reference to bra/ket subspace.

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
      int bra_subspace_index_, ket_subspace_index_;
      const SubspaceType* bra_subspace_ptr_;
      const SubspaceType* ket_subspace_ptr_;
      int multiplicity_index_;
    };

  // SectorDirection -- sector direction specifier
  //
  //   kCanonical : specifies bra_sector_index <= ket_sector_index
  //   kBoth : specifies that both directions are allowed
  //
  // Note: It is up to the derived class constructor to accept a
  // sector direction argument and to implement this symmetry
  // constraint, if so decided.

  enum class SectorDirection {kCanonical,kBoth};

  // BaseSectors -- container to hold a set of sectors with reverse
  // lookup by sector labels
  //
  // Template arguments:
  //   tSpaceType (typename) : type for space

  template <typename tSpaceType>
    class BaseSectors
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      typedef tSpaceType SpaceType;
      typedef typename tSpaceType::SubspaceType SubspaceType;
      typedef BaseSector<SubspaceType> SectorType;

      ////////////////////////////////////////////////////////////////
      // subspace lookup and retrieval
      ////////////////////////////////////////////////////////////////

      const SectorType& GetSector(int i) const
      // Given sector index, return reference to sector itself.
      {
        return sectors_[i];
      };

      bool ContainsSector(int bra_subspace_index, int ket_subspace_index, int multiplicity_index=1) const
      // Given the labels for a sector, returns whether or not the sector
      // is found within the the sector set.
      {
        typename SectorType::KeyType key(bra_subspace_index,ket_subspace_index,multiplicity_index);
        return lookup_.count(key);
      };

      int LookUpSectorIndex(int bra_subspace_index, int ket_subspace_index, int multiplicity_index=1) const
      // Given the labels for a sector, look up its index within the
      // sector set.
      //
      // If no such labels are found, basis::kNone is returned.
      {

        // PREVIOUSLY: trap failed lookup with assert for easier debugging
        // assert(ContainsSector(bra_subspace_index,ket_subspace_index,multiplicity_index));
        // return lookup_.at(key);

        typename SectorType::KeyType key(bra_subspace_index,ket_subspace_index,multiplicity_index);
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
      // Return number of sectors within sector set.
      {
        return sectors_.size();
      };

      protected:

      ////////////////////////////////////////////////////////////////
      // sector push (for initial construction)
      ////////////////////////////////////////////////////////////////

      void PushSector(const SectorType& sector)
      // Create indexing information (in both directions, index <->
      // labels) for a sector.
      {
        lookup_[sector.Key()] = sectors_.size(); // index for lookup
        sectors_.push_back(sector);  // save sector
      };

      ////////////////////////////////////////////////////////////////
      // internal storage
      ////////////////////////////////////////////////////////////////

      // sectors (accessible by index)
      std::vector<SectorType> sectors_;

      // sector index lookup by subspace indices
#ifdef INDEXING_HASH
      std::unordered_map<typename SectorType::KeyType,int,boost::hash<typename SectorType::KeyType>> lookup_;
#else
      std::map<typename SectorType::KeyType,int> lookup_;
#endif

    };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
