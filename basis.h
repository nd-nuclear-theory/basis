/************************************************************//**
  @file basis.h

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

    So notice that the hierarchy of dependencies of template
    definitions

        state <- subspace -> space -> sectors

    does not quite parallel the mathematical hierarchy

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

    BASIS_HASH: In lookup tables, replace std::map with
      std::unordered_map, using boost hash extensions.

  Library dependences:

    boost -- only if BASIS_HASH enabled

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 11/21/15 (mac): Created (indexing_base.h), abstracted from code in indexing_lsjt.
  + 3/5/16 (mac): Update header comment.
  + 3/9/16 (mac): Relax protection of BaseSpace indexing containers
    - from private to protected to support sp3rlib.
  + 3/11/16 (mac): Add subspace size() method and space
    - TotalDimension() and ContainsSubspace() methods.
  + 5/3/16 (mac): Add boost-style hashing option ad hoc on subspace labels.
  + 6/8/16 (mac): Extract to new basis module:
    - Rename file from indexing_base.h to indexing.h.
    - Change namespace from shell to basis.
    - Systematize switching between map and unordered_map via BASIS_HASH
      directive.
    - Remove deprecated subspace method Dimension().
  + 6/16/16 (mac): Fix missing typename qualifier.
  + 6/24/16 (mac): Add direction enum class.
  + 6/30/16 (mac):
    - Flesh out documentation of member functions.
    - Add utility member functions ContainsState, ContainsSector,
      IsDiagonal.
  + 7/5/16 (mac): Add pass-through typedef StateLabelsType to BaseState.
  + 7/13/16 (mac): Clean up accessor naming: Add labels() accessor to
    subspace and state, deprecating subspace GetSubspaceLabels() and
    state GetStateLabels().  Add subspace() accessor to state, deprecating
    Subspace().
  + 7/25/16 (mac): Add utility member function IsUpperTriangle.
  + 10/25/16 (mac):
    - Implement return of flag value for failed lookup.
    - Provide sector lookup by sector key and refactor ContainsSector and
      LookUpSectorIndex.
  + 11/22/16 (mac):
    - Fix wrong return in broken-out ContainsSector.
    - Reexpand (remultiply?) ContainsSector and LookUpSectorIndex.
  + 2/17/17 (mac): Add generic DebugStr for BaseSectors.
  + 4/4/17 (mac):
    - Rename from indexing.h to basis.h.
    - Rename compilation flag INDEXING_HASH to BASIS_HASH.
  + 6/6/7 (mac): Disable deprecated member functions by default
    (BASIS_ALLOW_DEPRECATED).
  + 6/11/17 (mac): Rename TotalDimension to Dimension.
  + 08/11/17 (pjf): Emit warnings if deprecated member functions are used.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.
  + 05/27/19 (pjf): Modify BaseSpace and BaseSectors to fix fragility problems:
    - Move BaseSpace subspaces_ and lookup_ into shared_ptrs; copying a
      space is now a lightweight operation.
    - Store copy of bra and ket spaces inside BaseSectors.
    - Construct BaseSector at access time, dereferencing from local copy of
      space.
    - Add BaseSpace::EmplaceSubspace() for in-place addition of subspaces.
  + 06/24/21 (pjf): Modify to allow unlimited nesting of spaces:
    - Extract labels to BaseLabeling, with an explicit (empty) specialization
      for labels of type void (i.e. no labels).
    - Add new template parameter to BaseSpace to define a LabelsType for spaces;
      this allows instantiating a space of spaces.
    - Add dimension() accessor to BaseSubspace and BaseSpace.
    - Store offsets to subspaces and total dimension directly in BaseSpace.
  + 06/25/21 (pjf):
    - Add accessor for subspace offset in BaseSpace.
    - Fix initialization of dimension_ in BaseSpace.
    - Add BaseDegenerateSpace as friend of BaseSpace.
  + 06/28/21 (pjf): Fix friend declaration in BaseSpace.
  + 07/02/21 (pjf):
    - Replace `typedef`s with alias declarations (`using`).
    - Add constructors for BaseSubspace and BaseSpace which accept labels.
  + 07/03/21 (pjf):
    - Remove BaseLabeling and replace with inheritance from partial
      template specialization.
    - Mark labels_ as private and const to prevent modification
      after construction.
    - Mark constructors as protected to prevent direct initialization.
  + 07/04/21 (pjf):
    - Add tDerivedSubspaceType and tStateType as template arguments to
      BaseSubspace.
    - Add GetState accessors to BaseSubspace.
  + 08/08/21 (pjf):
    - Modify BaseSector and BaseSectors to allow bra and ket space types to
      differ.
    - Add additional template argument to BaseSectors so that custom sector type
      can be used, instead of default instantiation of BaseSector.
  + 08/16/21 (pjf):
    - Switch from storing shared_ptr<vector> to vector<shared_ptr> in BaseSpace.
    - Add tDerivedSpaceType as template argument to BaseSpace for CRTP with
      std::enable_shared_from_this.
    - Replace raw pointers with shared_ptr in BaseSector.
  + 08/17/21 (pjf):
    - Revert to storing instances of sector in BaseSectors rather than just keys.
    - Only copy space into BaseSectors if it is not managed by a shared_ptr.
****************************************************************/

#ifndef BASIS_BASIS_H_
#define BASIS_BASIS_H_

#include <cassert>
#include <cstddef>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

// for DebugStr
#include <iomanip>
#include <sstream>

// emit warnings on deprecated
#include "mcutils/deprecated.h"

#ifdef BASIS_HASH
#include <unordered_map>
#include "boost/functional/hash.hpp"
#else
#include <map>
#endif

namespace basis {

  ////////////////////////////////////////////////////////////////
  // generic subspace mixin
  ////////////////////////////////////////////////////////////////

  /// BaseLabeling -- holds labeling and accessors for subspaces
  ///
  /// Users are not expected to directly inherit from this class; instead,
  /// use BaseSubspace and BaseSpace
  ///
  /// Template arguments:
  ///   tSubspaceLabelsType : tuple for subspace labels, e.g., std::tuple<int,int,int,int,int>
  ///
  /// Note: Even if only a single integer label is needed, we must use
  /// tuple<int> (as opposed to plain int) to make the two forms of the
  /// state constructor syntactically distinct.

  static constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();
  // Flag value for missing target in index lookups.

  ////////////////////////////////////////////////////////////////
  // generic subspace
  ////////////////////////////////////////////////////////////////

  /// BaseSubspace -- holds indexing of states within a symmetry
  /// subspace
  ///
  /// The derived class is expected to set up a constructor and
  /// friendlier accessors for the individual labels.
  ///
  /// The derived class passes itself as the first template argument, as an
  /// example of the Curiously Recurring Template Pattern (CRTP). This allows
  /// BaseSubspace access to the full type information of the derived class.
  ///
  /// Template arguments:
  ///   tDerivedSubspaceType : type of the derived subspace class (for CRTP)
  ///   tSubspaceLabelsType : tuple for subspace labels, e.g., std::tuple<int,int,int,int,int>
  ///   tStateType : type (possibly incomplete) for states
  ///   tStateLabelsType : tuple for state labels, e.g., std::tuple<int>
  ///
  /// Note: Even if only a single integer label is needed, we must use
  /// tuple<int> (as opposed to plain int) to make the two forms of the
  /// state constructor syntactically distinct.

  template <
      typename tDerivedSubspaceType, typename tSubspaceLabelsType,
      typename tStateType, typename tStateLabelsType
    >
  class BaseSubspace
      : public std::enable_shared_from_this<tDerivedSubspaceType>
  {

    public:

    ////////////////////////////////////////////////////////////////
    //  common typedefs
    ////////////////////////////////////////////////////////////////

    using LabelsType = tSubspaceLabelsType;
    using SubspaceLabelsType = tSubspaceLabelsType;
    using StateLabelsType = tStateLabelsType;

    protected:

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    /// default constructor
    //   Implicitly invoked by derived class.
    BaseSubspace() : dimension_(0) {}

    // pass-through constructor accepting labels
    explicit BaseSubspace(const SubspaceLabelsType& labels)
      : dimension_{0}, labels_{labels}
    {}

    public:

    ////////////////////////////////////////////////////////////////
    // retrieval
    ////////////////////////////////////////////////////////////////

    const SubspaceLabelsType& labels() const
    /// Return the labels of the subspace itself.
    {
      return labels_;
    }

#ifdef BASIS_ALLOW_DEPRECATED
    DEPRECATED("use labels() instead")
    const SubspaceLabelsType& GetSubspaceLabels() const
    /// Return the labels of the subspace itself.
    ///
    /// DEPRECATED in favor of labels().
    ///
    /// To fix: grep GetSubspaceLabels -R --include="*.cpp"
    {
      return labels_;
    }
#endif

    const StateLabelsType& GetStateLabels(std::size_t index) const
    /// Retrieve the labels of a state within the subspace, given its
    /// index within the subspace.
    ///
    /// Note: It is not normally expected that this member function
    /// should be used directly.  The idea is instead to work with the
    /// *state* type (derived from BaseState) associated with this
    /// subspace.  One would rather instantiate a *state*, identified
    /// by this index, then query the state for various specific
    /// labels, as they are needed, using accessors provided by the
    /// state.
    {
      return state_table_[index];
    }

    bool ContainsState(const StateLabelsType& state_labels) const
    /// Given the labels for a state, returns whether or not the state
    /// is found within the subspace.
    {
      return lookup_.count(state_labels);
    };

    std::size_t LookUpStateIndex(const StateLabelsType& state_labels) const
    /// Given the labels for a state, look up its index within the
    /// subspace.
    ///
    /// If no such labels are found, basis::kNone is returned.
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

    tStateType GetState(std::size_t index) const;
    /// Given the index for a state, construct the corresponding state object.
    ///
    /// Defined below, outside the definition of the class, so that
    /// instantiation is deferred until tStateType is a complete type.

    tStateType GetState(const StateLabelsType& state_labels) const;
    /// Given the labels for a state, construct the corresponding state object.
    ///
    /// Defined below, outside the definition of the class, so that
    /// instantiation is deferred until tStateType is a complete type.

    ////////////////////////////////////////////////////////////////
    // size retrieval
    ////////////////////////////////////////////////////////////////

    std::size_t size() const
    /// Return the size of the subspace.
    {
      return dimension_;
    }

    std::size_t dimension() const
    /// Return the size of the subspace.
    {
      return dimension_;
    }

    protected:

    ////////////////////////////////////////////////////////////////
    // state label push (for initial construction)
    ////////////////////////////////////////////////////////////////

    void PushStateLabels(const StateLabelsType& state_labels)
    /// Create indexing information (in both directions, index <->
    /// labels) for a state.
    {
      lookup_[state_labels] = dimension_; // index for lookup
      state_table_.push_back(state_labels);  // save state
      dimension_++;
    };

    private:

    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

    // subspace properties
    const SubspaceLabelsType labels_;
    std::size_t dimension_;

    // state labels (accessible by index)
    std::vector<StateLabelsType> state_table_;

    // state index lookup by labels
#ifdef BASIS_HASH
    std::unordered_map<StateLabelsType,std::size_t,boost::hash<StateLabelsType>> lookup_;
#else
    std::map<StateLabelsType,std::size_t> lookup_;
#endif

  };

  // implementations of BaseSubspace::GetState; we do this outside the class
  // declaration to delay instantiation until tDerivedSubspaceType and
  // tStateType are complete types

  template <
      typename tDerivedSubspaceType, typename tSubspaceLabelsType,
      typename tStateType, typename tStateLabelsType
    >
  tStateType BaseSubspace<tDerivedSubspaceType,tSubspaceLabelsType,tStateType,tStateLabelsType>::GetState(
      std::size_t index
    ) const
  {
    return tStateType{*static_cast<const tDerivedSubspaceType*>(this), index};
  }

  template <
      typename tDerivedSubspaceType, typename tSubspaceLabelsType,
      typename tStateType, typename tStateLabelsType
    >
  tStateType BaseSubspace<tDerivedSubspaceType,tSubspaceLabelsType,tStateType,tStateLabelsType>::GetState(
      const tStateLabelsType& state_labels
    ) const
  {
    return tStateType{*static_cast<const tDerivedSubspaceType*>(this), state_labels};
  }

  ////////////////////////////////////////////////////////////////
  // generic state realized within subspace
  ////////////////////////////////////////////////////////////////

  /// BaseState -- realization of a state withinin a given subspace
  ///
  /// The derived class is expected to set up a constructor and
  /// friendlier accessors for the individual labels.
  ///
  /// The space (and the indexing it provides) is *not* copied into the
  /// state but rather stored by pointer reference.  It should
  /// therefore exist for the lifetime of the state object.
  ///
  /// Template arguments:
  ///   tSubspaceType : subspace type in which this state lives

  template <typename tSubspaceType>
    class BaseState
    {

      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using SubspaceType = tSubspaceType;
      using StateLabelsType = typename SubspaceType::StateLabelsType;

      protected:

      ////////////////////////////////////////////////////////////////
      // general constructors
      ////////////////////////////////////////////////////////////////

      /// default constructor -- disabled
      BaseState();

      // copy constructor -- synthesized

      // constructors

      BaseState(const SubspaceType& subspace, std::size_t index)
        /// Construct state, given index within subspace.
        : subspace_ptr_(&subspace), index_(index)
      {
        assert(ValidIndex());
      }

      BaseState(const SubspaceType& subspace, const StateLabelsType& state_labels)
        /// Construct state, by reverse lookup on labels within subspace.
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

      public:

      ////////////////////////////////////////////////////////////////
      // retrieval
      ////////////////////////////////////////////////////////////////

      const SubspaceType& subspace() const
      /// Return reference to subspace in which this state lies.
      {return *subspace_ptr_;}

#ifdef BASIS_ALLOW_DEPRECATED
      DEPRECATED("use subspace() instead")
      const SubspaceType& Subspace() const
      /// Return reference to subspace in which this state lies.
      ///
      /// DEPRECATED in favor of subspace().
      {return *subspace_ptr_;}
#endif

      const StateLabelsType& labels() const
      /// Return labels of this state.
      ///
      /// Note: It is not normally expected that this member function
      /// should be used directly (see related comment for
      /// BaseSubspace::GetStateLabels).  Rather, derived types will
      /// provide accessors for convenient access to individual labels.
      {
        return subspace().GetStateLabels(index());
      }

#ifdef BASIS_ALLOW_DEPRECATED
      DEPRECATED("use labels() instead")
      const StateLabelsType& GetStateLabels() const
      /// Return labels of this state.
      ///
      /// DEPRECATED in favor of labels().
      ///
      /// Note: It is not normally expected that this member function
      /// should be used directly (see related comment for
      /// BaseSubspace::GetStateLabels).  Rather, derived types will
      /// provide accessors for convenient access to individual labels.
      {
        return Subspace().GetStateLabels(index());
      }
#endif

      std::size_t index() const
      /// Retrieve integer index of state within subspace.
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

      std::size_t ValidIndex() const
      /// Verify whether or not state indexing lies within allowed
      /// dimension.
      ///
      /// For use on construction (or possibly with iteration).
      {
        return index() < subspace().size();
      }

      private:

      ////////////////////////////////////////////////////////////////
      // private storage
      ////////////////////////////////////////////////////////////////

      const SubspaceType* subspace_ptr_;  ///< subspace in which state lies
      std::size_t index_;   ///< 0-based index within space

    };



  ////////////////////////////////////////////////////////////////
  // generic space
  ////////////////////////////////////////////////////////////////

  /// BaseSpace -- container to hold subspaces, with reverse lookup by
  /// subspace labels
  ///
  /// Template arguments:
  ///   tSubspaceType (typename) : type for subspace
  ///   tSpaceLabelsType (typename, optional) : type for space labels

  // declare BaseSpace template
  template <typename tDerivedSpaceType, typename tSubspaceType, typename tSpaceLabelsType = void>
    class BaseSpace;

  // specialize class for case with no labels
  template <typename tDerivedSpaceType, typename tSubspaceType>
    class BaseSpace<tDerivedSpaceType, tSubspaceType, void>
        : public std::enable_shared_from_this<tDerivedSpaceType>
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using SubspaceType = tSubspaceType;

      protected:

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSpace()
        : dimension_{0}, subspace_ptrs_{}, lookup_{}
      {}

      public:

      ////////////////////////////////////////////////////////////////
      // subspace lookup and retrieval
      ////////////////////////////////////////////////////////////////

      bool ContainsSubspace(const typename SubspaceType::LabelsType& subspace_labels) const
      /// Given the labels for a subspace, returns whether or not the
      /// subspace is found within the space.
      {
        return lookup_.count(subspace_labels);
      }

      std::size_t LookUpSubspaceIndex(
          const typename SubspaceType::LabelsType& subspace_labels
        ) const
      /// Given the labels for a subspace, look up its index within the
      /// space.
      ///
      /// If no such labels are found, basis::kNone is returned.
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

      const SubspaceType& LookUpSubspace(
          const typename SubspaceType::LabelsType& subspace_labels
        ) const
      /// Given the labels for a subspace, retrieve a reference to the
      /// subspace.
      ///
      /// If no such labels are found, an exception will result
      /// (enforced by LookUpSubspaceIndex).
      {

        std::size_t subspace_index = LookUpSubspaceIndex(subspace_labels);
        assert(subspace_index!=kNone);
        return *(subspace_ptrs_.at(subspace_index));
      };

      const SubspaceType& GetSubspace(std::size_t i) const
      /// Given the index for a subspace, return a reference to the
      /// subspace.
      {
        return *(subspace_ptrs_.at(i));
      };

      std::size_t GetSubspaceOffset(std::size_t i) const
      /// Given the index for a subspace, return the offset of the subspace
      /// in the full space indexing.
      {
        return subspace_offsets_.at(i);
      }

      ////////////////////////////////////////////////////////////////
      // size retrieval
      ////////////////////////////////////////////////////////////////

      std::size_t size() const
      /// Return the number of subspaces within the space.
      {
        return subspace_ptrs_.size();
      };

      std::size_t dimension() const
      // Return the total dimension of all subspaces within the space.
      {
        return dimension_;
      }

#ifdef BASIS_ALLOW_DEPRECATED
      DEPRECATED("use dimension() instead")
      std::size_t Dimension() const
      /// Return the total dimension of all subspaces within the space.
      {
        std::size_t dimension = 0;
        for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
          dimension += GetSubspace(subspace_index).size();
        return dimension;
      }
#endif

      protected:

      ////////////////////////////////////////////////////////////////
      // subspace push (for initial construction)
      ////////////////////////////////////////////////////////////////

      template<typename T, typename std::enable_if_t<std::is_same_v<std::decay_t<T>, SubspaceType>>* = nullptr>
      void PushSubspace(T&& subspace)
      /// Create indexing information (in both directions, index <->
      /// labels) for a subspace.
      {
        const std::size_t index = subspace_ptrs_.size();  // index for lookup
        subspace_offsets_.push_back(dimension_);
        subspace_ptrs_.push_back(
            std::make_shared<const SubspaceType>(std::forward<T>(subspace))
          );  // save subspace
        lookup_[subspace.labels()] = index;  // index for lookup
        dimension_ += subspace.dimension();
      };

      template <class... Args>
      void EmplaceSubspace(Args&&... args)
      /// Create indexing information (in both directions, index <->
      /// labels) for a subspace.
      {
        const std::size_t index = subspace_ptrs_.size();  // index for lookup
        subspace_offsets_.push_back(dimension_);
        subspace_ptrs_.push_back(
            std::make_shared<const SubspaceType>(std::forward<Args>(args)...)
          );
        const SubspaceType& subspace = *(subspace_ptrs_.back());
        lookup_[subspace.labels()] = index;
        dimension_ += subspace.dimension();
      }

      private:

      // allow BaseDegenerateSpace to access private members to override
      // PushSubspace and EmplaceSubspace
      template<typename, typename, typename> friend class BaseDegenerateSpace;

      ////////////////////////////////////////////////////////////////
      // internal storage
      ////////////////////////////////////////////////////////////////

      /// subspaces (accessible by index)
      std::vector<std::shared_ptr<const SubspaceType>> subspace_ptrs_;
      std::vector<std::size_t> subspace_offsets_;
      std::size_t dimension_;

      // subspace index lookup by labels
#ifdef BASIS_HASH
      std::unordered_map<
          typename SubspaceType::LabelsType,std::size_t,
          boost::hash<typename SubspaceType::LabelsType>
        > lookup_;
#else
      std::map<typename SubspaceType::LabelsType,std::size_t> lookup_;
#endif
    };

  // inherit from BaseSpace and add labels
  template <typename tDerivedSpaceType, typename tSubspaceType, typename tSpaceLabelsType>
    class BaseSpace : public BaseSpace<tDerivedSpaceType, tSubspaceType, void>
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using LabelsType = tSpaceLabelsType;
      using SpaceLabelsType = tSpaceLabelsType;

      protected:

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSpace() = default;

      // pass-through constructor accepting labels
      explicit BaseSpace(const SpaceLabelsType& labels)
        : labels_{labels}
      {}

      public:

      ////////////////////////////////////////////////////////////////
      // retrieval
      ////////////////////////////////////////////////////////////////

      const SpaceLabelsType& labels() const
      /// Return the labels of the subspace itself.
      {
        return labels_;
      }

      private:

      ////////////////////////////////////////////////////////////////
      // internal storage
      ////////////////////////////////////////////////////////////////

      /// space labels
      const SpaceLabelsType labels_;
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


  // declare BaseSectors template
  template <
      typename tBraSubspaceType, typename tKetSubspaceType = tBraSubspaceType,
      bool = std::is_same_v<tBraSubspaceType, tKetSubspaceType>
    >
    class BaseSector;

  /// BaseSector -- provide storage of information for a single sector

  // Here we specialize the template for the case where the bra and ket
  // subspaces are of different types. For the case where they are the same
  // (below), we will inherit from this specialization.
  template <typename tBraSubspaceType, typename tKetSubspaceType>
    class BaseSector<tBraSubspaceType, tKetSubspaceType, false>
    // Store indexing and subspace reference information for a sector.
    //
    // Allows for multiplicity index on sector, for symmetry sectors
    // of groups with outer multiplicities.
    {

      ////////////////////////////////////////////////////////////////
      // typedefs
      ////////////////////////////////////////////////////////////////

      public:
      using BraSubspaceType = tBraSubspaceType;
      using KetSubspaceType = tKetSubspaceType;
      typedef std::tuple<std::size_t,std::size_t,std::size_t> KeyType;

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSector(
          std::size_t bra_subspace_index, std::size_t ket_subspace_index,
          const BraSubspaceType& bra_subspace, const KetSubspaceType& ket_subspace,
          std::size_t multiplicity_index=1
        )
        : bra_subspace_index_(bra_subspace_index),
          ket_subspace_index_(ket_subspace_index),
          multiplicity_index_(multiplicity_index),
          bra_subspace_ptr_(bra_subspace.shared_from_this()),
          ket_subspace_ptr_(ket_subspace.shared_from_this())
      {}

      ////////////////////////////////////////////////////////////////
      // accessors
      ////////////////////////////////////////////////////////////////

      inline KeyType Key() const
      // Return tuple key identifying sector for sorting/lookup
      // purposes.
      {
        return KeyType(bra_subspace_index(),ket_subspace_index(),multiplicity_index());
      }

      std::size_t bra_subspace_index() const {return bra_subspace_index_;}
      std::size_t ket_subspace_index() const {return ket_subspace_index_;}
      // Return integer index of bra/ket subspace.

      const BraSubspaceType& bra_subspace() const {return *bra_subspace_ptr_;}
      const KetSubspaceType& ket_subspace() const {return *ket_subspace_ptr_;}
      // Return reference to bra/ket subspace.

      std::size_t multiplicity_index() const {return multiplicity_index_;}
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
      std::size_t bra_subspace_index_, ket_subspace_index_;
      std::size_t multiplicity_index_;
      const std::shared_ptr<const BraSubspaceType> bra_subspace_ptr_;
      const std::shared_ptr<const KetSubspaceType> ket_subspace_ptr_;
    };

  // Here we specialize to the case where the bra and ket subspaces have the
  // same type. We inherit all the functionality from the case where they
  // differ, and add an additional type alias.
  template<typename tSubspaceType>
    class BaseSector<tSubspaceType, tSubspaceType, true>
      : public BaseSector<tSubspaceType, tSubspaceType, false>
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using SubspaceType = tSubspaceType;

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      inline BaseSector(
          std::size_t bra_subspace_index, std::size_t ket_subspace_index,
          const SubspaceType& bra_subspace, const SubspaceType& ket_subspace,
          std::size_t multiplicity_index=1
        )
        : BaseSector<tSubspaceType, tSubspaceType, false>{
            bra_subspace_index, ket_subspace_index,
            bra_subspace,ket_subspace, multiplicity_index
        }
        {}
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
  //   tBraSpaceType (typename) : type for bra space
  //   tKetSpaceType (typename) : type for ket space
  //   tSectorType (typename) : type for sector
  //   (bool): true if tBraSpaceType and tKetSpaceType are the same

  // declare BaseSectors template
  template <
      typename tBraSpaceType, typename tKetSpaceType = tBraSpaceType,
      typename tSectorType
        = BaseSector<typename tBraSpaceType::SubspaceType, typename tKetSpaceType::SubspaceType>,
      bool = std::is_same_v<tBraSpaceType, tKetSpaceType>
    >
    class BaseSectors;

  // Here we specialize the template for the case where the bra and ket spaces
  // are of different types. For the case where they are the same (below), we
  // will inherit from this specialization.
  template<typename tBraSpaceType, typename tKetSpaceType, typename tSectorType>
    class BaseSectors<tBraSpaceType, tKetSpaceType, tSectorType, false>
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using BraSpaceType = tBraSpaceType;
      using KetSpaceType = tKetSpaceType;

      using BraSubspaceType = typename tBraSpaceType::SubspaceType;
      using KetSubspaceType = typename tKetSpaceType::SubspaceType;
      using SectorType = tSectorType;


      protected:

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSectors() = default;
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      template<
          typename T, typename U,
          typename std::enable_if_t<std::is_same_v<std::decay_t<T>, BraSpaceType>>* = nullptr,
          typename std::enable_if_t<std::is_same_v<std::decay_t<U>, KetSpaceType>>* = nullptr
        >
      BaseSectors(T&& bra_space, U&& ket_space)
        : bra_space_ptr_(bra_space.weak_from_this().lock()),
          ket_space_ptr_(ket_space.weak_from_this().lock())
      {
        // if bra_space or ket_space is not managed by a shared_ptr, make copy
        // or move it to a locally-created shared_ptr.
        //
        // Note: this is probably safe, since the spaces themselves store
        // vectors of shared_ptrs to subspaces, so a copy only involves copying
        // those vectors and maps
        if (!bra_space_ptr_)
          bra_space_ptr_
            = std::make_shared<const BraSpaceType>(std::forward<T>(bra_space));
        if (!ket_space_ptr_)
        {
          // if bra and ket spaces are identical, we can share the new copy we
          // just made for the bra_space_ptr_
          //
          // note: this also ensures that if bra_space and ket_space are the
          // same object, then if we just moved from bra_space we don't try to
          // move from it again (which would be UB)
          if (&ket_space == &bra_space)
            ket_space_ptr_ = bra_space_ptr_;
          else
            ket_space_ptr_
              = std::make_shared<const KetSpaceType>(std::forward<U>(ket_space));
        }
      }

      public:

      ////////////////////////////////////////////////////////////////
      // space accessors
      ////////////////////////////////////////////////////////////////

      const BraSpaceType& bra_space() const
      {
        return *bra_space_ptr_;
      }

      const KetSpaceType& ket_space() const
      {
        return *ket_space_ptr_;
      }

      ////////////////////////////////////////////////////////////////
      // sector lookup and retrieval
      ////////////////////////////////////////////////////////////////

      const SectorType& GetSector(std::size_t sector_index) const
      // Given sector index, return reference to sector itself.
      {
        return sectors_.at(sector_index);
      };

      bool ContainsSector(typename SectorType::KeyType key) const
      // Given the labels for a sector, returns whether or not the sector
      // is found within the the sector set.
      {
        return lookup_.count(key);
      };

      bool ContainsSector(
          std::size_t bra_subspace_index,
          std::size_t ket_subspace_index,
          std::size_t multiplicity_index=1
        ) const
      // Given the labels for a sector, returns whether or not the sector
      // is found within the the sector set.
      {
        return ContainsSector(
            typename SectorType::KeyType{bra_subspace_index,ket_subspace_index,multiplicity_index}
          );
      };

      std::size_t LookUpSectorIndex(const typename SectorType::KeyType& key) const
      // Given the key for a sector, look up its index within the
      // sector set.
      //
      // If no such labels are found, basis::kNone is returned.
      {
        auto pos = lookup_.find(key);
        if (pos==lookup_.end())
          return kNone;
        else
          return pos->second;
      };

      std::size_t LookUpSectorIndex(
          std::size_t bra_subspace_index,
          std::size_t ket_subspace_index,
          std::size_t multiplicity_index=1
        ) const
      // Given the labels for a sector, look up its index within the
      // sector set.
      //
      // If no such labels are found, basis::kNone is returned.
      {
        return LookUpSectorIndex(
            typename SectorType::KeyType{bra_subspace_index,ket_subspace_index,multiplicity_index}
          );
      }

      ////////////////////////////////////////////////////////////////
      // size retrieval
      ////////////////////////////////////////////////////////////////

      std::size_t size() const
      // Return number of sectors within sector set.
      {
        return sectors_.size();
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
      // sector push (for initial construction)
      ////////////////////////////////////////////////////////////////

      template<
          typename T,
          typename std::enable_if_t<std::is_same_v<std::decay_t<T>,SectorType>>* = nullptr
        >
      void PushSector(T&& sector)
      // Create indexing information (in both directions, index <->
      // labels) for a sector.
      {
        const std::size_t index = sectors_.size();
        sectors_.push_back(std::forward<T>(sector));  // save sector
        lookup_[sector.Key()] = index; // index for lookup
      };

      template<typename... Args>
      void EmplaceSector(Args&&... args)
      // Create indexing information (in both directions, index <->
      // labels) for a sector.
      {
        const std::size_t index = sectors_.size();
        sectors_.emplace_back(std::forward<Args>(args)...);  // save sector
        const SectorType& sector = sectors_.back();
        lookup_[sector.Key()] = index; // index for lookup
      };

      inline void PushSector(
          std::size_t bra_subspace_index,
          std::size_t ket_subspace_index,
          std::size_t multiplicity_index=1
        )
      // Create indexing information (in both directions, index <->
      // labels) for a sector, given indices.
      {
        EmplaceSector(
            bra_subspace_index, ket_subspace_index,
            bra_space_ptr_->GetSubspace(bra_subspace_index),
            ket_space_ptr_->GetSubspace(ket_subspace_index),
            multiplicity_index
          );  // save sector
      }

      inline void PushSector(const typename SectorType::KeyType& key)
      // Create indexing information (in both directions, index <->
      // labels) for a sector, given key.
      {
        const auto& [bra_subspace_index,ket_subspace_index,multiplicity_index] = key;
        PushSector(bra_subspace_index, ket_subspace_index, multiplicity_index);
      };

      private:
      ////////////////////////////////////////////////////////////////
      // internal storage
      ////////////////////////////////////////////////////////////////

      // spaces
      std::shared_ptr<const BraSpaceType> bra_space_ptr_;
      std::shared_ptr<const KetSpaceType> ket_space_ptr_;

      // sectors (accessible by index)
      std::vector<SectorType> sectors_;

      // sector index lookup by subspace indices
#ifdef BASIS_HASH
      std::unordered_map<
          typename SectorType::KeyType,std::size_t,
          boost::hash<typename SectorType::KeyType>
        > lookup_;
#else
      std::map<typename SectorType::KeyType,std::size_t> lookup_;
#endif

    };

  // Here we specialize to the case where the bra and ket spaces have the same
  // type. We inherit all the functionality from the case where they differ, and
  // add some additional type aliases and constructors.
  template<typename tSpaceType, typename tSectorType>
    class BaseSectors<tSpaceType, tSpaceType, tSectorType, true>
      : public BaseSectors<tSpaceType, tSpaceType, tSectorType, false>
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using SpaceType = tSpaceType;
      using SubspaceType = typename SpaceType::SubspaceType;

      protected:

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSectors() = default;

      template<
          typename T,
          typename std::enable_if_t<std::is_same_v<std::decay_t<T>, SpaceType>>* = nullptr
        >
      inline explicit BaseSectors(T&& space)
        : BaseSectors{std::forward<T>(space), std::forward<T>(space)}
      {}

      template<
          typename T, typename U,
          typename std::enable_if_t<std::is_same_v<std::decay_t<T>, SpaceType>>* = nullptr,
          typename std::enable_if_t<std::is_same_v<std::decay_t<U>, SpaceType>>* = nullptr
        >
      BaseSectors(T&& bra_space, U&& ket_space)
        : BaseSectors<tSpaceType, tSpaceType, tSectorType, false>{
              std::forward<T>(bra_space), std::forward<U>(ket_space)
            }
      {}

    };

  template <typename tBraSpaceType, typename tKetSpaceType, typename tSectorType>
    std::string BaseSectors<tBraSpaceType, tKetSpaceType, tSectorType, false>::DebugStr() const
    {
      std::ostringstream os;
      for (std::size_t sector_index=0; sector_index<size(); ++sector_index)
        {
          const SectorType& sector = GetSector(sector_index);

          os << "  sector " << sector_index
             << "  bra index " << sector.bra_subspace_index()
             << " labels " << sector.bra_subspace().LabelStr()
             << " dim " << sector.bra_subspace().size()
             << "  ket index " << sector.ket_subspace_index()
             << " labels " << sector.ket_subspace().LabelStr()
             << " dim " << sector.ket_subspace().size()
             << "  multiplicity index " << sector.multiplicity_index()
             << std::endl;
        }
      return os.str();
    }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
