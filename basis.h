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

  Library dependencies:

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
  + 08/18/21 (pjf): Add std::vector-like size management to BaseSubspace and
    BaseSpace.
  + 09/22/21 (pjf):
    - Improve reserve() in BaseSubspace and BaseSpace.
    - Add iterator support to BaseSubspace, BaseSpace, and BaseSectors.
  + 09/24/21 (pjf):
    - Fix use-after-move in BaseSectors::PushBack.
    - Add indexing offsets for BaseSectors, and num_elements accessor to
      both BaseSectors and BaseSector.
    - Fix BaseSectors::DebugStr for size/dimension distinction.
  + 09/30/21 (pjf): Add BaseSpace::full_size() accessor.
  + 11/04/21 (pjf):
    - Revert to shared pointer to vector in BaseSpace.
    - Remove std::enable_shared_from_this.
    - Add GetSubspacePtr() accessor to BaseSpace using shared_ptr aliasing.
    - Make BaseSector take shared pointers instead of const references.
  + 03/25/22 (pjf):
    - Modify BaseSubspace to allow definition with `void` labels.
  + 04/02/22 (pjf):
    - Defer initialization of shared_ptr in BaseSpace until
      PushSubspace/EmplaceSubspace/reserve.
    - Make labels in BaseSubspace and BaseSpace non-const, to avoid deleting
      copy and move constructors.
  + 04/03/22 (pjf): Normalize all constructors so that fields get initialized
    correctly.
  + 04/12/22 (pjf): Allow BaseSectors between non-direct subspace of space.
  + 02/03/23 (mac): Add alias StateType to BaseSubspace.
  + 05/09/23 (pjf): Fix comparison of pointers in BaseSectors.
  + 05/12/23 (pjf):
    - Fix iterator types.
    - Fix constness of member shared_ptr in BaseSector.
    - Use std::addressof for comparison of objects in BaseSectors.
****************************************************************/

#ifndef BASIS_BASIS_H_
#define BASIS_BASIS_H_

#include <cassert>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// for DebugStr
#include <iosfwd>
#include <sstream>
#include <string>

#ifdef BASIS_ALLOW_DEPRECATED
// emit warnings on deprecated
#include "mcutils/deprecated.h"
#endif

#ifdef BASIS_HASH
#include <unordered_map>
#include <boost/container_hash/hash.hpp>
#else
#include <map>
#endif

namespace basis {

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

  // declare BaseSubspace template
  template <
      typename tDerivedSubspaceType, typename tSubspaceLabelsType,
      typename tStateType, typename tStateLabelsType
    >
  class BaseSubspace;

  // specialize class for case with no labels
  template <
      typename tDerivedSubspaceType,
      typename tStateType, typename tStateLabelsType
    >
  class BaseSubspace<tDerivedSubspaceType, void, tStateType, tStateLabelsType>
  {

    public:

    ////////////////////////////////////////////////////////////////
    //  common typedefs
    ////////////////////////////////////////////////////////////////

    using StateLabelsType = tStateLabelsType;
    using StateType = tStateType;

    protected:

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    /// default constructor
    //   Implicitly invoked by derived class.
    BaseSubspace()
      : dimension_{}
    {}

    public:

    ////////////////////////////////////////////////////////////////
    // iterators
    ////////////////////////////////////////////////////////////////

    struct iterator
    {
      using iterator_category = std::input_iterator_tag;  // we don't satisfy the requirements of forward_iterator_tag because we return temporaries
      using difference_type = std::make_signed_t<std::size_t>;
      using value_type = tStateType;
      using pointer = void*;
      // using pointer = value_type*;  // n.b. we can't return a pointer to a temporary
      using reference = value_type;  // n.b. this can't be a reference since value_type is always a temporary

      iterator(const tDerivedSubspaceType* subspace_ptr, std::size_t index)
        : subspace_ptr_{subspace_ptr}, index_{index}
      {}

      reference operator*() const { return subspace_ptr_->GetState(index_); }
      reference operator[](std::size_t index) const { return subspace_ptr_->GetState(index_+index); }
      iterator& operator++() { index_++; return *this; }
      iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }
      iterator& operator--() { index_--; return *this; }
      iterator operator--(int) { iterator tmp = *this; --(*this); return tmp; }
      iterator& operator+=(difference_type n) { index_ += n; return *this; }
      iterator& operator-=(difference_type n) { index_ -= n; return *this; }
      friend iterator operator+(iterator it, difference_type n) { it += n; return it; }
      friend iterator operator-(iterator it, difference_type n) { it -= n; return it; }
      friend difference_type operator-(const iterator& lhs, const iterator& rhs)
      {
        return (lhs.index_ - rhs.index_);
      }
      friend bool operator==(const iterator& lhs, const iterator& rhs)
      {
        return (lhs.subspace_ptr_==rhs.subspace_ptr_) && (lhs.index_==rhs.index_);
      }
      friend bool operator!=(const iterator& lhs, const iterator& rhs)
      {
        return (lhs.subspace_ptr_!=rhs.subspace_ptr_) || (lhs.index_!=rhs.index_);
      }
      friend bool operator<(const iterator& lhs, const iterator& rhs)
      {
        return (lhs.index_<rhs.index_);
      }
      friend bool operator>(const iterator& lhs, const iterator& rhs)
      {
        return (lhs.index_>rhs.index_);
      }
      friend bool operator<=(const iterator& lhs, const iterator& rhs)
      {
        return (lhs.index_<=rhs.index_);
      }
      friend bool operator>=(const iterator& lhs, const iterator& rhs)
      {
        return (lhs.index_>=rhs.index_);
      }

      private:
      const tDerivedSubspaceType* subspace_ptr_;
      std::size_t index_;
    };

    using value_type = tStateType;
    using const_reference = const tStateType&;
    using const_iterator = iterator;
    using difference_type = typename iterator::difference_type;
    using size_type = std::make_unsigned_t<difference_type>;

    iterator begin() const
    {
      return iterator{static_cast<const tDerivedSubspaceType*>(this), 0};
    }
    iterator end() const
    {
      return iterator{static_cast<const tDerivedSubspaceType*>(this), size()};
    }
    iterator cbegin() const
    {
      return iterator{static_cast<const tDerivedSubspaceType*>(this), 0};
    }
    iterator cend() const
    {
      return iterator{static_cast<const tDerivedSubspaceType*>(this), size()};
    }

    ////////////////////////////////////////////////////////////////
    // retrieval
    ////////////////////////////////////////////////////////////////

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

    StateType GetState(std::size_t index) const;
    /// Given the index for a state, construct the corresponding state object.
    ///
    /// Defined below, outside the definition of the class, so that
    /// instantiation is deferred until tStateType is a complete type.

    StateType GetState(const StateLabelsType& state_labels) const;
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
    // state label storage management (for initial construction)
    ////////////////////////////////////////////////////////////////

    inline void reserve(std::size_t new_cap)
    /// Reserve storage for labels.
    {
      state_table_.reserve(new_cap);
      lookup_.reserve(new_cap);
    }

    inline std::size_t capacity() const noexcept
    /// Get size of reserved storage for labels.
    {
      return state_table_.capacity();
    }

    inline void shrink_to_fit()
    /// Shrink reserved storage for labels to fit contents.
    {
      state_table_.shrink_to_fit();
    }

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
      typename tDerivedSubspaceType,
      typename tStateType, typename tStateLabelsType
    >
  tStateType BaseSubspace<tDerivedSubspaceType,void,tStateType,tStateLabelsType>::GetState(
      std::size_t index
    ) const
  {
    return tStateType{*static_cast<const tDerivedSubspaceType*>(this), index};
  }

  template <
      typename tDerivedSubspaceType,
      typename tStateType, typename tStateLabelsType
    >
  tStateType BaseSubspace<tDerivedSubspaceType,void,tStateType,tStateLabelsType>::GetState(
      const tStateLabelsType& state_labels
    ) const
  {
    return tStateType{*static_cast<const tDerivedSubspaceType*>(this), state_labels};
  }

  // inherit from BaseSubspace and add labels
  template <
      typename tDerivedSubspaceType, typename tSubspaceLabelsType,
      typename tStateType, typename tStateLabelsType
    >
  class BaseSubspace : public BaseSubspace<tDerivedSubspaceType, void, tStateType, tStateLabelsType>
  {

    public:

    ////////////////////////////////////////////////////////////////
    //  common typedefs
    ////////////////////////////////////////////////////////////////

    using LabelsType = tSubspaceLabelsType;
    using SubspaceLabelsType = tSubspaceLabelsType;

    protected:

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    /// default constructor
    //   Implicitly invoked by derived class.
    BaseSubspace() = default;

    // pass-through constructor accepting labels
    template<
        typename T = SubspaceLabelsType,
        std::enable_if_t<std::is_constructible_v<SubspaceLabelsType, T>>* = nullptr
      >
    explicit BaseSubspace(T&& labels)
      : BaseSubspace<tDerivedSubspaceType, void, tStateType, tStateLabelsType>{},
        labels_{std::forward<T>(labels)}
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

    private:

    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

    // subspace properties
    SubspaceLabelsType labels_;
  };

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
        : dimension_{0}
      {}

      public:

      ////////////////////////////////////////////////////////////////
      // iterators
      ////////////////////////////////////////////////////////////////

      struct iterator
      {
        using iterator_category = std::random_access_iterator_tag;
        using difference_type = std::make_signed_t<std::size_t>;
        using value_type = const SubspaceType;
        using pointer = const SubspaceType*;
        using reference = const SubspaceType&;

        iterator(const tDerivedSpaceType* space_ptr, std::size_t index)
          : space_ptr_{space_ptr}, index_{index}
        {}

        reference operator*() const { return space_ptr_->GetSubspace(index_); }
        pointer operator->() const { return &(*this); }
        reference operator[](std::size_t index) const { return space_ptr_->GetSubspace(index_+index); }
        iterator& operator++() { index_++; return *this; }
        iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }
        iterator& operator--() { index_--; return *this; }
        iterator operator--(int) { iterator tmp = *this; --(*this); return tmp; }
        iterator& operator+=(difference_type n) { index_ += n; return *this; }
        iterator& operator-=(difference_type n) { index_ -= n; return *this; }
        friend iterator operator+(iterator it, difference_type n) { it += n; return it; }
        friend iterator operator-(iterator it, difference_type n) { it -= n; return it; }
        friend difference_type operator-(const iterator& lhs, const iterator& rhs)
        {
          return (lhs.index_ - rhs.index_);
        }
        friend bool operator==(const iterator& lhs, const iterator& rhs)
        {
          return (lhs.space_ptr_==rhs.space_ptr_) && (lhs.index_==rhs.index_);
        }
        friend bool operator!=(const iterator& lhs, const iterator& rhs)
        {
          return (lhs.space_ptr_!=rhs.space_ptr_) || (lhs.index_!=rhs.index_);
        }
        friend bool operator<(const iterator& lhs, const iterator& rhs)
        {
          return (lhs.index_<rhs.index_);
        }
        friend bool operator>(const iterator& lhs, const iterator& rhs)
        {
          return (lhs.index_>rhs.index_);
        }
        friend bool operator<=(const iterator& lhs, const iterator& rhs)
        {
          return (lhs.index_<=rhs.index_);
        }
        friend bool operator>=(const iterator& lhs, const iterator& rhs)
        {
          return (lhs.index_>=rhs.index_);
        }

        private:
        const tDerivedSpaceType* space_ptr_;
        std::size_t index_;
      };

      using value_type = tSubspaceType;
      using const_reference = const tSubspaceType&;
      using const_iterator = iterator;
      using difference_type = typename iterator::difference_type;
      using size_type = std::make_unsigned_t<difference_type>;

      iterator begin()  const { return iterator(static_cast<const tDerivedSpaceType*>(this), 0); }
      iterator end()    const { return iterator(static_cast<const tDerivedSpaceType*>(this), size()); }
      iterator cbegin() const { return iterator(static_cast<const tDerivedSpaceType*>(this), 0); }
      iterator cend()   const { return iterator(static_cast<const tDerivedSpaceType*>(this), size()); }

      ////////////////////////////////////////////////////////////////
      // subspace lookup and retrieval
      ////////////////////////////////////////////////////////////////

      bool ContainsSubspace(const typename SubspaceType::LabelsType& subspace_labels) const
      /// Given the labels for a subspace, returns whether or not the
      /// subspace is found within the space.
      {
        return (subspaces_ptr_) && lookup_.count(subspace_labels);
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

        // short-circuit if subspaces_ptr_ is nullptr
        if (!subspaces_ptr_)
          return kNone;

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
      /// If no such labels are found, an exception will result.
      {
        if (!subspaces_ptr_)
          throw std::out_of_range("empty space");
        std::size_t subspace_index = LookUpSubspaceIndex(subspace_labels);
        if (subspace_index==kNone)
          throw std::out_of_range("key not found");
        return subspaces_ptr_->at(subspace_index);
      };

      const SubspaceType& GetSubspace(std::size_t i) const
      /// Given the index for a subspace, return a reference to the
      /// subspace.
      {
        if (!subspaces_ptr_)
          throw std::out_of_range("empty space");
        return subspaces_ptr_->at(i);
      };

      std::shared_ptr<const SubspaceType> GetSubspacePtr(std::size_t i) const
      /// Given the index for a subspace, return a shared pointer to the
      /// subspace.
      {
        if (!subspaces_ptr_)
          throw std::out_of_range("empty space");
        return std::shared_ptr<const SubspaceType>(
            subspaces_ptr_, &(subspaces_ptr_->at(i))
          );
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
        if (!subspaces_ptr_)
          return 0;
        return subspaces_ptr_->size();
      };

      inline std::size_t full_size() const { return size(); }
      /// Get the total number of subspaces, considering degeneracies.

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
    // subspace storage management (for initial construction)
    ////////////////////////////////////////////////////////////////

    inline void reserve(std::size_t new_cap)
    /// Reserve storage for subspaces and labels.
    {
      if (!subspaces_ptr_)
        subspaces_ptr_ = std::make_shared<std::vector<SubspaceType>>();
      subspaces_ptr_->reserve(new_cap);
      subspace_offsets_.reserve(new_cap);
      lookup_.reserve(new_cap);
    }

    inline std::size_t capacity() const noexcept
    /// Get current storage reserved for subspaces.
    {
      if (!subspaces_ptr_)
        return 0;
      return subspaces_ptr_->capacity();
    }

    inline void shrink_to_fit()
    /// Shrink reserved storage for labels and offsets to fit contents.
    {
      if (subspaces_ptr_)
        subspaces_ptr_->shrink_to_fit();
      subspace_offsets_.shrink_to_fit();
    }

      ////////////////////////////////////////////////////////////////
      // subspace push (for initial construction)
      ////////////////////////////////////////////////////////////////

      template<typename T, typename std::enable_if_t<std::is_same_v<std::decay_t<T>, SubspaceType>>* = nullptr>
      void PushSubspace(T&& subspace)
      /// Create indexing information (in both directions, index <->
      /// labels) for a subspace.
      {
        if (!subspaces_ptr_)
          subspaces_ptr_ = std::make_shared<std::vector<SubspaceType>>();
        subspace_offsets_.push_back(dimension_);
        lookup_[subspace.labels()] = size();  // index for lookup
        dimension_ += subspace.dimension();
        subspaces_ptr_->push_back(std::forward<T>(subspace));  // save subspace
      };

      template <class... Args>
      void EmplaceSubspace(Args&&... args)
      /// Create indexing information (in both directions, index <->
      /// labels) for a subspace.
      {
        if (!subspaces_ptr_)
          subspaces_ptr_ = std::make_shared<std::vector<SubspaceType>>();
        const std::size_t index = size();  // index for lookup
        subspace_offsets_.push_back(dimension_);
        subspaces_ptr_->emplace_back(std::forward<Args>(args)...);
        const SubspaceType& subspace = subspaces_ptr_->back();
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
      std::shared_ptr<std::vector<SubspaceType>> subspaces_ptr_;
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
      template<
          typename T = SpaceLabelsType,
          std::enable_if_t<std::is_constructible_v<SpaceLabelsType, T>>* = nullptr
        >
      explicit BaseSpace(T&& labels)
        : BaseSpace<tDerivedSpaceType, tSubspaceType,void>{},
          labels_{std::forward<T>(labels)}
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
      SpaceLabelsType labels_;
    };

  ////////////////////////////////////////////////////////////////
  // type trait helpers
  //
  // for use with SFINAE
  ////////////////////////////////////////////////////////////////

  template<typename, typename, typename = void>
  struct is_subspace
      : std::false_type
  {};

  template<typename Subspace, typename Space>
  struct is_subspace<
      Subspace, Space,
      std::enable_if_t<std::is_same_v<Subspace, typename Space::SubspaceType>>
    >
      : std::true_type
  {};

  template<typename Subspace, typename Space>
  struct is_subspace<
      Subspace, Space,
      std::enable_if_t<
          is_subspace<Subspace, typename Space::SubspaceType>::value
        >
    >
      : std::true_type
  {};

  template<typename T, typename U>
  constexpr bool is_subspace_v = is_subspace<T, U>::value;


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

      BaseSector()
        : bra_subspace_index_{}, ket_subspace_index_{}, multiplicity_index_{}
      {}
      // default constructor

      BaseSector(
          std::size_t bra_subspace_index, std::size_t ket_subspace_index,
          std::shared_ptr<const BraSubspaceType> bra_subspace_ptr,
          std::shared_ptr<const KetSubspaceType> ket_subspace_ptr,
          std::size_t multiplicity_index=1
        )
        : bra_subspace_index_{bra_subspace_index},
          ket_subspace_index_{ket_subspace_index},
          multiplicity_index_{multiplicity_index},
          bra_subspace_ptr_{std::move(bra_subspace_ptr)},
          ket_subspace_ptr_{std::move(ket_subspace_ptr)}
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

      std::shared_ptr<const BraSubspaceType> bra_subspace_ptr() const
      {
        return bra_subspace_ptr_;
      }
      std::shared_ptr<const KetSubspaceType> ket_subspace_ptr() const
      {
        return ket_subspace_ptr_;
      }
      // Return shared pointer to bra/ket subspace.

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

      inline std::size_t num_elements() const
      {
        return (bra_subspace().dimension() * ket_subspace().dimension());
      }

      private:
      std::size_t bra_subspace_index_, ket_subspace_index_;
      std::size_t multiplicity_index_;
      std::shared_ptr<const BraSubspaceType> bra_subspace_ptr_;
      std::shared_ptr<const KetSubspaceType> ket_subspace_ptr_;
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

      BaseSector() = default;

      inline BaseSector(
          std::size_t bra_subspace_index, std::size_t ket_subspace_index,
          std::shared_ptr<const SubspaceType> bra_subspace_ptr,
          std::shared_ptr<const SubspaceType> ket_subspace_ptr,
          std::size_t multiplicity_index=1
        )
        : BaseSector<tSubspaceType, tSubspaceType, false>{
            bra_subspace_index, ket_subspace_index,
            std::move(bra_subspace_ptr), std::move(ket_subspace_ptr),
            multiplicity_index
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

      using SectorType = tSectorType;
      using BraSubspaceType = typename SectorType::BraSubspaceType;
      static_assert(
          is_subspace_v<BraSubspaceType,BraSpaceType>,
          "BraSubspaceType must be subspace of BraSpaceType"
        );
      using KetSubspaceType = typename SectorType::KetSubspaceType;
      static_assert(
          is_subspace_v<KetSubspaceType,KetSpaceType>,
          "KetSubspaceType must be subspace of KetSpaceType"
        );


      protected:

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSectors()
        : num_elements_{}
      {}
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      BaseSectors(
            const std::shared_ptr<const BraSpaceType>& bra_space_ptr,
            const std::shared_ptr<const KetSpaceType>& ket_space_ptr
          )
          : bra_space_ptr_{bra_space_ptr},
            ket_space_ptr_{ket_space_ptr},
            num_elements_{}
      {}

      BaseSectors(
            std::shared_ptr<const BraSpaceType>&& bra_space_ptr,
            std::shared_ptr<const KetSpaceType>&& ket_space_ptr
          )
          : bra_space_ptr_{std::move(bra_space_ptr)},
            ket_space_ptr_{std::move(ket_space_ptr)},
            num_elements_{}
      {}

      template<
          typename T, typename U,
          typename std::enable_if_t<std::is_same_v<std::decay_t<T>, BraSpaceType>>* = nullptr,
          typename std::enable_if_t<std::is_same_v<std::decay_t<U>, KetSpaceType>>* = nullptr
        >
      BaseSectors(T&& bra_space, U&& ket_space)
        // : bra_space_ptr_(bra_space.weak_from_this().lock()),
        //   ket_space_ptr_(ket_space.weak_from_this().lock())
        : num_elements_{}
      {
        // if bra_space or ket_space is not managed by a shared_ptr, make copy
        // or move it to a locally-created shared_ptr.
        //
        // Note: this is probably safe, since the spaces themselves store
        // vectors of shared_ptrs to subspaces, so a copy only involves copying
        // those vectors and maps
        // if (!bra_space_ptr_)
          bra_space_ptr_
            = std::make_shared<const BraSpaceType>(std::forward<T>(bra_space));
        // if (!ket_space_ptr_)
        {
          // if bra and ket spaces are identical, we can share the new copy we
          // just made for the bra_space_ptr_
          //
          // note: this also ensures that if bra_space and ket_space are the
          // same object, then if we just moved from bra_space we don't try to
          // move from it again (which would be UB)
          if (std::addressof(ket_space) == std::addressof(bra_space))
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

      std::shared_ptr<const BraSpaceType> bra_space_ptr() const
      {
        return bra_space_ptr_;
      }

      std::shared_ptr<const KetSpaceType> ket_space_ptr() const
      {
        return ket_space_ptr_;
      }

      ////////////////////////////////////////////////////////////////
      // iterators
      ////////////////////////////////////////////////////////////////

      using value_type = SectorType;
      using const_reference = const SectorType&;
      using iterator = decltype(std::declval<const std::vector<SectorType>>().begin());
      using const_iterator = decltype(std::declval<const std::vector<SectorType>>().cbegin());
      using difference_type = typename iterator::difference_type;
      using size_type = std::make_unsigned_t<difference_type>;

      iterator begin()  const { return sectors_.begin(); }
      iterator end()    const { return sectors_.end(); }
      const_iterator cbegin() const { return sectors_.cbegin(); }
      const_iterator cend()   const { return sectors_.cend(); }

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

      std::size_t GetSectorOffset(std::size_t i) const
      /// Given the index for a sector, return the offset of the sector in
      /// the full matrix indexing.
      {
        return sector_offsets_.at(i);
      }

      ////////////////////////////////////////////////////////////////
      // size retrieval
      ////////////////////////////////////////////////////////////////

      std::size_t size() const
      // Return number of sectors within sector set.
      {
        return sectors_.size();
      };

      std::size_t num_elements() const
      // Return number of matrix elements within sector set
      {
        return num_elements_;
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
        sector_offsets_.push_back(num_elements());
        lookup_[sector.Key()] = size(); // index for lookup
        num_elements_ += sector.num_elements();
        sectors_.push_back(std::forward<T>(sector));  // save sector
      };

      template<typename... Args>
      void EmplaceSector(Args&&... args)
      // Create indexing information (in both directions, index <->
      // labels) for a sector.
      {
        const std::size_t index = sectors_.size();
        sector_offsets_.push_back(num_elements());
        sectors_.emplace_back(std::forward<Args>(args)...);  // save sector
        const SectorType& sector = sectors_.back();
        lookup_[sector.Key()] = index; // index for lookup
        num_elements_ += sector.num_elements();
      };

      template<
          typename T = SectorType,
          std::enable_if_t<std::is_same_v<typename T::BraSubspaceType,typename BraSpaceType::SubspaceType>>* = nullptr,
          std::enable_if_t<std::is_same_v<typename T::KetSubspaceType,typename KetSpaceType::SubspaceType>>* = nullptr
        >
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
            bra_space_ptr_->GetSubspacePtr(bra_subspace_index),
            ket_space_ptr_->GetSubspacePtr(ket_subspace_index),
            multiplicity_index
          );  // save sector
      }

      template<
          typename T = SectorType,
          std::enable_if_t<std::is_same_v<typename T::BraSubspaceType,typename BraSpaceType::SubspaceType>>* = nullptr,
          std::enable_if_t<std::is_same_v<typename T::KetSubspaceType,typename KetSpaceType::SubspaceType>>* = nullptr
        >
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

      // indexing
      std::vector<std::size_t> sector_offsets_;
      std::size_t num_elements_;

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
      // TODO (pjf): Fix up template type here so that SubspaceType only
      // gets defined if unambiguous
      // template<
      //     typename T = tSectorType,
      //     std::enable_if_t<std::is_same_v<typename T::BraSubspaceType, typename T::KetSubspaceType>>* = nullptr
      //   >
      static_assert(
          std::is_same_v<typename tSectorType::BraSubspaceType, typename tSectorType::KetSubspaceType>,
          "same-Space BaseSectors currently allowed with same-Subspace BaseSector"
        );
      using SubspaceType = typename tSectorType::BraSubspaceType;

      protected:

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSectors() = default;

      inline explicit BaseSectors(std::shared_ptr<SpaceType> space_ptr)
        : BaseSectors{space_ptr, std::move(space_ptr)}
      {}
      // construct from shared_ptr to space

      inline explicit BaseSectors(
          std::shared_ptr<SpaceType> bra_space_ptr,
          std::shared_ptr<SpaceType> ket_space_ptr
        )
        : BaseSectors<tSpaceType, tSpaceType, tSectorType, false>{
              std::move(bra_space_ptr), std::move(ket_space_ptr)
            }
      {}

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
             << " size " << sector.bra_subspace().size()
             << " dim " << sector.bra_subspace().dimension()
             << "  ket index " << sector.ket_subspace_index()
             << " labels " << sector.ket_subspace().LabelStr()
             << " size " << sector.ket_subspace().size()
             << " dim " << sector.ket_subspace().dimension()
             << "  multiplicity index " << sector.multiplicity_index()
             << "  elements " << sector.num_elements()
             << std::endl;
        }
      return os.str();
    }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
