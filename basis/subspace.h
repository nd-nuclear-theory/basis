/****************************************************************

  subspace.h

  Language: C++17

  Mark A. Caprio
  University of Notre Dame
  Patrick J. Fasano
  Argonne National Laboratory

  + 09/25/24: Split basis.h into separate headers.
  + 09/27/24: Move implementations into impl namespace to hide template
    parameters which should be invisible to library users.

****************************************************************/

#ifndef BASIS_SUBSPACE_H_
#define BASIS_SUBSPACE_H_

#include <cassert>
#include <cstddef>
#include <iterator>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef BASIS_ALLOW_DEPRECATED
// emit warnings on deprecated
#include "mcutils/deprecated.h"
#endif

#include "basis/map.h"
#include "basis/type_traits.h"

namespace basis {

namespace impl {

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

    using LabelsType = void;
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
      using iterator_category = std::input_iterator_tag;  // we don't satisfy the requirements of LegacyForwardIterator because we return temporaries
      using iterator_concept = std::random_access_iterator_tag;  // C++20 iterators have weaker requirements
      using difference_type = std::make_signed_t<std::size_t>;
      using value_type = tStateType;
      using pointer = void*;
      // using pointer = value_type*;  // n.b. we can't return a pointer to a temporary
      using reference = value_type;  // n.b. this can't be a reference since value_type is always a temporary

      iterator() = default;

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
      friend iterator operator+(difference_type n, iterator it) { it += n; return it; }
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
#ifdef BASIS_HASH
      lookup_.reserve(new_cap);
#endif
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
    basis::map<StateLabelsType,std::size_t> lookup_;

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

}  // namespace impl

template <
    typename tDerivedSubspaceType, typename tSubspaceLabelsType,
    typename tStateType, typename tStateLabelsType
  >
using BaseSubspace
  = impl::BaseSubspace<tDerivedSubspaceType, tSubspaceLabelsType, tStateType, tStateLabelsType>;

}  // namespace basis

#endif  // BASIS_SUBSPACE_H_
