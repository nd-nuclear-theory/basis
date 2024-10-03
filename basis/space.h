/****************************************************************

  space

  Language: C++17

  Mark A. Caprio
  University of Notre Dame
  Patrick J. Fasano
  Argonne National Laboratory

  + 09/25/24: Split basis.h into separate headers.
  + 09/27/24: Move implementations into impl namespace to hide template
    parameters which should be invisible to library users.

****************************************************************/

#ifndef BASIS_SPACE_H_
#define BASIS_SPACE_H_

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

      using LabelsType = void;
      using SubspaceType = tSubspaceType;
      using SubspaceLabelsType = typename SubspaceType::SubspaceLabelsType;

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
#ifdef BASIS_HASH
      lookup_.reserve(new_cap);
#endif
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
      template<typename, typename, bool, typename...> friend class BaseProductSpace;

      ////////////////////////////////////////////////////////////////
      // internal storage
      ////////////////////////////////////////////////////////////////

      /// subspaces (accessible by index)
      std::shared_ptr<std::vector<SubspaceType>> subspaces_ptr_;
      std::vector<std::size_t> subspace_offsets_;
      std::size_t dimension_;

      // subspace index lookup by labels
      basis::map<typename SubspaceType::LabelsType,std::size_t> lookup_;
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

}  // namespace impl

template <typename tDerivedSpaceType, typename tSubspaceType, typename tSpaceLabelsType = void>
using BaseSpace = impl::BaseSpace<tDerivedSpaceType, tSubspaceType, tSpaceLabelsType>;


}  // namespace basis

#endif  // BASIS_SPACE_H_
