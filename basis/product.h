/****************************************************************

  Define template base classes for indexing quantum mechanical product
  states arranged into subspaces.

  Language: C++20

  Patrick J. Fasano
  Argonne National Laboratory

  + 09/20/24: Created.

****************************************************************/

#ifndef BASIS_PRODUCT_H_
#define BASIS_PRODUCT_H_

#include <memory>
#include <tuple>
#include <utility>

#include "basis.h"

namespace basis {

namespace impl {

////////////////////////////////////////////////////////////////
// indexing helper functions
////////////////////////////////////////////////////////////////
template<typename... Ts>
constexpr auto ind2sub(std::size_t index, Ts... s)
{
  std::array<std::size_t, sizeof...(Ts)> out;
  auto shape = std::array<std::size_t, sizeof...(Ts)>{std::size_t(s)...};
  for (int i=sizeof...(Ts)-1; i>=0; --i)
  {
    out[i] = index % shape[i];
    index -= out[i];
    index /= shape[i];
  }
  return out;
}

template<typename... Ts>
constexpr auto sub2ind(std::array<std::size_t, sizeof...(Ts)> index, Ts... s)
{
  std::size_t out = 0;
  auto shape = std::array<std::size_t, sizeof...(Ts)>{std::size_t(s)...};
  for (int i=0; i<sizeof...(Ts); ++i)
  {
    out = out*shape[i] + index[i];
  }
  return out;
}

template<
    typename tDerivedSubspaceType, typename tStateType, bool bWithLabels,
    typename tSequence, typename... tFactorSubspaceTypes
  >
class BaseProductSubspace;

template<
    typename tDerivedSubspaceType, typename tStateType,
    std::size_t... Is, typename... tFactorSubspaceTypes
  >
class BaseProductSubspace<tDerivedSubspaceType, tStateType, false, std::index_sequence<Is...>, tFactorSubspaceTypes...>
  // : public BaseSubspace<
  //     tDerivedSubspaceType,
  //     std::tuple<typename tFactorSubspaceTypes::SubspaceLabelsType...>,
  //     tStateType,
  //     std::tuple<typename tFactorSubspaceTypes::StateLabelsType...>
  //   >
{
 private:
  ////////////////////////////////////////////////////////////////
  // private (convenience) typedefs
  ////////////////////////////////////////////////////////////////
  // using BaseSubspaceType
  //   = BaseSubspace<
  //       tDerivedSubspaceType, std::tuple<typename tFactorSubspaceTypes::SubspaceLabelsType...>,
  //       tStateType, std::tuple<typename tFactorSubspaceTypes::StateLabelsType...>
  //     >;

 public:
  ////////////////////////////////////////////////////////////////
  //  common typedefs
  ////////////////////////////////////////////////////////////////

  using StateLabelsType = std::tuple<typename tFactorSubspaceTypes::StateLabelsType...>;
  using StateLabelsRef = std::tuple<const typename tFactorSubspaceTypes::StateLabelsType&...>;

  using StateMultiIndexType = std::array<std::size_t, sizeof...(tFactorSubspaceTypes)>;
  using StateType = tStateType;


 protected:
  // default constructor
  //   Implicitly invoked by derived class.
  BaseProductSubspace() = default;

  // pass-through constructor accepting labels
  //   n.b. the dummy argument "BST" is SFINAE black magic
  //   see: https://stackoverflow.com/a/13401982
  explicit BaseProductSubspace(
      std::shared_ptr<const tFactorSubspaceTypes>... factor_subspace_ptrs
    )
    : //BaseSubspaceType{std::make_tuple(factor_subspace_ptrs->labels()...)},
      factor_subspace_ptrs_{std::move(factor_subspace_ptrs)...}
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

  template<std::size_t I>
  auto GetFactorSubspacePtr() const
  {
    return std::get<I>(factor_subspace_ptrs_);
  }

  template<std::size_t I>
  const auto& GetFactorSubspace() const
  {
    return *std::get<I>(factor_subspace_ptrs_);
  }

  StateMultiIndexType IndexToMultiIndex(std::size_t index) const
  {
    return (
        (index != kNone)
        ? ind2sub(index, GetFactorSubspace<Is>().size()...)
        : StateMultiIndexType{((void)Is, kNone)...}  // (void) to suppress unused warning
      );
  }

  std::size_t MultiIndexToIndex(const StateMultiIndexType& indices) const
  {
    return (
        ((std::get<Is>(indices) != kNone) && ...)
        ? sub2ind(indices, GetFactorSubspace<Is>().size()...)
        : kNone
      );
  }

  StateLabelsRef GetStateLabels(const StateMultiIndexType& indices) const
  {
    // construct tuple from state labels
    return std::forward_as_tuple(GetFactorSubspace<Is>().GetStateLabels(indices[Is])...);
  }

  StateLabelsRef GetStateLabels(std::size_t index) const
  {
    return GetStateLabels(IndexToMultiIndex(index));
  }

  bool ContainsState(const StateLabelsType& state_labels) const
  /// Given the labels for a state, returns whether or not the state
  /// is found within the subspace.
  {
    return (GetFactorSubspace<Is>().ContainsState(std::get<Is>(state_labels)) && ...);
  }

  StateMultiIndexType LookUpStateMultiIndex(const StateLabelsType& state_labels) const
  /// Given the labels for a state, look up its index within the
  /// subspace.
  ///
  /// If no such labels are found, basis::kNone is returned.
  {
    return StateMultiIndexType{
        GetFactorSubspace<Is>().LookUpStateIndex(std::get<Is>(state_labels))...
      };
  }

  std::size_t LookUpStateIndex(const StateLabelsType& state_labels) const
  /// Given the labels for a state, look up its index within the
  /// subspace.
  ///
  /// If no such labels are found, basis::kNone is returned.
  {
    return MultiIndexToIndex(LookUpStateMultiIndex(state_labels));
  }

  StateType GetState(std::size_t index) const;
  /// Given the index for a state, construct the corresponding state object.
  ///
  /// Defined below, outside the definition of the class, so that
  /// instantiation is deferred until tStateType is a complete type.
  StateType GetState(const StateMultiIndexType& indices) const;
  /// Given the multi-index for a state, construct the corresponding state object.
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
  {
    // fold size over multiplication
    return (GetFactorSubspace<Is>().size() * ...);
  }

  std::size_t dimension() const
  {
    // fold dimension over multiplication
    return (GetFactorSubspace<Is>().dimension() * ...);
  }

 private:
  std::tuple<std::shared_ptr<const tFactorSubspaceTypes>...> factor_subspace_ptrs_;
};

// implementations of BaseSubspace::GetState; we do this outside the class
// declaration to delay instantiation until tDerivedSubspaceType and
// tStateType are complete types

template<
    typename tDerivedSubspaceType, typename tStateType,
    std::size_t... Is, typename... tFactorSubspaceTypes
  >
tStateType
BaseProductSubspace<tDerivedSubspaceType,tStateType,false,std::index_sequence<Is...>,tFactorSubspaceTypes...>::GetState(
    std::size_t index
  ) const
{
  return tStateType{*static_cast<const tDerivedSubspaceType*>(this), index};
}

template<
    typename tDerivedSubspaceType, typename tStateType,
    std::size_t... Is, typename... tFactorSubspaceTypes
  >
tStateType
BaseProductSubspace<tDerivedSubspaceType,tStateType,false,std::index_sequence<Is...>,tFactorSubspaceTypes...>::GetState(
    const StateMultiIndexType& indices
  ) const
{
  return tStateType{*static_cast<const tDerivedSubspaceType*>(this), indices};
}

template<
    typename tDerivedSubspaceType, typename tStateType,
    std::size_t... Is, typename... tFactorSubspaceTypes
  >
tStateType
BaseProductSubspace<tDerivedSubspaceType,tStateType,false,std::index_sequence<Is...>,tFactorSubspaceTypes...>::GetState(
    const StateLabelsType& state_labels
  ) const
{
  return tStateType{*static_cast<const tDerivedSubspaceType*>(this), state_labels};
}


// inherit from BaseSubspace and add labels
template<
    typename tDerivedSubspaceType, typename tStateType,
    std::size_t... Is, typename... tFactorSubspaceTypes
  >
class BaseProductSubspace<tDerivedSubspaceType, tStateType, true, std::index_sequence<Is...>, tFactorSubspaceTypes...>
  : public BaseProductSubspace<tDerivedSubspaceType, tStateType, false, std::index_sequence<Is...>, tFactorSubspaceTypes...>
{
 public:
  ////////////////////////////////////////////////////////////////
  //  common typedefs
  ////////////////////////////////////////////////////////////////

  using LabelsType =
      std::tuple<typename tFactorSubspaceTypes::SubspaceLabelsType...>;
  using SubspaceLabelsType = LabelsType;
  using LabelsRef =
      std::tuple<const typename tFactorSubspaceTypes::SubspaceLabelsType&...>;

 protected:
  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////

  /// default constructor
  //   Implicitly invoked by derived class.
  BaseProductSubspace() = default;

  // inherit constructors
  using BaseProductSubspace<
      tDerivedSubspaceType, tStateType, false, std::index_sequence<Is...>, tFactorSubspaceTypes...
    >::BaseProductSubspace;

  ////////////////////////////////////////////////////////////////
  // retrieval
  ////////////////////////////////////////////////////////////////

  LabelsRef labels() const
  /// Return the labels of the subspace itself.
  {
    // fold size over multiplication
    return std::forward_as_tuple(BaseProductSubspace::template GetFactorSubspace<Is>().labels()...);
  }
};

////////////////////////////////////////////////////////////////
// generic product state realized within subspace
////////////////////////////////////////////////////////////////

/// BaseProductState -- realization of a state withinin a given subspace
template<typename tSubspaceType>
class BaseProductState
  : public BaseState<tSubspaceType>
{
 private:
  ////////////////////////////////////////////////////////////////
  // private (convenience) typedefs
  ////////////////////////////////////////////////////////////////
  using BaseStateType = BaseState<tSubspaceType>;

 public:
  ////////////////////////////////////////////////////////////////
  // common typedefs
  ////////////////////////////////////////////////////////////////

  using SubspaceType = typename BaseStateType::SubspaceType;
  using StateLabelsType = typename BaseStateType::StateLabelsType;
  using StateMultiIndexType = typename SubspaceType::StateMultiIndexType;

 protected:
  ////////////////////////////////////////////////////////////////
  // general constructors
  ////////////////////////////////////////////////////////////////

  // default constructor -- disabled
  BaseProductState() = delete;

  // copy constructor -- synthesized

  // constructors

  BaseProductState(const SubspaceType& subspace, const StateMultiIndexType& indices)
    // Construct state, given index within subspace.
    : BaseStateType(subspace,subspace.MultiIndexToIndex(indices)), multiindex_{indices}
  {}

  BaseProductState(const SubspaceType& subspace, std::size_t index)
    // Construct state, given index within subspace.
    : BaseStateType(subspace,index), multiindex_{subspace.IndexToMultiIndex(index)}
  {}

  BaseProductState(const SubspaceType& subspace, const StateLabelsType& state_labels)
    // Construct state, by reverse lookup on labels within subspace.
    : BaseStateType(subspace,state_labels)
  {
    multiindex_ = BaseStateType::subspace().IndexToMultiIndex(BaseStateType::index());
  }

 public:
    ////////////////////////////////////////////////////////////////
    // retrieval
    ////////////////////////////////////////////////////////////////

  const StateMultiIndexType& indices() const { return multiindex_; }

  using BaseStateType::index;
  template<std::size_t I>
  decltype(auto) index() const
  {
    return std::get<I>(indices());
  }

  using BaseStateType::subspace;
  template<std::size_t I>
  decltype(auto) subspace() const
  {
    return BaseStateType::subspace().template GetFactorSubspace<I>();
  }

  using BaseStateType::labels;
  template<std::size_t I>
  decltype(auto) labels() const
  {
    return subspace<I>().GetStateLabels(index<I>());
  }

 private:
  StateMultiIndexType multiindex_;
};

////////////////////////////////////////////////////////////////
// generic space
////////////////////////////////////////////////////////////////

/// BaseProductSpace -- container to hold product subspaces, with reverse lookup
/// by subspace labels
///
/// Note: This class behaves more like BaseSubspace than BaseSpace: product
/// subspaces are not stored, but rather constructed "on the fly" by
/// GetSubspace(). Of course, the BaseProductSubspace type is very lightweight,
/// since it only contains a tuple of factor subspaces.
///
/// Template arguments: tSubspaceType (typename) : type for subspace
/// tSpaceLabelsType (typename, optional) : type for space labels

template<std::size_t Depth, std::size_t I=0>
struct NestedLoop
{
  template<typename F, typename... Xs>
  inline constexpr void operator()(F&& f, const std::array<std::size_t,Depth>& bounds, Xs... xs)
  {
    for (std::size_t i = 0; i < bounds[I]; ++i)
    {
      NestedLoop<Depth,I+1>{}(std::forward<F>(f), bounds, xs..., i);
    }
  }
};

template<std::size_t Depth>
struct NestedLoop<Depth,Depth>
{
  template<typename F, typename... Xs>
  inline constexpr void operator()(F&& f, const std::array<std::size_t,Depth>& bounds, Xs... xs)
  {
    f(xs...);
  }
};

// declare BaseSpace template
template <
  typename tDerivedSpaceType, typename tSubspaceType, bool bWithLabels,
  typename tSequence, typename... tFactorSpaceTypes
  >
class BaseProductSpace;

// specialize class for case with no labels
template <typename tDerivedSpaceType, typename tSubspaceType, std::size_t... Is, typename... tFactorSpaceTypes>
class BaseProductSpace<tDerivedSpaceType, tSubspaceType, false, std::index_sequence<Is...>, tFactorSpaceTypes...>
  : public BaseSpace<tDerivedSpaceType, tSubspaceType, void>
{
 public:
  ////////////////////////////////////////////////////////////////
  // common typedefs
  ////////////////////////////////////////////////////////////////

  using SubspaceType = tSubspaceType;
  using SubspaceMultiIndexType = std::array<std::size_t, sizeof...(tFactorSpaceTypes)>;
  using SubspaceLabelsType = std::tuple<typename tFactorSpaceTypes::SubspaceLabelsType...>;
  using SubspaceLabelsRef
    = std::tuple<const typename tFactorSpaceTypes::SubspaceLabelsType&...>;

 protected:
  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////

  BaseProductSpace() = default;

  explicit BaseProductSpace(std::shared_ptr<const tFactorSpaceTypes>... factor_spaces)
    : factor_space_ptrs_{std::move(factor_spaces)...}
  {
    auto bounds = std::array<std::size_t,sizeof...(Is)>{this->GetFactorSpace<Is>().size()...};
    NestedLoop<sizeof...(tFactorSpaceTypes)>{}(
        [&,this](auto... indices) {
            this->subspace_offsets_.push_back(this->dimension_);
            this->dimension_ += (this->GetFactorSpace<Is>().GetSubspace(indices).dimension() * ...);
          },
          bounds
      );
    if ((GetFactorSpace<Is>().dimension() * ...) != this->dimension_)
      throw std::logic_error("dimension mismatch");
  }

 public:
  ////////////////////////////////////////////////////////////////
  // retrieval
  ////////////////////////////////////////////////////////////////

  template<std::size_t I>
  auto GetFactorSpacePtr() const
  {
    return std::get<I>(factor_space_ptrs_);
  }

  template<std::size_t I>
  const auto& GetFactorSpace() const
  {
    return *std::get<I>(factor_space_ptrs_);
  }

  SubspaceMultiIndexType IndexToMultiIndex(std::size_t index) const
  {
    return (
        (index != kNone)
        ? ind2sub(index, GetFactorSpace<Is>().size()...)
        : SubspaceMultiIndexType{((void)Is, kNone)...}  // (void) to suppress unused warning
      );
  }

  std::size_t MultiIndexToIndex(const SubspaceMultiIndexType& indices) const
  {
    return (
        ((std::get<Is>(indices) != kNone) && ...)
        ? sub2ind(indices, GetFactorSpace<Is>().size()...)
        : kNone
      );
  }

  bool ContainsSubspace(const SubspaceLabelsType& subspace_labels) const
  /// Given the labels for a state, returns whether or not the state
  /// is found within the subspace.
  {
    return (GetFactorSpace<Is>().ContainsSubspace(std::get<Is>(subspace_labels)) && ...);
  }

  SubspaceMultiIndexType LookUpSubspaceMultiIndex(const SubspaceLabelsType& subspace_labels) const
  /// Given the labels for a subspace, look up its multi-index within the
  /// space.
  ///
  /// If no such labels are found, basis::kNone is returned.
  {
    return SubspaceMultiIndexType{
        GetFactorSpace<Is>().LookUpSubspaceIndex(std::get<Is>(subspace_labels))...
      };
  }

  std::size_t LookUpSubspaceIndex(const SubspaceLabelsType& subspace_labels) const
  /// Given the labels for a state, look up its index within the
  /// subspace.
  ///
  /// If no such labels are found, basis::kNone is returned.
  {
    return MultiIndexToIndex(LookUpSubspaceMultiIndex(subspace_labels));
  }

  SubspaceType LookUpSubspace(const SubspaceLabelsType& subspace_labels) const;
  /// Given the labels for a subspace, construct the corresponding subspace object.
  ///
  /// Defined below, outside the definition of the class, so that
  /// instantiation is deferred until SubspaceType is a complete type.

  SubspaceType GetSubspace(const SubspaceMultiIndexType& indices) const;
  /// Given the multi-index for a subspace, construct the corresponding subspace object.
  ///
  /// Defined below, outside the definition of the class, so that
  /// instantiation is deferred until SubspaceType is a complete type.
  SubspaceType GetSubspace(std::size_t index) const;
  /// Given the index for a subspace, construct the corresponding subspace object.
  ///
  /// Defined below, outside the definition of the class, so that
  /// instantiation is deferred until SubspaceType is a complete type.

  ////////////////////////////////////////////////////////////////
  // size retrieval
  ////////////////////////////////////////////////////////////////

  std::size_t size() const
  {
    return this->subspace_offsets_.size();
  }

 private:
  std::tuple<std::shared_ptr<const tFactorSpaceTypes>...> factor_space_ptrs_;
};

// implementations of BaseSpace::GetSubspace; we do this outside the class
// declaration to delay instantiation until tDerivedSubspaceType and
// tStateType are complete types

template <typename tDerivedSpaceType, typename tSubspaceType, std::size_t... Is, typename... tFactorSpaceTypes>
tSubspaceType BaseProductSpace<tDerivedSpaceType,tSubspaceType,false,std::index_sequence<Is...>, tFactorSpaceTypes...>::GetSubspace(
    const SubspaceMultiIndexType& indices
  ) const
{
  return tSubspaceType{
        GetFactorSpace<Is>().GetSubspacePtr(std::get<Is>(indices))...
      };
}

template <typename tDerivedSpaceType, typename tSubspaceType, std::size_t... Is, typename... tFactorSpaceTypes>
tSubspaceType BaseProductSpace<tDerivedSpaceType,tSubspaceType,false,std::index_sequence<Is...>, tFactorSpaceTypes...>::GetSubspace(
    std::size_t index
  ) const
{
  auto indices = IndexToMultiIndex(index);
  return GetSubspace(indices);
}

template <typename tDerivedSpaceType, typename tSubspaceType, std::size_t... Is, typename... tFactorSpaceTypes>
tSubspaceType BaseProductSpace<tDerivedSpaceType,tSubspaceType,false,std::index_sequence<Is...>, tFactorSpaceTypes...>::LookUpSubspace(
    const SubspaceLabelsType& subspace_labels
  ) const
{
  auto indices = LookUpSubspaceMultiIndex(subspace_labels);
  return GetSubspace(indices);
}


// inherit from BaseSpace and add labels
template <typename tDerivedSpaceType, typename tSubspaceType, std::size_t... Is, typename... tFactorSpaceTypes>
class BaseProductSpace<tDerivedSpaceType, tSubspaceType, true, std::index_sequence<Is...>, tFactorSpaceTypes...>
  : public BaseProductSpace<tDerivedSpaceType, tSubspaceType, false, std::index_sequence<Is...>, tFactorSpaceTypes...>
{
 private:

  ////////////////////////////////////////////////////////////////
  // private (convenience) typedefs
  ////////////////////////////////////////////////////////////////
  using BaseProductSpaceType = BaseProductSpace<tDerivedSpaceType, tSubspaceType, false, std::index_sequence<Is...>, tFactorSpaceTypes...>;

 public:

  ////////////////////////////////////////////////////////////////
  // common typedefs
  ////////////////////////////////////////////////////////////////

  using LabelsType = std::tuple<typename tFactorSpaceTypes::LabelsType...>;
  using SpaceLabelsType = LabelsType;
  using LabelsRef = std::tuple<const typename tFactorSpaceTypes::LabelsType&...>;

 protected:

  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////

  BaseProductSpace() = default;
  using BaseProductSpaceType::BaseProductSpaceType;

 public:

  ////////////////////////////////////////////////////////////////
  // retrieval
  ////////////////////////////////////////////////////////////////

  LabelsRef labels() const
  /// Return the labels of the subspace itself.
  {
    return std::forward_as_tuple(BaseProductSpace::template GetFactorSpace<Is>().labels()...);
  }
};

}  // namespace impl

template<
    typename tDerivedSubspaceType, typename tStateType,
    typename... tFactorSubspaceTypes
  >
using BaseProductSubspace
  = impl::BaseProductSubspace<
      tDerivedSubspaceType, tStateType,
      !(std::is_void_v<typename tFactorSubspaceTypes::LabelsType> || ...),
      std::index_sequence_for<tFactorSubspaceTypes...>,
      tFactorSubspaceTypes...
    >;

template<typename tSubspaceType>
using BaseProductState = impl::BaseProductState<tSubspaceType>;

template<
    typename tDerivedSpaceType, typename tSubspaceType,
    typename... tFactorSpaceTypes
  >
using BaseProductSpace
  = impl::BaseProductSpace<
      tDerivedSpaceType, tSubspaceType,
      !(std::is_void_v<typename tFactorSpaceTypes::LabelsType> || ...),
      std::index_sequence_for<tFactorSpaceTypes...>,
      tFactorSpaceTypes...
    >;
}  // namespace basis

#endif  // BASIS_PRODUCT_H_
