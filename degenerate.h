/**************************************************************
  @file degenerate.h

  Defines basis containers in which each state has associated with it
  a substate degeneracy.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 06/06/17 (mac): Created as multibasis.h, abstracted from code in spncci
    branching_u3s.
  + 06/11/17 (mac): Rename TotalFullDimension to FullDimension.
  + 06/17/17 (mac): Rename to degenerate.h.  Rename multiplicity to degeneracy.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.
  + 06/24/21 (pjf):
    - Add dimension() accessor to BaseDegenerateSubspace.
    - Deprecate full_dimension() in BaseDegenerateSubspace and FullDimension()
      in BaseDegenerateSpace.
  + 06/25/21 (pjf):
    - Complete renaming of multiplicity to degeneracy.
    - Add GetStateOffset() and GetStateDegeneracy() to BaseDegenerateSubspace.
    - Repurpose BaseDegenerateSpace to fit its name:
      * Allow multiple degenerate subspaces within BaseDegenerateSpace.
      * Add subspace degeneracy.
      * Override PushSubspace and EmplaceSubspace to account for degeneracies.
  + 06/28/21 (pjf): Pass tSpaceLabelsType through BaseDegenerateSpace
    to BaseSpace.
  + 06/29/21 (pjf): Fix GetSubspaceOffset accessor.
  + 07/02/21 (pjf):
    - Add private convenience typedefs to access base class.
    - Replace `typedef`s with alias declarations (`using`).
    - Add constructors for BaseSubspace and BaseSpace which accept labels.
  + 07/03/21 (pjf): Update access specifiers.
  + 07/04/21 (pjf): Add additional template parameters to BaseDegenerateSubspace
    for pass-through to BaseSubspace.
  + 08/16/21 (pjf): Add additional template parameter to BaseDegenerateSpace for
    pass-through to BaseSpace.
  + 09/24/21 (pjf):
    - Create BaseDegenerateSector for sectors between subspaces in
      a BaseDegenerateSpace.
    - Make all degeneracy indices unsigned int.
****************************************************************/

#ifndef BASIS_DEGENERATE_H_
#define BASIS_DEGENERATE_H_

#include <cassert>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "basis.h"

namespace basis {


  ////////////////////////////////////////////////////////////////
  // generic subspace
  ////////////////////////////////////////////////////////////////

  // BaseDegenerateSubspace -- holds indexing of states within a symmetry
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

  template <
      typename tDerivedSubspaceType, typename tSubspaceLabelsType,
      typename tStateType, typename tStateLabelsType
    >
  class BaseDegenerateSubspace
    : public BaseSubspace<tDerivedSubspaceType,tSubspaceLabelsType,tStateType,tStateLabelsType>
  {

    private:

    ////////////////////////////////////////////////////////////////
    // private (convenience) typedefs
    ////////////////////////////////////////////////////////////////
    using BaseSubspaceType
      = BaseSubspace<tDerivedSubspaceType,tSubspaceLabelsType,tStateType,tStateLabelsType>;

    public:

    ////////////////////////////////////////////////////////////////
    //  common typedefs
    ////////////////////////////////////////////////////////////////

    // these typedefs apparently need to be repeated, not inherited
    // from BaseSubspace (they are not recognized as inherited types
    // below)

    using SubspaceLabelsType = typename BaseSubspaceType::SubspaceLabelsType;
    using StateLabelsType = typename BaseSubspaceType::StateLabelsType;

    protected:

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    // default constructor
    //   Implicitly invoked by derived class.
    BaseDegenerateSubspace() : full_dimension_(0) {}

    // pass-through constructor accepting labels
    explicit BaseDegenerateSubspace(const SubspaceLabelsType& labels)
      : full_dimension_{0}, BaseSubspaceType{labels}
    {}

    public:

    ////////////////////////////////////////////////////////////////
    // accessors for substate information
    ////////////////////////////////////////////////////////////////

    const std::vector<std::size_t>& state_offsets() const
    {
      return state_offsets_;
    }

    const std::vector<unsigned int>& state_degeneracies() const
    {
      return state_degeneracies_;
    }

    std::size_t GetStateOffset(std::size_t i, int degeneracy_index=1) const
    {
      assert(degeneracy_index <= state_degeneracies_.at(i));
      return state_offsets()[i]+(degeneracy_index-1);
    }

    unsigned int GetStateDegeneracy(std::size_t i) const
    {
      return state_degeneracies_.at(i);
    }

    std::size_t dimension() const
    {
      return full_dimension_;
    }

#ifdef BASIS_ALLOW_DEPRECATED
      DEPRECATED("use dimension() instead")
    std::size_t full_dimension() const
    {
      return full_dimension_;
    }
#endif

    protected:

    ////////////////////////////////////////////////////////////////
    // state label push (for initial construction)
    ////////////////////////////////////////////////////////////////

    void PushStateLabels(const StateLabelsType& state_labels, int degeneracy)
    // Create indexing information (in both directions, index <->
    // labels) for a state.
    {
      // push state itself
      BaseSubspaceType::PushStateLabels(state_labels);

      // push substate information
      state_offsets_.push_back(full_dimension_);
      state_degeneracies_.push_back(degeneracy);
      full_dimension_ += degeneracy;
    };

    private:

    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

    // degeneracy counting information
    std::vector<std::size_t> state_offsets_;  // offset to given state's starting substate
    std::vector<unsigned int> state_degeneracies_;  // given state's number of substates
    std::size_t full_dimension_;  // total number of substates

  };

  ////////////////////////////////////////////////////////////////
  // generic state realized within subspace
  ////////////////////////////////////////////////////////////////

  // BaseDegenerateState -- realization of a state withinin a given subspace
  //
  // The derived class is expected to set up a constructor and
  // friendlier accessors for the individual labels.
  //
  // The space (and the indexing it provides) is *not* copied into the
  // state but rather stored by pointer reference.  It should
  // therefore exist for the lifetime of the state object.
  //
  // Template arguments:
  //   tSubspaceType : subspace type in which this state lives

  template <typename tSubspaceType>
  class BaseDegenerateState
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

    protected:

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    // default constructor -- disabled
    BaseDegenerateState();

    // copy constructor -- synthesized

    // constructors

    BaseDegenerateState(const SubspaceType& subspace, std::size_t index)
      // Construct state, given index within subspace.
      : BaseStateType(subspace,index)
    {
    }

    BaseDegenerateState(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state, by reverse lookup on labels within subspace.
      : BaseStateType(subspace,state_labels)
      {
      }

    public:

    ////////////////////////////////////////////////////////////////
    // retrieval of substate information
    ////////////////////////////////////////////////////////////////

    std::size_t offset() const
    {
      return BaseStateType::subspace().state_offsets()[BaseStateType::index()];
    }

    unsigned int degeneracy() const
    {
      return BaseStateType::subspace().state_degeneracies()[BaseStateType::index()];
    }


  };



  ////////////////////////////////////////////////////////////////
  // generic space
  ////////////////////////////////////////////////////////////////

  // BaseDegenerateSpace -- container to hold subspaces, with reverse lookup by
  // subspace labels
  //
  // Template arguments:
  //   tSubspaceType (typename) : type for subspace

  template <typename tDerivedSpaceType, typename tSubspaceType, typename tSpaceLabelsType = void>
  class BaseDegenerateSpace
    : public BaseSpace<tDerivedSpaceType, tSubspaceType, tSpaceLabelsType>
  {

    private:

    ////////////////////////////////////////////////////////////////
    // private (convenience) typedefs
    ////////////////////////////////////////////////////////////////
    using BaseSpaceType = BaseSpace<tDerivedSpaceType, tSubspaceType,tSpaceLabelsType>;

    public:

    ////////////////////////////////////////////////////////////////
    // common typedefs
    ////////////////////////////////////////////////////////////////

    using SubspaceType = typename BaseSpaceType::SubspaceType;
    using SpaceLabelsType = typename BaseSpaceType::SpaceLabelsType;

    protected:

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    // default constructor
    //   Implicitly invoked by derived class.
    BaseDegenerateSpace() = default;

    // pass-through constructor accepting labels
    //   n.b. this template here is SFINAE black magic
    template<typename T = SpaceLabelsType>
    explicit BaseDegenerateSpace(const T& labels)
      : BaseSpaceType{labels}
    {}

    public:

    ////////////////////////////////////////////////////////////////
    // size retrieval
    ////////////////////////////////////////////////////////////////

#ifdef BASIS_ALLOW_DEPRECATED
    DEPRECATED("use dimension() instead")
    std::size_t FullDimension() const
    // Return the total dimension of all subspaces within the space,
    // taking into account substate degeneracies.
    {
      std::size_t full_dimension = 0;
      for (std::size_t subspace_index=0; subspace_index<BaseSpace<tSubspaceType>::size(); ++subspace_index)
        full_dimension += BaseSpace<tSubspaceType>::GetSubspace(subspace_index).full_dimension();
      return full_dimension;
    }
#endif

    ////////////////////////////////////////////////////////////////
    // subspace lookup and retrieval
    ////////////////////////////////////////////////////////////////

    unsigned int GetSubspaceDegeneracy(std::size_t i) const
    /// Given the index for a subspace, return the degeneracy of the
    /// subspace within the space.
    {
      return subspace_degeneracies_.at(i);
    };

    int GetSubspaceOffset(std::size_t i, int degeneracy_index) const
    /// Given the index for a subspace and the degeneracy index, return
    /// the offset of the subspace (with given degeneracy index) within
    /// the space.
    {
      assert(degeneracy_index <= subspace_degeneracies_.at(i));
      return (BaseSpaceType::GetSubspaceOffset(i)
        +(degeneracy_index-1)*this->GetSubspace(i).dimension());
    }

    protected:

    ////////////////////////////////////////////////////////////////
    // subspace push (for initial construction)
    ////////////////////////////////////////////////////////////////

    template<typename T, typename std::enable_if_t<std::is_same_v<std::decay_t<T>, SubspaceType>>* = nullptr>
    void PushSubspace(T&& subspace) = delete;
    // Prevent use of PushSubspace without degeneracy.

    template<typename T, typename std::enable_if_t<std::is_same_v<std::decay_t<T>, SubspaceType>>* = nullptr>
    void PushSubspace(T&& subspace, int degeneracy)
    /// Create indexing information (in both directions, index <->
    /// labels) for a subspace.
    {
      this->subspace_offsets_.push_back(this->dimension_);  // save offset
      (this->lookup_)[subspace.labels()] = BaseSpaceType::size();  // save index
      subspace_degeneracies_.push_back(degeneracy);  // save degeneracy
      this->dimension_ += subspace.dimension()*degeneracy;  // save dimension
      this->subspace_ptrs_.push_back(
          std::make_shared<const SubspaceType>(std::forward<T>(subspace))
        );  // save space
    };

    template <class... Args>
    void EmplaceSubspace(Args&&... args) = delete;
    // Prevent use of EmplaceSubspace without degeneracy.

    template <class... Args>
    void EmplaceSubspace(Args&&... args, int degeneracy)
    /// Create indexing information (in both directions, index <->
    /// labels) for a subspace.
    {
      const std::size_t index = BaseSpaceType::size();  // index for lookup
      this->subspace_offsets_.push_back(this->dimension_);  // save offset
      this->subspace_ptrs_.push_back(
          std::make_shared<SubspaceType>(std::forward<Args>(args)...)
        );  // construct/emplace space
      const SubspaceType& subspace = *(this->subspace_ptrs_.back());
      (this->lookup_)[subspace.labels()] = index;  // save index for lookup
      subspace_degeneracies_.push_back(degeneracy);  // save degeneracy
      this->dimension_ += subspace.dimension()*degeneracy;  // save dimension
    }

  private:

    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

    // degeneracy counting information
    std::vector<unsigned int> subspace_degeneracies_;  // given subspace's number of sub-subspaces

  };

  ////////////////////////////////////////////////////////////////
  // generic sector -- sector between subspaces which have both
  //   an index and degeneracy, i.e. live within a BaseDegenerateSpace
  ////////////////////////////////////////////////////////////////

  // BaseDegenerateSector

  template<typename tBraSubspaceType, typename tKetSubspaceType = tBraSubspaceType>
  class BaseDegenerateSector
    : public BaseSector<tBraSubspaceType, tKetSubspaceType>
  {
    private:

    ////////////////////////////////////////////////////////////////
    // private (convenience) typedefs
    ////////////////////////////////////////////////////////////////
    using BaseSectorType = BaseSector<tBraSubspaceType, tKetSubspaceType>;

    public:

    ////////////////////////////////////////////////////////////////
    // common typedefs
    ////////////////////////////////////////////////////////////////

    using BraSubspaceType = typename BaseSectorType::BraSubspaceType;
    using KetSubspaceType = typename BaseSectorType::KetSubspaceType;

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    BaseDegenerateSector(
        std::size_t bra_subspace_index, std::size_t ket_subspace_index,
        unsigned int bra_subspace_degeneracy, unsigned int ket_subspace_degeneracy,
        const BraSubspaceType& bra_subspace, const KetSubspaceType& ket_subspace,
        std::size_t multiplicity_index=1
      )
      : BaseSectorType{
            bra_subspace_index, ket_subspace_index,
            bra_subspace, ket_subspace,
            multiplicity_index
          },
        bra_subspace_degeneracy_{bra_subspace_degeneracy},
        ket_subspace_degeneracy_{ket_subspace_degeneracy}
    {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline unsigned int bra_subspace_degeneracy() const
    {
      return bra_subspace_degeneracy_;
    }

    inline unsigned int ket_subspace_degeneracy() const
    {
      return ket_subspace_degeneracy_;
    }

    inline std::size_t num_elements() const
    {
      return (bra_subspace_degeneracy_ * BaseSectorType::bra_subspace().dimension())
        * (ket_subspace_degeneracy_ * BaseSectorType::ket_subspace().dimension());
    }

    private:

    unsigned int bra_subspace_degeneracy_, ket_subspace_degeneracy_;

  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
}  // namespace basis

#endif  // BASIS_DEGENERATE_H_
