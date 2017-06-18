/**************************************************************
  @file degenerate.h

  Defines basis containers in which each state has associated with it
  a substate degeneracy.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 6/6/17 (mac): Created as multibasis.h, abstracted from code in spncci
    branching_u3s.
  + 6/11/17 (mac): Rename TotalFullDimension to FullDimension.
  + 6/17/17 (mac): Rename to degenerate.h.  Rename multiplicity to degeneracy.

****************************************************************/

#ifndef BASIS_DEGENERATE_H_
#define BASIS_DEGENERATE_H_

#include "basis/basis.h"

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

  template <typename tSubspaceLabelsType, typename tStateLabelsType>
    class BaseDegenerateSubspace
    : public BaseSubspace<tSubspaceLabelsType,tStateLabelsType>
  {

    public:

    ////////////////////////////////////////////////////////////////
    //  common typedefs
    ////////////////////////////////////////////////////////////////

    // these typedefs apparently need to be repeated, not inherited
    // from BaseSubspace (they are not recognized as inherited types
    // below)

    typedef tSubspaceLabelsType SubspaceLabelsType;
    typedef tStateLabelsType StateLabelsType;

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    // default constructor
    //   Implicitly invoked by derived class.
    BaseDegenerateSubspace() : full_dimension_(0) {}

    ////////////////////////////////////////////////////////////////
    // accessors for substate information
    ////////////////////////////////////////////////////////////////

    const std::vector<int>& state_offsets() const
    {
      return state_offsets_;
    }

    const std::vector<int>& state_multiplicities() const
    {
      return state_multiplicities_;
    }

    int full_dimension() const
    {
      return full_dimension_;
    }

    protected:

    ////////////////////////////////////////////////////////////////
    // state label push (for initial construction)
    ////////////////////////////////////////////////////////////////

    void PushStateLabels(const StateLabelsType& state_labels, int degeneracy)
    // Create indexing information (in both directions, index <->
    // labels) for a state.
    {
      // push state itself
      BaseSubspace<tSubspaceLabelsType,tStateLabelsType>::PushStateLabels(state_labels);

      // push substate information
      state_offsets_.push_back(full_dimension_);
      state_multiplicities_.push_back(degeneracy);
      full_dimension_ += degeneracy;
    };

    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

    // degeneracy counting information
    std::vector<int> state_offsets_;  // offset to given state's starting substate
    std::vector<int> state_multiplicities_;  // given state's number of substates
    int full_dimension_;  // total number of substates

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
      BaseDegenerateState();

      // copy constructor -- synthesized

      // constructors

      BaseDegenerateState(const SubspaceType& subspace, int index)
        // Construct state, given index within subspace.
        : BaseState<tSubspaceType>(subspace,index)
      {
      }

      BaseDegenerateState(const SubspaceType& subspace, const StateLabelsType& state_labels)
        // Construct state, by reverse lookup on labels within subspace.
        : BaseState<tSubspaceType>(subspace,state_labels)
        {
        }

      ////////////////////////////////////////////////////////////////
      // retrieval of substate information
      ////////////////////////////////////////////////////////////////

      int offset() const
      {
        return BaseState<tSubspaceType>::subspace().state_offsets()[BaseState<tSubspaceType>::index()];
      }

      int degeneracy() const
      {
        return BaseState<tSubspaceType>::subspace().state_multiplicities()[BaseState<tSubspaceType>::index()];
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

  template <typename tSubspaceType>
    class BaseDegenerateSpace
    : public BaseSpace<tSubspaceType>
    {

      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      typedef tSubspaceType SubspaceType;

      ////////////////////////////////////////////////////////////////
      // size retrieval
      ////////////////////////////////////////////////////////////////

      int FullDimension() const
      // Return the total dimension of all subspaces within the space,
      // taking into account substate multiplicities.
      {
        int full_dimension = 0;
        for (int subspace_index=0; subspace_index<BaseSpace<tSubspaceType>::size(); ++subspace_index)
          full_dimension += BaseSpace<tSubspaceType>::GetSubspace(subspace_index).full_dimension();
        return full_dimension;
      }

    };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
