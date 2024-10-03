/****************************************************************

  state.h

  Language: C++17

  Mark A. Caprio
  University of Notre Dame
  Patrick J. Fasano
  Argonne National Laboratory

  + 09/25/24: Split basis.h into separate headers.
  + 09/27/24: Move implementations into impl namespace to hide template
    parameters which should be invisible to library users.

****************************************************************/

#ifndef BASIS_STATE_H_
#define BASIS_STATE_H_

#include <cassert>
#include <cstddef>
#include <type_traits>
#include <utility>

#ifdef BASIS_ALLOW_DEPRECATED
// emit warnings on deprecated
#include "mcutils/deprecated.h"
#endif

#include "basis/type_traits.h"

namespace basis {

namespace impl {

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
      using LabelsType = typename SubspaceType::StateLabelsType;
      using StateLabelsType = typename SubspaceType::StateLabelsType;

      protected:

      ////////////////////////////////////////////////////////////////
      // general constructors
      ////////////////////////////////////////////////////////////////

      /// default constructor -- disabled
      BaseState() = delete;

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

}  // namespace impl

template <typename tSubspaceType>
using BaseState = impl::BaseState<tSubspaceType>;

}  // namespace basis

#endif  // BASIS_STATE_H_
