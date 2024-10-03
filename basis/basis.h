/**************************************************************/
/**
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

    BASIS_HASH: In lookup tables, replace std::map with std::unordered_map.
    BASIS_BOOST_HASH: Use Boost.ContainerHash extensions for lookup tables.
    BASIS_STD_HASH: Use STL hashing (std::hash) for lookup tables.

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
  + 03/05/24 (mac): Add alias LabelsType to BaseState.
  + 09/25/24 (pjf):
    - Split basis.h into separate headers.
    - Add BASIS_BOOST_HASH and BASIS_STD_HASH directives.
****************************************************************/

#ifndef BASIS_BASIS_H_
#define BASIS_BASIS_H_

#include "basis/map.h"
#include "basis/subspace.h"
#include "basis/state.h"
#include "basis/space.h"
#include "basis/sector.h"

#endif  // BASIS_BASIS_H_
