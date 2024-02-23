/************************************************************//**
  @file hypersector_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <cstdlib>
#include <algorithm>
#include <iostream>

#include "hypersector.h"

#include "am/am.h"
#include "basis.h"
#include "lsjt_scheme.h"


////////////////////////////////////////////////////////////////
// example hypersectors implementation
////////////////////////////////////////////////////////////////

// acts on RelativeSpaceLSJT, with operators also labeled by
// RelativeSpaceLSJT

namespace basis {

  // declaration

  class RelativeHypersectorsLSJT
    : public BaseHypersectors<RelativeSpaceLSJT,RelativeSpaceLSJT>
  {

    public:

    // constructor

    RelativeHypersectorsLSJT() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    RelativeHypersectorsLSJT(
        const RelativeSpaceLSJT& space,
        const RelativeSpaceLSJT& operator_space,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

  };

  // definition

  RelativeHypersectorsLSJT::RelativeHypersectorsLSJT(
      const RelativeSpaceLSJT& space,
      const RelativeSpaceLSJT& operator_space,
      basis::SectorDirection sector_direction
    )
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        for (int operator_subspace_index=0; operator_subspace_index<operator_space.size(); ++operator_subspace_index)
          {
            // enforce canonical ordering
            if (
                (sector_direction == basis::SectorDirection::kCanonical)
                && !(bra_subspace_index<=ket_subspace_index)
              )
              continue;

            // retrieve subspaces
            const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
            const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);
            const OperatorSubspaceType& operator_subspace = operator_space.GetSubspace(operator_subspace_index);

            // verify angular momentum, isosopin, and parity selection rules
            bool allowed = true;
            allowed &= am::AllowedTriangle(ket_subspace.L(),operator_subspace.L(),bra_subspace.L());
            allowed &= am::AllowedTriangle(ket_subspace.S(),operator_subspace.S(),bra_subspace.S());
            allowed &= am::AllowedTriangle(ket_subspace.J(),operator_subspace.J(),bra_subspace.J());
            allowed &= am::AllowedTriangle(ket_subspace.T(),operator_subspace.T(),bra_subspace.T());
            allowed &= ((ket_subspace.g()+operator_subspace.g()+bra_subspace.g())%2==0);

            // push hypersector
            if (allowed)
              PushHypersector(
                  bra_subspace_index,ket_subspace_index,operator_subspace_index
                );
          }
  }

} // namespace


void TestHyprersectors()
{

  std::cout << "Relative space" << std::endl;
  int N_max = 2;
  int J_max = 2;
  basis::RelativeSpaceLSJT space(N_max,J_max);
  std::cout << space.DebugStr();

  std::cout << "Relative operators" << std::endl;
  int N0_max = 0;
  int J0_max = 0;
  basis::RelativeSpaceLSJT operator_space(N0_max,J0_max);
  std::cout << operator_space.DebugStr();

  // then set up allowed hypersectors
  std::cout << "Relative operator hypersectors" << std::endl;
  basis::RelativeHypersectorsLSJT hypersectors(space,operator_space);
  std::cout << hypersectors.DebugStr();

  // allocate zero operator
  basis::OperatorHyperblocks<float> matrices;
  basis::SetHyperoperatorToZero(hypersectors,matrices);

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  TestHyprersectors();

  // termination
  return EXIT_SUCCESS;
}
