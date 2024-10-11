/************************************************************//**
  @file degenerate_test.cpp

  Language: C++11

  Mark A. Caprio, Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include "basis/basis.h"
#include "basis/product.h"
#include "basis/nlj_orbital.h"

#include "am/am.h"


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

class TestProductState;

class TestProductSubspace
  : public basis::BaseProductSubspace<
      TestProductSubspace, TestProductState, basis::OrbitalSubspacePN, basis::OrbitalSubspacePN
    >
{
 public:
  using BaseProductSubspaceType = basis::BaseProductSubspace<
      TestProductSubspace, TestProductState, basis::OrbitalSubspacePN, basis::OrbitalSubspacePN
    >;

  TestProductSubspace() = default;

  explicit TestProductSubspace(
      std::shared_ptr<const basis::OrbitalSubspacePN> subspace1_ptr,
      std::shared_ptr<const basis::OrbitalSubspacePN> subspace2_ptr

    )
    : BaseProductSubspace{std::move(subspace1_ptr), std::move(subspace2_ptr)}
  {}
};

class TestProductSpace
  : public basis::BaseProductSpace<TestProductSpace, TestProductSubspace, basis::OrbitalSpacePN, basis::OrbitalSpacePN>
{
 public:
  using BaseProductSpaceType = basis::BaseProductSpace<
      TestProductSpace, TestProductSubspace, basis::OrbitalSpacePN, basis::OrbitalSpacePN
    >;

  TestProductSpace() = default;

  explicit TestProductSpace(
      std::shared_ptr<const basis::OrbitalSpacePN> space1_ptr,
      std::shared_ptr<const basis::OrbitalSpacePN> space2_ptr

    )
    : BaseProductSpaceType{std::move(space1_ptr), std::move(space2_ptr)}
  {}
};

void Test()
{
  int N_max = 2;
  int J_max = 3;

  auto space1 = std::make_shared<basis::OrbitalSpacePN>(N_max);
  auto space2 = std::make_shared<basis::OrbitalSpacePN>(N_max);

  auto subspace = TestProductSubspace(space1->GetSubspacePtr(0),space2->GetSubspacePtr(0));
  fmt::print("{}\n", subspace.dimension());

  auto space = TestProductSpace(space1,space2);
  fmt::print("{}\n", space.dimension());
}

int main(int argc, char **argv)
{
  Test();

  // termination
  return EXIT_SUCCESS;
}
