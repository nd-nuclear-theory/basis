/****************************************************************
  jjjt_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <iomanip>
#include <iostream>


#include "jjjt_scheme.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void TestTwoBody()
{

  ////////////////////////////////////////////////////////////////
  // two-body basis tests
  ////////////////////////////////////////////////////////////////

  std::cout << "Two-body basis" << std::endl;

  // example subspace
  std::cout << "Example subspace" << std::endl;
  basis::TwoBodySubspaceJJJT subspace(0,1,0,basis::Rank::kTwoBody,6);  // JTg Nmax
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();


  // spaces and sectors

  // first set up space

  std::cout << "Two-body space" << std::endl;
  int Nmax = 2;
  basis::TwoBodySpaceJJJT space(basis::Rank::kTwoBody,Nmax);
  std::cout << space.DebugStr();

  // then set up allowed sectors
  std::cout << "Two-body operator sectors" << std::endl;
  int J0 = 2;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  int T0 = 0;
  int g0 = 0;
  basis::TwoBodySectorsJJJT sectors(space,J0,T0,g0);

  std::cout << " J0 " << J0 << " T0 " << T0 << " g0 " << g0 << std::endl;
  for (std::size_t sector_index=0; sector_index < sectors.size(); ++sector_index)
    {
      std::size_t bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
      const basis::TwoBodySubspaceJJJT& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
      std::size_t ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
      const basis::TwoBodySubspaceJJJT& ket_subspace = sectors.GetSector(sector_index).ket_subspace();

      std::cout
        << " sector "
        << std::setw(3) << sector_index
        << "     "
        << " index "
        << std::setw(3) << bra_subspace_index
        << " JTg "
        << std::setw(3) << bra_subspace.J()
        << std::setw(3) << bra_subspace.T()
        << std::setw(3) << bra_subspace.g()
        // << std::setw(3) << bra_subspace.Nmax()
        << " dim "
        << std::setw(3) << bra_subspace.size()
        << "     "
        << " index "
        << std::setw(3) << ket_subspace_index
        << " JTg "
        << std::setw(3) << ket_subspace.J()
        << std::setw(3) << ket_subspace.T()
        << std::setw(3) << ket_subspace.g()
        // << std::setw(3) << ket_subspace.Nmax()
        << " dim "
        << std::setw(3) << ket_subspace.size()
        << std::endl;
    }
}

void TestTwoBodyN()
{

  ////////////////////////////////////////////////////////////////
  // two-body basis tests -- subspaced by N
  ////////////////////////////////////////////////////////////////

  std::cout << "Two-body basis -- subspaced by N" << std::endl;

  // example subspace
  std::cout << "Example subspace" << std::endl;
  basis::TwoBodySubspaceJJJTN subspace(0,1,0,6,basis::Rank::kTwoBody,6);  // JTgN Nmax
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();


  // spaces and sectors

  // first set up space

  std::cout << "Two-body space" << std::endl;
  int Nmax = 2;
  basis::TwoBodySpaceJJJTN space(basis::Rank::kTwoBody,Nmax);
  std::cout << space.DebugStr();

  // then set up allowed sectors
  std::cout << "Two-body operator sectors" << std::endl;
  int J0 = 2;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  int T0 = 0;
  int g0 = 0;
  basis::TwoBodySectorsJJJTN sectors(space,J0,T0,g0);

  std::cout << " J0 " << J0 << " T0 " << T0 << " g0 " << g0 << std::endl;
  for (std::size_t sector_index=0; sector_index < sectors.size(); ++sector_index)
    {
      std::size_t bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
      const basis::TwoBodySubspaceJJJTN& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
      std::size_t ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
      const basis::TwoBodySubspaceJJJTN& ket_subspace = sectors.GetSector(sector_index).ket_subspace();

      std::cout
        << " sector "
        << std::setw(3) << sector_index
        << "     "
        << " index "
        << std::setw(3) << bra_subspace_index
        << " JTg-N "
        << std::setw(3) << bra_subspace.J()
        << std::setw(3) << bra_subspace.T()
        << std::setw(3) << bra_subspace.g()
        << std::setw(3) << bra_subspace.N()
        << " dim "
        << std::setw(3) << bra_subspace.size()
        << "     "
        << " index "
        << std::setw(3) << ket_subspace_index
        << " JTg-N "
        << std::setw(3) << ket_subspace.J()
        << std::setw(3) << ket_subspace.T()
        << std::setw(3) << ket_subspace.g()
        << std::setw(3) << ket_subspace.N()
        << " dim "
        << std::setw(3) << ket_subspace.size()
        << std::endl;
    }
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  TestTwoBody();
  TestTwoBodyN();

  // termination
  return 0;
}
