/****************************************************************
  jjjttz_test.cpp

  Zhou Zhou
  University of Notre Dame

****************************************************************/

#include <iomanip>
#include <iostream>


#include "basis/jjjttz_scheme.h"
#include "basis/basis.h"
#include "basis/many_body.h"


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
  basis::TwoBodySubspaceJJJTTz subspace(0,0,0,0,basis::Rank::kTwoBody,4);  // JTTzg Nmax
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();


  // spaces and sectors

  // first set up space

  std::cout << "Two-body space" << std::endl;
  int Nmax = 4;
  basis::TwoBodySpaceJJJTTz space(basis::Rank::kTwoBody,Nmax);
  std::cout << space.DebugStr();

  for (const auto& subspace : space)
  {
    std::cout << "subspace: " << subspace.LabelStr() << std::endl;
    std::cout << subspace.DebugStr();
  }

  // then set up allowed sectors
  std::cout << "Two-body operator sectors" << std::endl;
  int J0 = 0;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  int g0 = 0;
  int Tz0 = 0;
  basis::TwoBodySectorsJJJTTz sectors(space,J0,g0,Tz0);

  std::cout << " J0 " << J0 << " g0 " << g0 << " Tz0 " << Tz0 << std::endl;
  for (std::size_t sector_index=0; sector_index < sectors.size(); ++sector_index)
    {
      std::size_t bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
      const basis::TwoBodySubspaceJJJTTz& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
      std::size_t ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
      const basis::TwoBodySubspaceJJJTTz& ket_subspace = sectors.GetSector(sector_index).ket_subspace();

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
        << std::setw(3) << bra_subspace.Tz()
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
        << std::setw(3) << ket_subspace.Tz()
        // << std::setw(3) << ket_subspace.Nmax()
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

  // termination
  return 0;
}
