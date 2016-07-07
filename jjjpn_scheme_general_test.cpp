/****************************************************************
  jjjpn_scheme_general_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <iomanip>
#include <iostream>

#include "jjjpn_scheme_general.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void test_orbitals_Nmax(const std::string& filename)
{

  ////////////////////////////////////////////////////////////////
  // orbital tests (should come after suborbital tests are passed)
  ////////////////////////////////////////////////////////////////

  std::cout << "Orbitals -- Nmax scheme" << std::endl;

  // set up space
  std::cout << "Space" << std::endl;
  int Nmax = 4;
  basis::OrbitalSpacePN space(Nmax);
  std::cout << space.DebugStr();

  // check subspaces
  std::cout << "Subspaces" << std::endl;
  for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
    {
      const basis::OrbitalSubspacePN& subspace = space.GetSubspace(subspace_index);

      std::cout << " index " << subspace_index
                << " species " << int(subspace.orbital_species())
                << std::endl;

      std::cout << subspace.DebugStr();
    }

  // check file output
  std::ofstream os(filename.c_str());
  os << space.OrbitalDefinitionStr();

}

void test_two_body_Nmax()
{

//  ////////////////////////////////////////////////////////////////
//  // two-body basis tests
//  ////////////////////////////////////////////////////////////////
//
//  std::cout << "Two-body basis" << std::endl;
//
//  // example subspace
//  std::cout << "Example subspace" << std::endl;
//  basis::TwoBodySubspaceJJJT subspace(0,1,0,6);
//  std::cout << subspace.DebugStr();
//
//
//  // spaces and sectors
//
//  // first set up space
//
//  std::cout << "Two-body space" << std::endl;
//  int Nmax = 2;
//  basis::TwoBodySpaceJJJT space(Nmax);
//  std::cout << space.DebugStr();
//
//  // then set up allowed sectors
//  std::cout << "Two-body operator sectors" << std::endl;
//  int J0 = 2;  // try: J0=0 for interaction, J0=2 for quadrupole operator
//  int T0 = 0;
//  int g0 = 0;
//  basis::TwoBodySectorsJJJT sectors(space,J0,T0,g0);
//
//  std::cout << " J0 " << J0 << " T0 " << T0 << " g0 " << g0 << std::endl;
//  for (int sector_index=0; sector_index < sectors.size(); ++sector_index)
//    {
//      int bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
//      const basis::TwoBodySubspaceJJJT& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
//      int ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
//      const basis::TwoBodySubspaceJJJT& ket_subspace = sectors.GetSector(sector_index).ket_subspace();
//
//      std::cout 
//        << " sector " 
//        << std::setw(3) << sector_index 
//        << "     "
//        << " index "
//        << std::setw(3) << bra_subspace_index
//        << " JTg "
//        << std::setw(3) << bra_subspace.J() 
//        << std::setw(3) << bra_subspace.T() 
//        << std::setw(3) << bra_subspace.g()
//        // << std::setw(3) << bra_subspace.Nmax()
//        << " dim "
//        << std::setw(3) << bra_subspace.size()
//        << "     "
//        << " index "
//        << std::setw(3) << ket_subspace_index
//        << " JTg "
//        << std::setw(3) << ket_subspace.J() 
//        << std::setw(3) << ket_subspace.T() 
//        << std::setw(3) << ket_subspace.g()
//        // << std::setw(3) << ket_subspace.Nmax()
//        << " dim "
//        << std::setw(3) << ket_subspace.size()
//        << std::endl;
//    }
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  std::string filename("test/jjjpn_scheme_general_test_orbitals_Nmax04.dat");
  test_orbitals_Nmax(filename);
  test_two_body_Nmax();

  // termination
  return 0;
}
