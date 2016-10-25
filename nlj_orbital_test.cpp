/****************************************************************
  jjjpnorb_scheme_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <iomanip>
#include <iostream>

#include "nlj_orbital.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void TestOrbitalsNmax(const std::string& filename)
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
  std::ofstream os(filename);
  os << basis::OrbitalDefinitionStr(space.OrbitalInfo(),true);

}

void TestOrbitalsRead(const std::string& filename) {
  std::cout << "Read Orbitals" << std::endl;

  std::ifstream is(filename);
  std::vector<basis::OrbitalPNInfo> states = basis::ParseOrbitalPNStream(is,true);
  basis::OrbitalSpacePN space(states);
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
}


void TestLJOrbitalsNmax(const std::string& filename)
{
  std::cout << "Orbitals -- Nmax scheme -- lj-subspaces" << std::endl;

  // set up space
  std::cout << "Space" << std::endl;
  int Nmax = 4;
  basis::OrbitalSpaceLJPN space(Nmax);
  std::cout << space.DebugStr();

  // check subspaces
  std::cout << "Subspaces" << std::endl;
  for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
    {
      const basis::OrbitalSubspaceLJPN& subspace = space.GetSubspace(subspace_index);

      std::cout << " index " << subspace_index
                << " species " << int(subspace.orbital_species())
                << " l " << int(subspace.l())
                << " j " << subspace.j().Str()
                << std::endl;

      std::cout << subspace.DebugStr();
    }

  // check file output
  std::ofstream os(filename);
  os << basis::OrbitalDefinitionStr(space.OrbitalInfo(),true);

}

void TestLJOrbitalsRead(const std::string& filename) {
  std::cout << "Read Orbitals -- lj-subspaces" << std::endl;

  std::ifstream is(filename);
  std::vector<basis::OrbitalPNInfo> states = basis::ParseOrbitalPNStream(is,true);
  basis::OrbitalSpaceLJPN space(states);
  std::cout << space.DebugStr();

  // check subspaces
  std::cout << "Subspaces" << std::endl;
  for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
    {
      const basis::OrbitalSubspaceLJPN& subspace = space.GetSubspace(subspace_index);

      std::cout << " index " << subspace_index
                << " species " << int(subspace.orbital_species())
                << " l " << int(subspace.l())
                << " j " << subspace.j().Str()
                << std::endl;

      std::cout << subspace.DebugStr();
    }
}

void TestLJSectors() {
  const int width = 3;
  std::cout << "Sectors -- Nmax scheme -- lj-subspaces" << std::endl;

  // set up space
  std::cout << "Space" << std::endl;
  int Nmax = 4;
  basis::OrbitalSpaceLJPN space(Nmax);
  std::cout << space.DebugStr();

  std::cout << "Sectors -- all-to-all" << std::endl;
  basis::OrbitalSectorsLJPN sectors(space);
  std::cout << sectors.DebugStr();

  std::cout << "Sectors -- l0max = 2" << std::endl;
  basis::OrbitalSectorsLJPN constrained_sectors(space, 2, 0);
  std::cout << constrained_sectors.DebugStr();

  std::cout << "Bra space" << std::endl;
  basis::OrbitalSpaceLJPN bra_space(Nmax+2);
  std::cout << bra_space.DebugStr();
  std::cout << "Ket space" << std::endl;
  basis::OrbitalSpaceLJPN ket_space(Nmax);
  std::cout << ket_space.DebugStr();

  std::cout << "Sectors -- all-to-all distinct spaces" << std::endl;
  basis::OrbitalSectorsLJPN distinct_space_sectors(bra_space, ket_space);
  std::cout << distinct_space_sectors.DebugStr();

  std::cout << "Sectors -- distinct spaces, l0max = 2" << std::endl;
  basis::OrbitalSectorsLJPN constrained_distinct_space_sectors(bra_space, ket_space, 2, 0);
  std::cout << constrained_distinct_space_sectors.DebugStr();
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  std::string filename("test/jjjpnorb_scheme_test_orbitals_Nmax04.dat");
  TestOrbitalsNmax(filename);
  TestOrbitalsRead(filename);
  std::string filename2("test/ljpn_scheme_test_orbitals_Nmax04.dat");
  TestLJOrbitalsNmax(filename2);
  TestLJOrbitalsRead(filename2);

  TestLJSectors();


  // termination
  return 0;
}
