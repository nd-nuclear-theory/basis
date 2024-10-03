/****************************************************************
  m_scheme_test.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <iostream>

#include "basis/m_scheme.h"
#include "basis/nlj_orbital.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void TestStatesNmax()
{

  ////////////////////////////////////////////////////////////////
  // orbital tests (should come after suborbital tests are passed)
  ////////////////////////////////////////////////////////////////

  std::cout << "Orbitals -- Nmax scheme" << std::endl;

  // set up space
  std::cout << "Space" << std::endl;
  int Nmax = 4;
  basis::OrbitalSpacePN orbital_space(Nmax);
  basis::SingleParticleSpacePN space(orbital_space);
  std::cout << space.DebugStr();

  // check subspaces
  std::cout << "Subspaces" << std::endl;
  for (std::size_t subspace_index=0; subspace_index<space.size(); ++subspace_index)
    {
      const basis::SingleParticleSubspacePN& subspace = space.GetSubspace(subspace_index);

      std::cout << " index " << subspace_index
                << " species " << int(subspace.orbital_species())
                << std::endl;

      std::cout << subspace.DebugStr();
    }

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  TestStatesNmax();

  // termination
  return 0;
}
