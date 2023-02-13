/****************************************************************
  oscillator_orbital_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <iomanip>
#include <iostream>
#include <limits>

#include "oscillator_orbital.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void TestOscillator()
{
  ////////////////////////////////////////////////////////////////
  // relative basis tests
  ////////////////////////////////////////////////////////////////

  std::cout << "Oscillator orbital basis" << std::endl;

  // subspace construction

  // Normally subspaces are always constructed as part of a space.  We would not
  // construct a standalone subspace.  But, for testing purposes, it is
  // convenient to do so, since the code for subspace and state can be written
  // and tested independently of the code for space and/or sectors.

  const basis::OscillatorOrbitalSubspace subspace(2, 4);  // l=2, Nmax=4
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();

  // state construction

  // Let us give the state type a more thorough workout.  Note that the state
  // type was already used in the implementation of the subspace's DebugStr()
  // above.

  // illustrate construction by index vs. state labels

  // construct by index -- invalid cases
  
  if (false)
    {
      // Index -1 should be trapped and lead to assertion failure.
      //
      // Note: Negative index cast to size_t will be large positive index, by
      // two-complement.

      std::size_t test_index = (std::size_t)(-1);
      std::cout << "Constructing with index " << test_index << std::endl;
      const basis::OscillatorOrbitalState state_from_index(subspace, test_index);
    }

  if (false)
    {
      // But a very large negative integer could wrap and become a valid index again!
      //
      // std::size_t test_index = (std::size_t)(-18446744073709551615);
      
      std::size_t max_size = std::numeric_limits<std::size_t>().max();
      std::size_t test_index = -max_size;  // should be equivalent to +1
      std::cout << "Constructing with index " << test_index << std::endl;
      const basis::OscillatorOrbitalState state_from_index(subspace, test_index);
    }

  // construct by index -- valid case
  std::size_t test_index = 1;
  std::cout << "Constructing with index " << test_index << std::endl;
  const basis::OscillatorOrbitalState state_from_index(subspace, test_index);
  std::cout << "  label string " << state_from_index.LabelStr() << std::endl;
  std::cout << "  quantum numbers"
            << " l " << state_from_index.l()
            << " g " << state_from_index.g()
            << " n " << state_from_index.n()
            << " N " << state_from_index.N()
            << std::endl;

  // construct by state labels
  int n = 0;
  std::cout << "Constructing with n " << n << std::endl;
  const basis::OscillatorOrbitalState state_from_labels(subspace, basis::OscillatorOrbitalState::StateLabelsType(n));
  std::cout << "  label string " << state_from_labels.LabelStr() << std::endl;
  std::cout << "  quantum numbers"
            << " l " << state_from_labels.l()
            << " g " << state_from_labels.g()
            << " n " << state_from_labels.n()
            << " N " << state_from_labels.N()
            << std::endl;
  
  // text equality operator
  std::cout << "  equality test " << (state_from_index == state_from_index) << std::endl;
  std::cout << "  equality test " << (state_from_index == state_from_labels) << std::endl;
  

  //  index-based looping
  // r (std::size_t k=0; k<subspace.size(); ++k)
  // {
  //   basis::RelativeStateLSJT state(subspace,k);
  //   std::cout << "index " << state.index() << " N " << state.N() << std::endl;
  // };
  // 
  // d::cout << "State construction by index" << std::endl;
  // r (int N=0; N<=6; N+=2)
  // {
  //   basis::RelativeStateLSJT state(subspace,N/2);
  //   std::cout << "N " << N << " index " << state.index() << " N " << state.N() << " n " << state.n() << std::endl;
  // }
  // 
  // d::cout << "State lookup" << std::endl;
  // r (int N=0; N<=6; N+=2)
  // {
  //   basis::RelativeStateLSJT state(subspace,basis::RelativeSubspaceLSJT::StateLabelsType(N));
  //   std::cout << "N " << N << " index " << state.index() << " N " << state.N() << " n " << state.n() << std::endl;
  // }
  // 
  // 
  //  spaces and sectors
  // 
  //  first set up space
  // 
  // d::cout << "Relative space" << std::endl;
  // t N_max = 2;
  // t J_max = 3;
  // sis::RelativeSpaceLSJT space(N_max,J_max);
  // d::cout << space.DebugStr();
  // 
  //  then set up allowed sectors
  // d::cout << "Relative operator sectors" << std::endl;
  // t J0 = 0;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  // t T0 = 1;
  // t g0 = 0;
  // sis::RelativeSectorsLSJT sectors(space,J0,T0,g0);
  // 
  // 
  // d::cout << " J0 " << J0 << " T0 " << T0 << " g0 " << g0 << std::endl;
  // r (std::size_t sector_index=0; sector_index < sectors.size(); ++sector_index)
  // {
  //   std::size_t bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
  //   const basis::RelativeSubspaceLSJT& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
  //   std::size_t ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
  //   const basis::RelativeSubspaceLSJT& ket_subspace = sectors.GetSector(sector_index).ket_subspace();
  // 
  //   std::cout
  //     << " sector "
  //     << std::setw(3) << sector_index
  //     << "     "
  //     << " index "
  //     << std::setw(3) << bra_subspace_index
  //     << " LSJTg "
  //     << std::setw(3) << bra_subspace.L()
  //     << std::setw(3) << bra_subspace.S()
  //     << std::setw(3) << bra_subspace.J()
  //     << std::setw(3) << bra_subspace.T()
  //     << std::setw(3) << bra_subspace.g()
  //     << " dim "
  //     << std::setw(3) << bra_subspace.size()
  //     << "     "
  //     << " index "
  //     << std::setw(3) << ket_subspace_index
  //     << " LSJTg "
  //     << std::setw(3) << ket_subspace.L()
  //     << std::setw(3) << ket_subspace.S()
  //     << std::setw(3) << ket_subspace.J()
  //     << std::setw(3) << ket_subspace.T()
  //     << std::setw(3) << ket_subspace.g()
  //     << " dim "
  //     << std::setw(3) << ket_subspace.size()
  //     << " "
  //     << (sectors.GetSector(sector_index).IsDiagonal() ? "*" : " ")
  //     << std::endl;
  // }

}


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  TestOscillator();

  // termination
  return 0;
}
