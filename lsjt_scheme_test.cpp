/****************************************************************
  lsjt_scheme_test.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/

#include <iomanip>
#include <iostream>

#include "lsjt_scheme.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // relative state tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {
      std::cout << "Relative subspace" << std::endl;

      // subspace construction

      // RelativeSubspaceLSJT(0,0,0,0,0,0);  // should violate assertion due to T
      // RelativeSubspaceLSJT(0,0,0,1,0,7);  // should violate assertion due to Nmax
      basis::RelativeSubspaceLSJT subspace(0,0,0,1,0,6);  // LSJTg Nmax


      // index-based looping
      for (int k=0; k<subspace.size(); ++k)
	{
	  basis::RelativeStateLSJT state(subspace,k);
	  std::cout << "index " << state.index() << " Nr " << state.Nr() << std::endl;
	};

      std::cout << "State construction by index" << std::endl;
      for (int Nr=0; Nr<=6; Nr+=2)
	{
	  basis::RelativeStateLSJT state(subspace,Nr/2);
	  std::cout << "Nr " << Nr << " index " << state.index() << " Nr " << state.Nr() << std::endl;
	}

      std::cout << "State lookup" << std::endl;
      for (int Nr=0; Nr<=6; Nr+=2)
	{
	  basis::RelativeStateLSJT state(subspace,basis::RelativeSubspaceLSJT::StateLabelsType(Nr));
	  std::cout << "Nr " << Nr << " index " << state.index() << " Nr " << state.Nr() << std::endl;
	}


      // spaces and sectors

      // first set up relative spaces

      std::cout << "Relative space" << std::endl;
      int Nmax = 2;
      basis::RelativeSpaceLSJT space(Nmax);
      std::cout << space.Str();

      //      for (int s=0; s<space.size(); ++s)
      //	{
      //	  const RelativeSubspaceLSJT& subspace = space.GetSubspace(s);
      //	  std::cout 
      //	    << std::setw(3) << s 
      //	    << std::setw(3) << subspace.L() 
      //	    << std::setw(3) << subspace.S() 
      //	    << std::setw(3) << subspace.J() 
      //	    << std::setw(3) << subspace.T() 
      //	    << std::setw(3) << subspace.g()
      //	    << std::setw(3) << subspace.Nmax()
      //	    << std::setw(3) << subspace.size()
      //	    << std::endl;
      //	}

      // then set up allowed sectors
      std::cout << "Relative interaction sectors" << std::endl;
      int J0 = 0;  // also try J0=2 for quadrupole operator
      int g0 = 0;
      basis::RelativeSectorsLSJT sectors(space,J0,g0);

      for (int sector_index=0; sector_index < sectors.size(); ++sector_index)
	{
	  int bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
	  const basis::RelativeSubspaceLSJT& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
	  int ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
	  const basis::RelativeSubspaceLSJT& ket_subspace = sectors.GetSector(sector_index).ket_subspace();

	  std::cout 
	    << " sector " 
	    << std::setw(3) << sector_index 
	    << "     "
	    << " index "
	    << std::setw(3) << bra_subspace_index
	    << " LSJTg "
	    << std::setw(3) << bra_subspace.L() 
	    << std::setw(3) << bra_subspace.S() 
	    << std::setw(3) << bra_subspace.J() 
	    << std::setw(3) << bra_subspace.T() 
	    << std::setw(3) << bra_subspace.g()
	    // << std::setw(3) << bra_subspace.Nmax()
	    << " dim "
	    << std::setw(3) << bra_subspace.size()
	    << "     "
	    << " index "
	    << std::setw(3) << ket_subspace_index
	    << " LSJTg "
	    << std::setw(3) << ket_subspace.L() 
	    << std::setw(3) << ket_subspace.S() 
	    << std::setw(3) << ket_subspace.J() 
	    << std::setw(3) << ket_subspace.T() 
	    << std::setw(3) << ket_subspace.g()
	    // << std::setw(3) << ket_subspace.Nmax()
	    << " dim "
	    << std::setw(3) << ket_subspace.size()
	    << std::endl;
	}

    }

  ////////////////////////////////////////////////////////////////
  // relative-cm state tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {
      std::cout << "Spurious relative state" << std::endl;

      // subspace construction

      // "bad" case to test validity assertion
      if (false)
        {
          basis::RelativeCMSubspaceLSJT bad_subspace(
              1,1,  // Ncm, lcm
              2,1,3,0,0,   // L, S, J, T, g
              6  // Nmax
            );
        }

      basis::RelativeCMSubspaceLSJT subspace(
          1,1,  // Ncm, lcm
          2,1,3,1,0,   // L, S, J, T, g
          6  // Nmax
        );

      // index-based looping
      const int dimension = subspace.size();
      for (int k = 0; k < dimension; ++k)
	{
	  basis::RelativeCMStateLSJT state(subspace,k);
	  std::cout << state.index() 
		    << " " << state.N() 
		    << " " << state.Nr() << " " << state.lr() 
		    << std::endl;
	};
    }

  ////////////////////////////////////////////////////////////////
  // two-body state tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {
      std::cout << "Two-body subspace" << std::endl;

      // subspace construction

      basis::TwoBodySubspaceLSJT subspace(0,0,0,1,0,6);

      // index-based looping
      for (int k = 0; k < subspace.size(); ++k)
	{
	  basis::TwoBodyStateLSJT state(subspace,k);
	  std::cout << state.index() 
		    << " " << state.N() 
		    << " " << state.N1() << " " << state.l1() << " " << state.N2() << " " << state.l2() 
		    << std::endl;
	};


      // spaces and sectors

      // first set up relative spaces

      std::cout << "Two-body space" << std::endl;
      int Nmax = 3;
      basis::TwoBodySpaceLSJT space(Nmax);
      std::cout << space.Str();

      //      for (int s=0; s<space.size(); ++s)
      //	{
      //	  const TwoBodySubspaceLSJT& subspace = space.GetSubspace(s);
      //	  std::cout 
      //	    << std::setw(3) << s 
      //	    << std::setw(3) << subspace.L() 
      //	    << std::setw(3) << subspace.S() 
      //	    << std::setw(3) << subspace.J() 
      //	    << std::setw(3) << subspace.T() 
      //	    << std::setw(3) << subspace.g()
      //	    << std::setw(3) << subspace.Nmax()
      //	    << std::setw(3) << subspace.size()
      //	    << std::endl;
      //	}

      // then set up allowed sectors
      std::cout << "Two-body interaction sectors" << std::endl;
      int J0 = 0;  // also try J0=2 for quadrupole operator
      int g0 = 0;
      basis::TwoBodySectorsLSJT sectors(space,J0,g0);

      for (int sector_index=0; sector_index < sectors.size(); ++sector_index)
	{
	  int bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
	  const basis::TwoBodySubspaceLSJT& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
	  int ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
	  const basis::TwoBodySubspaceLSJT& ket_subspace = sectors.GetSector(sector_index).ket_subspace();

	  std::cout 
	    << " sector " 
	    << std::setw(3) << sector_index 
	    << "     "
	    << " index "
	    << std::setw(3) << bra_subspace_index
	    << " LSJTg "
	    << std::setw(3) << bra_subspace.L() 
	    << std::setw(3) << bra_subspace.S() 
	    << std::setw(3) << bra_subspace.J() 
	    << std::setw(3) << bra_subspace.T() 
	    << std::setw(3) << bra_subspace.g()
	    // << std::setw(3) << bra_subspace.Nmax()
	    << " dim "
	    << std::setw(3) << bra_subspace.size()
	    << "     "
	    << " index "
	    << std::setw(3) << ket_subspace_index
	    << " LSJTg "
	    << std::setw(3) << ket_subspace.L() 
	    << std::setw(3) << ket_subspace.S() 
	    << std::setw(3) << ket_subspace.J() 
	    << std::setw(3) << ket_subspace.T() 
	    << std::setw(3) << ket_subspace.g()
	    // << std::setw(3) << ket_subspace.Nmax()
	    << " dim "
	    << std::setw(3) << ket_subspace.size()
	    << std::endl;
	}
    }

  
  // termination
  return 0;
}
