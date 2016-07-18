/****************************************************************
  lsjt_scheme_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <iomanip>
#include <iostream>

#include "lsjt_scheme.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void TestRelative()
{
  ////////////////////////////////////////////////////////////////
  // relative basis tests
  ////////////////////////////////////////////////////////////////

  std::cout << "Relative basis" << std::endl;

  // subspace construction

  // RelativeSubspaceLSJT(0,0,0,0,0,0);  // should violate assertion due to T
  // RelativeSubspaceLSJT(0,0,0,1,0,7);  // should violate assertion due to Nmax
  basis::RelativeSubspaceLSJT subspace(0,0,0,1,0,6);  // LSJTg N_max
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();


  // index-based looping
  for (int k=0; k<subspace.size(); ++k)
    {
      basis::RelativeStateLSJT state(subspace,k);
      std::cout << "index " << state.index() << " N " << state.N() << std::endl;
    };

  std::cout << "State construction by index" << std::endl;
  for (int N=0; N<=6; N+=2)
    {
      basis::RelativeStateLSJT state(subspace,N/2);
      std::cout << "N " << N << " index " << state.index() << " N " << state.N() << std::endl;
    }

  std::cout << "State lookup" << std::endl;
  for (int N=0; N<=6; N+=2)
    {
      basis::RelativeStateLSJT state(subspace,basis::RelativeSubspaceLSJT::StateLabelsType(N));
      std::cout << "N " << N << " index " << state.index() << " N " << state.N() << std::endl;
    }


  // spaces and sectors

  // first set up space

  std::cout << "Relative space" << std::endl;
  int N_max = 2;
  int J_max = 3;
  basis::RelativeSpaceLSJT space(N_max,J_max);
  std::cout << space.DebugStr();

  // then set up allowed sectors
  std::cout << "Relative operator sectors" << std::endl;
  int J0 = 0;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  int T0 = 1;
  int g0 = 0;
  basis::RelativeSectorsLSJT sectors(space,J0,T0,g0);


  std::cout << " J0 " << J0 << " T0 " << T0 << " g0 " << g0 << std::endl;
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
        << " dim "
        << std::setw(3) << ket_subspace.size()
        << " "
        << (sectors.GetSector(sector_index).IsDiagonal() ? "*" : " ")
        << std::endl;
    }

}

void TestRelativeCM()
{
  ////////////////////////////////////////////////////////////////
  // relative-cm basis tests
  ////////////////////////////////////////////////////////////////

  std::cout << "Relative-c.m. basis" << std::endl;

  // example subspace
  std::cout << "Example subspace" << std::endl;
  basis::RelativeCMSubspaceLSJT subspace(0,0,0,0,0,2);  // LSJTg Nmax
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();

  // spaces and sectors

  // first set up space

  std::cout << "Relative-c.m. space" << std::endl;
  int Nmax = 2;
  basis::RelativeCMSpaceLSJT space(Nmax);
  std::cout << space.DebugStr();

  // then set up allowed sectors
  std::cout << "Relative-c.m. operator sectors" << std::endl;
  int J0 = 2;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  int T0 = 0;
  int g0 = 0;
  basis::RelativeCMSectorsLSJT sectors(space,J0,T0,g0);

  std::cout << " J0 " << J0 << " T0 " << T0 << " g0 " << g0 << std::endl;
  for (int sector_index=0; sector_index < sectors.size(); ++sector_index)
    {
      int bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
      const basis::RelativeCMSubspaceLSJT& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
      int ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
      const basis::RelativeCMSubspaceLSJT& ket_subspace = sectors.GetSector(sector_index).ket_subspace();

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

void TestRelativeCMN()
{
  ////////////////////////////////////////////////////////////////
  // relative-cm basis tests -- subspaced by N
  ////////////////////////////////////////////////////////////////

  std::cout << "Relative-c.m. basis -- subspaced by N" << std::endl;

  // spaces and sectors

  // first set up space

  std::cout << "Relative-c.m. space" << std::endl;
  int Nmax = 2;
  basis::RelativeCMSpaceLSJTN space(Nmax);
  std::cout << space.DebugStr();

  // then set up allowed sectors
  std::cout << "Relative-c.m. operator sectors" << std::endl;
  int J0 = 2;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  int T0 = 0;
  int g0 = 0;
  basis::RelativeCMSectorsLSJTN sectors(space,J0,T0,g0);

  std::cout << " J0 " << J0 << " T0 " << T0 << " g0 " << g0 << std::endl;
  for (int sector_index=0; sector_index < sectors.size(); ++sector_index)
    {
      int bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
      const basis::RelativeCMSubspaceLSJTN& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
      int ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
      const basis::RelativeCMSubspaceLSJTN& ket_subspace = sectors.GetSector(sector_index).ket_subspace();

      std::cout 
        << " sector " 
        << std::setw(3) << sector_index 
        << "     "
        << " index "
        << std::setw(3) << bra_subspace_index
        << " LSJTg-N "
        << std::setw(3) << bra_subspace.L() 
        << std::setw(3) << bra_subspace.S() 
        << std::setw(3) << bra_subspace.J() 
        << std::setw(3) << bra_subspace.T() 
        << std::setw(3) << bra_subspace.g()
        << std::setw(3) << bra_subspace.N()
        << " dim "
        << std::setw(3) << bra_subspace.size()
        << "     "
        << " index "
        << std::setw(3) << ket_subspace_index
        << " LSJTg-N "
        << std::setw(3) << ket_subspace.L() 
        << std::setw(3) << ket_subspace.S() 
        << std::setw(3) << ket_subspace.J() 
        << std::setw(3) << ket_subspace.T() 
        << std::setw(3) << ket_subspace.g()
        << std::setw(3) << ket_subspace.N()
        << " dim "
        << std::setw(3) << ket_subspace.size()
        << std::endl;
    }
}

void TestTwoBody()
{
  ////////////////////////////////////////////////////////////////
  // two-body basis tests
  ////////////////////////////////////////////////////////////////

  std::cout << "Two-body basis" << std::endl;

  // example subspace
  std::cout << "Example subspace" << std::endl;
  basis::TwoBodySubspaceLSJT subspace(0,0,0,0,0,2);  // LSJTg Nmax
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();


  // spaces and sectors

  // first set up space

  std::cout << "Two-body space" << std::endl;
  int Nmax = 2;
  basis::TwoBodySpaceLSJT space(Nmax);
  std::cout << space.DebugStr();

  // then set up allowed sectors
  std::cout << "Two-body operator sectors" << std::endl;
  int J0 = 2;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  int T0 = 0;
  int g0 = 0;
  basis::TwoBodySectorsLSJT sectors(space,J0,T0,g0);

  std::cout << " J0 " << J0 << " T0 " << T0 << " g0 " << g0 << std::endl;
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

void TestTwoBodyN1max()
{
  ////////////////////////////////////////////////////////////////
  // two-body basis tests
  ////////////////////////////////////////////////////////////////

  std::cout << "Two-body basis (N1max)" << std::endl;

  // example subspace
  std::cout << "Example subspace" << std::endl;
  basis::TwoBodySubspaceLSJT subspace(0,0,0,0,0,4,1);  // LSJTg N1max rank=1
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();


  // spaces and sectors

  // first set up space

  std::cout << "Two-body space" << std::endl;
  int N1max = 4;
  basis::TwoBodySpaceLSJT space(N1max,1);
  std::cout << space.DebugStr();

  // then set up allowed sectors
  std::cout << "Two-body operator sectors" << std::endl;
  int J0 = 0;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  int T0 = 0;
  int g0 = 0;
  basis::TwoBodySectorsLSJT sectors(space,J0,T0,g0);

  std::cout << " J0 " << J0 << " T0 " << T0 << " g0 " << g0 << std::endl;
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

void TestTwoBodyN()
{
  ////////////////////////////////////////////////////////////////
  // two-body basis tests -- subspaced by N
  ////////////////////////////////////////////////////////////////

  std::cout << "Two-body basis -- subspaced by N" << std::endl;

  // example subspace
  std::cout << "Example subspace" << std::endl;
  basis::TwoBodySubspaceLSJTN subspace(0,0,0,0,0,2,2);  // LSJTg Nmax
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();

  // spaces and sectors

  // first set up space

  std::cout << "Two-body space" << std::endl;
  int Nmax = 2;
  basis::TwoBodySpaceLSJTN space(Nmax);
  std::cout << space.DebugStr();

  // then set up allowed sectors
  std::cout << "Two-body operator sectors" << std::endl;
  int J0 = 2;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  int T0 = 0;
  int g0 = 0;
  basis::TwoBodySectorsLSJTN sectors(space,J0,T0,g0);

  std::cout << " J0 " << J0 << " T0 " << T0 << " g0 " << g0 << std::endl;
  for (int sector_index=0; sector_index < sectors.size(); ++sector_index)
    {
      int bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
      const basis::TwoBodySubspaceLSJTN& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
      int ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
      const basis::TwoBodySubspaceLSJTN& ket_subspace = sectors.GetSector(sector_index).ket_subspace();

      std::cout 
        << " sector " 
        << std::setw(3) << sector_index 
        << "     "
        << " index "
        << std::setw(3) << bra_subspace_index
        << " LSJTg-N "
        << std::setw(3) << bra_subspace.L() 
        << std::setw(3) << bra_subspace.S() 
        << std::setw(3) << bra_subspace.J() 
        << std::setw(3) << bra_subspace.T() 
        << std::setw(3) << bra_subspace.g()
        << std::setw(3) << bra_subspace.N()
        << " dim "
        << std::setw(3) << bra_subspace.size()
        << "     "
        << " index "
        << std::setw(3) << ket_subspace_index
        << " LSJTg-N "
        << std::setw(3) << ket_subspace.L() 
        << std::setw(3) << ket_subspace.S() 
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

  TestRelative();

  TestRelativeCM();

  TestRelativeCMN();

  TestTwoBody();
  TestTwoBodyN1max();
  TestTwoBodyN();

  // termination
  return 0;
}
