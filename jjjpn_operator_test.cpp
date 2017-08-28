/****************************************************************
  jjjpn_operator_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <iomanip>

#include "jjjpn_operator.h"


void WriteTest(const std::string& filename)
{
  ////////////////////////////////////////////////////////////////
  // defining operator
  ////////////////////////////////////////////////////////////////

  // since we set the operator to a "naive identity" operator, the
  // matrix elements are taken to be NAS

  // set up orbitals
  int orbital_Nmax = 4;
  basis::OrbitalSpacePN orbital_space(orbital_Nmax);

  // set up space
  int Nmax = 2;
  basis::TwoBodySpaceJJJPN space(
      orbital_space,
      basis::WeightMax(basis::Rank::kTwoBody,Nmax)
    );

  // set up operator containers
  basis::TwoBodySectorsJJJPN sectors;
  basis::OperatorBlocks<double> matrices;

  // populate operator containers
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  // enumerate sectors
  sectors = basis::TwoBodySectorsJJJPN(space,J0,g0,Tz0);
  basis::SetOperatorToIdentity(sectors,matrices);

  ////////////////////////////////////////////////////////////////
  // write test
  ////////////////////////////////////////////////////////////////

  // set up stream for output
  std::ostringstream os;

  // write matrices
  basis::WriteTwoBodyOperatorJJJPN(
      os,
      sectors,matrices,
      basis::NormalizationConversion::kNone,
      1
    );

  // dump to file
  std::ofstream ofile(filename.c_str());
  ofile << os.str();
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  std::string filename("test/jjjpnorb_operator_test_identity_Nmax02_NAS.dat");
  WriteTest(filename);

  // termination
  return 0;
}
