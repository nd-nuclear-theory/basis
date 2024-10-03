/****************************************************************
  jjjpn_operator_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <string>

#include "basis/jjjpn_operator.h"
#include "basis/jjjpn_scheme.h"
#include "basis/many_body.h"
#include "basis/nlj_orbital.h"
#include "basis/operator.h"


void CanonicalizationTest(const std::string& filename)
{
  ////////////////////////////////////////////////////////////////
  // set up operator storage
  ////////////////////////////////////////////////////////////////

  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  const basis::TwoBodySpaceJJJPN two_body_space(
      orbital_space,
      weight_max,
      shell::kH2SpaceOrdering.at(run_parameters.output_format)
    );
  const basis::TwoBodySectorsJJJPN two_body_sectors(
      two_body_space,
      J0, g0, Tz0
    );
  basis::OperatorBlocks<double> two_body_matrices;
  basis::SetOperatorToZero(two_body_sectors, two_body_matrices);

  ////////////////////////////////////////////////////////////////
  // test lookups
  ////////////////////////////////////////////////////////////////

  // lookup tests
  LookUpStateFromXPNLabels(orbital_space, two_body_space, 1, 1, 0);
  LookUpStateFromXPNLabels(orbital_space, two_body_space, 2, 1, 0);
  LookUpStateFromXPNLabels(orbital_space, two_body_space, 1, 29, 0);
  // LookUpStateFromXPNLabels(orbital_space, two_body_space, 29, 1, 0);  // FAILS -- ILLEGAL
}


void WriteTest(const std::string& filename)
{
  ////////////////////////////////////////////////////////////////
  // define operator
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
  // test write
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

  std::string filename("jjjpnorb_operator_test_identity_Nmax02_NAS.dat");
  CanonicalizationTest(filename);

  WriteTest(filename);

  // termination
  return 0;
}
