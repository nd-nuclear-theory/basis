/****************************************************************
  jjjpn_operator_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <string>

#include "basis/jjjttz_operator.h"
#include "basis/jjjttz_scheme.h"
#include "basis/many_body.h"
#include "basis/nlj_orbital.h"
#include "basis/operator.h"


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
  basis::TwoBodySpaceJJJTTz space(
      basis::Rank::kTwoBody,
      Nmax
    );

  // set up operator containers
  basis::TwoBodySectorsJJJTTz sectors;
  basis::OperatorBlocks<double> matrices;

  // populate operator containers
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  // enumerate sectors
  sectors = basis::TwoBodySectorsJJJTTz(space,J0,g0,Tz0);
  basis::SetOperatorToIdentity(sectors,matrices);

  ////////////////////////////////////////////////////////////////
  // write test
  ////////////////////////////////////////////////////////////////

  // set up stream for output
  std::ostringstream os;

  // write matrices
  basis::WriteTwoBodyOperatorJJJTTz(
      os,
      sectors,matrices,
      basis::NormalizationConversion::kNone
    );

  // dump to file
  std::ofstream ofile(filename.c_str());
  ofile << os.str();

  // test finding and saving matrix element
  basis::SetTwoBodyOperatorMatrixElementJJJTTz(
      space,
      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(1,0,1,0),
      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(1,0,1,0),
      basis::TwoBodyStateJJJTTz::StateLabelsType(0,HalfInt(1,2),1,HalfInt(1,2)),
      basis::TwoBodyStateJJJTTz::StateLabelsType(0,HalfInt(1,2),1,HalfInt(1,2)),
      sectors,
      matrices,
      100
    );
}

void StreamIOTest(const std::string& filenamei, const std::string& filenameo)
{
  std::ios_base::openmode mode_argumenti = std::ios_base::in;
  std::ifstream is(filenamei.c_str(), mode_argumenti);
  if (!is) {
    std::cout << "can't find file" << std::endl;
  }
  float a,b,c;
  is >> a;
  if (!is) {
    std::cout << "can't read number" << std::endl;
  }
  is >> b;
  if (!is) {
    std::cout << "can't read number" << std::endl;
  }
  is >> c;
  if (!is) {
    std::cout << "can't read number" << std::endl;
  }
  std::cout << a << std::endl << b << std::endl << c << std::endl;

  std::ios_base::openmode mode_argumento = std::ios_base::out;
  std::ofstream os(filenameo.c_str(), mode_argumento);
  os << a << " " << b << " " << c << std::endl;

  std::cout << int(HalfInt(1,2)+HalfInt(1,2)) << std::endl;
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  std::string filename("test/jjjttzorb_operator_test_identity_Nmax02.dat");
  WriteTest(filename);
  std::string filenamei("test/test.txt");
  std::string filenameo("test/outtest.txt");
  StreamIOTest(filenamei,filenameo);

  // termination
  return 0;
}
