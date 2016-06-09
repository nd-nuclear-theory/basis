/****************************************************************
  lsjt_operator_test.cpp

  Mark A. Caprio, University of Notre Dame.
  See interaction_lsjt.h for documentation.

****************************************************************/

#include <iomanip>
#include <sstream>

#include "lsjt_operator.h"

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

      std::cout << "Relative space" << std::endl;
      int Nmax = 2;
      basis::RelativeSpaceLSJT space(Nmax);
      int J0 = 0;  // also try J0=2 for quadrupole operator
      int g0 = 0;
      basis::RelativeSectorsLSJT sectors(space,J0,g0);
      basis::MatrixVector matrices;

      std::cout << "Zeros" << std::endl;
      basis::SetOperatorToZero(sectors,matrices);
      basis::WriteRelativeOperatorLSJT(std::cout,sectors,matrices);

      std::cout << "Identity" << std::endl;
      basis::SetOperatorToNaiveIdentityReduced(sectors,matrices);
      basis::WriteRelativeOperatorLSJT(std::cout,sectors,matrices);
      // std::cout << "...printed as matrices..." << std::endl;
      // WriteRelativeOperatorMatrices(std::cout,sectors,matrices,6,3);

      // writeout/readback test
      std::cout << "Readback" << std::endl;

      std::ostringstream os;
      basis::WriteRelativeOperatorLSJT(os,sectors,matrices);
      // std::cout << os.str() << std::endl; // debugging: inspect stream contents

      std::istringstream is(os.str());
      basis::MatrixVector matrices2;
      basis::ReadRelativeOperatorLSJT(is,sectors,matrices2);
      basis::WriteRelativeOperatorLSJT(std::cout,sectors,matrices2);

    }

  ////////////////////////////////////////////////////////////////
  // two-body state tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {

      std::cout << "Two-body space" << std::endl;
      int Nmax = 2;
      basis::TwoBodySpaceLSJT space(Nmax);
      int J0 = 0;  // also try J0=2 for quadrupole operator
      int g0 = 0;
      basis::TwoBodySectorsLSJT sectors(space,J0,g0);
      basis::MatrixVector matrices;

      std::cout << "Identity" << std::endl;
      basis::SetOperatorToNaiveIdentityReduced(sectors,matrices);
      basis::WriteTwoBodyOperatorLSJT(std::cout,sectors,matrices);
      //std::cout << "...printed as matrices..." << std::endl;
      //basis::WriteTwoBodyOperatorMatricesLSJT(std::cout,sectors,matrices,9,6);

    }


  // termination
  return 0;
}
