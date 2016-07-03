/****************************************************************
  lsjt_operator_test.cpp

  Mark A. Caprio, University of Notre Dame.
  See interaction_lsjt.h for documentation.

****************************************************************/

#include <iomanip>

#include "lsjt_operator.h"


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // relative basis tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {

      ////////////////////////////////////////////////////////////////
      // defining operator
      ////////////////////////////////////////////////////////////////

      std::cout << "Setup" << std::endl;

      // set up space
      int Nmax = 2;
      int Jmax = 3;
      basis::RelativeSpaceLSJT space(Nmax,Jmax);

      // set up operator containers
      //
      // These are vectors to store information for T0=0/1/2 components.
      std::vector<basis::RelativeSectorsLSJT> component_sectors(3);
      std::vector<basis::MatrixVector> component_matrices(3);

      // populate operator containers
      int J0 = 0;
      int g0 = 0;
      for (int T0=0; T0<=2; ++T0)
        // for each isospin component
        {

          // enumerate sectors
          component_sectors[T0] = basis::RelativeSectorsLSJT(space,J0,T0,g0);
          std::cout << " T0 " << T0 << " size " << component_sectors[T0].size() << std::endl;
          
          // populate matrices
          if (T0==0)
            basis::SetOperatorToIdentity(component_sectors[T0],component_matrices[T0]);
          else
            basis::SetOperatorToZero(component_sectors[T0],component_matrices[T0]);
        }

      ////////////////////////////////////////////////////////////////
      // write test
      ////////////////////////////////////////////////////////////////

      // set up stream for output
      //std::ofstream os("lsjt_operator_test_identity.dat");
      std::ostringstream os;

      // write header parameters
      basis::RelativeOperatorParametersLSJT parameters;
      parameters.version=1;
      parameters.J0=J0;
      parameters.g0=g0;
      parameters.T0_min=0;
      parameters.T0_max=2;
      parameters.Nmax=Nmax;
      parameters.Jmax=Jmax;
      basis::WriteRelativeOperatorParametersLSJT(os,parameters);

      // write matrices
      for (int T0=parameters.T0_min; T0<=parameters.T0_max; ++T0)
        {
          basis::WriteRelativeOperatorComponentLSJT(
              os,
              T0,
              component_sectors[T0],component_matrices[T0]
            );
        }

      // dump to terminal for inspection
      std::cout << os.str();

      ////////////////////////////////////////////////////////////////
      // read test
      ////////////////////////////////////////////////////////////////

      std::cout << "Readback test" << std::endl;

      // set up stream for readback
      std::istringstream is(os.str());

      // read header parameters
      basis::RelativeOperatorParametersLSJT parameters2;
      basis::ReadRelativeOperatorParametersLSJT(is,parameters2);
      // and inspect...
      basis::WriteRelativeOperatorParametersLSJT(std::cout,parameters2);

      // set up space
      basis::RelativeSpaceLSJT space2(parameters2.Nmax,parameters2.Jmax);

      // read matrices
      std::vector<basis::RelativeSectorsLSJT> component_sectors2;
      component_sectors2.resize(3);
      std::vector<basis::MatrixVector> component_matrices2;
      component_matrices2.resize(3);

      for (int T0=parameters2.T0_min; T0<=parameters2.T0_max; ++T0)
        // for each isospin component
        {
          // enumerate sectors
          component_sectors2[T0] = basis::RelativeSectorsLSJT(space2,J0,T0,g0);

          // read matrices
          basis::ReadRelativeOperatorComponentLSJT(
              is,
              T0,
              component_sectors2[T0],component_matrices2[T0]
            );
          // and inspect...
          basis::WriteRelativeOperatorComponentLSJT(
              std::cout,
              T0,
              component_sectors2[T0],component_matrices2[T0]
            );
        }
    }

  // termination
  return 0;
}
