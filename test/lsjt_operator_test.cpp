/****************************************************************
  lsjt_operator_test.cpp

  Mark A. Caprio, Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <iostream>

#include "basis/lsjt_operator.h"
#include "basis/jt_operator.h"
#include "basis/lsjt_scheme.h"
#include "basis/many_body.h"
#include "basis/operator.h"

////////////////////////////////////////////////////////////////
// test parts
////////////////////////////////////////////////////////////////

void WriteTestRelativeManual(const std::string& filename)
// "Manual" test of writing...  Supplanted by use of
// WriteRelativeOperatorLSJT().
{
  ////////////////////////////////////////////////////////////////
  // defining operator
  ////////////////////////////////////////////////////////////////

  std::cout << "Setup" << std::endl;

  // set operator and file header parameters
  basis::RelativeOperatorParametersLSJT operator_parameters;
  // operator tensorial parameters
  operator_parameters.J0=0;
  operator_parameters.g0=0;
  operator_parameters.symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
  operator_parameters.T0_min=0;
  operator_parameters.T0_max=2;
  // relative basis parameters
  operator_parameters.Nmax=2;
  operator_parameters.Jmax=operator_parameters.Nmax+1;

  // set up relative space
  basis::RelativeSpaceLSJT space(operator_parameters.Nmax,operator_parameters.Jmax);

  // set up operator containers
  //
  // These are arrays to store information for T0=0/1/2 components.
  std::array<basis::RelativeSectorsLSJT,3> component_sectors;
  std::array<basis::OperatorBlocks<double>,3> component_matrices;

  // populate operator containers
  //
  // Note: This can now instead be done automatically using
  // ConstructIdentityOperatorRelativeLSJT.

  for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      component_sectors[T0] = basis::RelativeSectorsLSJT(space,operator_parameters.J0,T0,operator_parameters.g0);
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
  std::ostringstream os;

  // write header parameters
  basis::WriteRelativeOperatorParametersLSJT(os,operator_parameters);

  // write matrices
  for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
    {
      basis::WriteRelativeOperatorComponentLSJT(
          os,
          T0,
          component_sectors[T0],component_matrices[T0]
        );
    }

  // dump to terminal for inspection
  std::cout << os.str();

  // dump to file
  std::ofstream ofile(filename.c_str());
  ofile << os.str();
}

void WriteTestRelative(const std::string& filename)
// Write using WriteRelativeOperatorLSJT().
//
// Actual contents of binary write:
//
//     % od -t x1 lsjt_operator_test_relative_identity_Nmax02.bin
//
//     0000000 01 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
//     0000020 02 00 00 00 00 00 00 00 02 00 00 00 03 00 00 00
//     0000040 00 00 00 00 00 00 f0 3f 00 00 00 00 00 00 00 00
//     0000060 00 00 00 00 00 00 f0 3f 00 00 00 00 00 00 f0 3f
//     0000100 00 00 00 00 00 00 00 00 00 00 00 00 00 00 f0 3f
//     0000120 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
//     0000140 00 00 00 00 00 00 f0 3f 00 00 00 00 00 00 f0 3f
//     *
//     0000240 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
//     *
//     0000440
//
//
//     % od -t u4 lsjt_operator_test_relative_identity_Nmax02.bin
//
//     0000000          1          0          0          0
//     0000020          2          0          2          3
//     0000040          0 1072693248          0          0
//     0000060          0 1072693248          0 1072693248
//     0000100          0          0          0 1072693248
//     0000120          0          0          0          0
//     0000140          0 1072693248          0 1072693248
//     *
//     0000240          0          0          0          0
//     *
//     0000440
//
//
//     % od -t fD lsjt_operator_test_relative_identity_Nmax02.bin
//
//     0000000                   5e-324                        0
//     0000020                   1e-323         6.365987374e-314
//     0000040                        1                        0
//     0000060                        1                        1
//     0000100                        0                        1
//     0000120                        0                        0
//     0000140                        1                        1
//     *
//     0000240                        0                        0
//     *
//     0000440
{
  ////////////////////////////////////////////////////////////////
  // defining operator
  ////////////////////////////////////////////////////////////////

  std::cout << "Setup" << std::endl;

  // set tensorial labels
  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0=0;
  operator_labels.g0=0;
  operator_labels.symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
  operator_labels.T0_min=0;
  operator_labels.T0_max=2;

  // set basis parameters
  int Nmax=2;
  int Jmax = Nmax+1;

  // set up relative space
  basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);

  // populate operator containers
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_component_matrices;
  basis::ConstructIdentityOperatorRelativeLSJT(
      operator_labels,
      relative_space,
      relative_component_sectors,
      relative_component_matrices
    );

  basis::WriteRelativeOperatorLSJT(
      filename,
      relative_space,
      operator_labels,relative_component_sectors,relative_component_matrices,
      true  // verbose
    );

}

void ReadTestRelativeManual(const std::string& filename)
// "Manual" test of reading...  Supplanted by use of
// ReadRelativeOperatorLSJT().
{
  ////////////////////////////////////////////////////////////////
  // read test
  ////////////////////////////////////////////////////////////////

  std::cout << "Readback test" << std::endl;

  // set up stream for readback
  std::ifstream is(filename.c_str());

  // read header parameters
  basis::RelativeOperatorParametersLSJT operator_parameters;
  basis::ReadRelativeOperatorParametersLSJT(is,operator_parameters);
  int J0=operator_parameters.J0;
  int g0=operator_parameters.g0;
  // and inspect...
  basis::WriteRelativeOperatorParametersLSJT(std::cout,operator_parameters);

  // set up space
  basis::RelativeSpaceLSJT space(operator_parameters.Nmax,operator_parameters.Jmax);

  // read matrices
  std::array<basis::RelativeSectorsLSJT,3> component_sectors;
  std::array<basis::OperatorBlocks<double>,3> component_matrices;

  for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
    // for each isospin component
    {
      // enumerate sectors
      component_sectors[T0] = basis::RelativeSectorsLSJT(space,J0,T0,g0);

      // read matrices
      basis::ReadRelativeOperatorComponentLSJT(
          is,
          T0,
          component_sectors[T0],component_matrices[T0]
        );
      // and inspect...
      basis::WriteRelativeOperatorComponentLSJT(
          std::cout,
          T0,
          component_sectors[T0],component_matrices[T0]
        );
    }
}


void IdentityTestOLD()
{
  ////////////////////////////////////////////////////////////////
  // construct relative identity operator
  ////////////////////////////////////////////////////////////////

  // define operator properties
  int Nmax_relative = 4;
  int Jmax_relative = Nmax_relative+1;

  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0 = 0;
  operator_labels.g0 = 0;
  operator_labels.T0_min = 0;
  operator_labels.T0_max = 2;
  operator_labels.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

  // define space and operator containers
  basis::RelativeSpaceLSJT relative_space(Nmax_relative,Jmax_relative);
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_component_matrices;

  // do construction
  ConstructIdentityOperatorRelativeLSJT(
      operator_labels,
      relative_space,relative_component_sectors,relative_component_matrices
    );

  // and now it would be nice to inspect the contents, wouldn't it,
  // ...

  // try out deletion
  basis::ClearOperatorJT(relative_component_sectors,relative_component_matrices);

}

void ReadTestRelative(const std::string& filename)
{
  basis::RelativeSpaceLSJT relative_space;
  basis::RelativeOperatorParametersLSJT operator_parameters;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_component_matrices;
  basis::ReadRelativeOperatorLSJT(
      filename,
      relative_space,
      operator_parameters,relative_component_sectors,relative_component_matrices,
      true  // verbose
    );
}

void WriteTestRelativeCMManual(const std::string& filename)
{
  ////////////////////////////////////////////////////////////////
  // defining operator
  ////////////////////////////////////////////////////////////////

  std::cout << "Setup" << std::endl;

  // set up space
  int Nmax = 2;
  basis::RelativeCMSpaceLSJT space(Nmax);

  // set up operator containers
  //
  // These are vectors to store information for T0=0/1/2 components.
  std::array<basis::RelativeCMSectorsLSJT,3> component_sectors;
  std::array<basis::OperatorBlocks<double>,3> component_matrices;

  // populate operator containers
  int J0 = 0;
  int g0 = 0;
  int T0_min=0;
  int T0_max=0;
  for (int T0=T0_min; T0<=T0_max; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      component_sectors[T0] = basis::RelativeCMSectorsLSJT(space,J0,T0,g0);
      std::cout << " T0 " << T0 << " size " << component_sectors[T0].size() << std::endl;

      // populate matrices
      if (T0==0)
        basis::SetOperatorToIdentity(component_sectors[T0],component_matrices[T0]);
      else
        basis::SetOperatorToZero(component_sectors[T0],component_matrices[T0]);
    }

  // set up stream for output
  std::ostringstream os;

  // write matrices
  for (int T0=T0_min; T0<=T0_max; ++T0)
    {
      basis::WriteRelativeCMOperatorComponentLSJT(
          os,
          T0,
          component_sectors[T0],component_matrices[T0]
        );
    }

  // dump to file
  std::ofstream ofile(filename.c_str());
  ofile << os.str();

}

void WriteTestRelativeCM(const std::string& filename)
// Write using WriteRelativeOperatorLSJT().
{
  ////////////////////////////////////////////////////////////////
  // defining operator
  ////////////////////////////////////////////////////////////////

  std::cout << "Setup" << std::endl;

  // set tensorial labels
  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0=0;
  operator_labels.g0=0;
  operator_labels.symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
  operator_labels.T0_min=0;
  operator_labels.T0_max=2;

  // set basis parameters
  int Nmax=2;
  int Jmax = Nmax+1;

  // set up relative space
  basis::RelativeCMSpaceLSJT relative_cm_space(Nmax);

  // populate operator containers
  std::array<basis::RelativeCMSectorsLSJT,3> relative_cm_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_cm_component_matrices;
  basis::ConstructIdentityOperatorRelativeCMLSJT(
      operator_labels,
      relative_cm_space,
      relative_cm_component_sectors,
      relative_cm_component_matrices
    );

  basis::WriteRelativeCMOperatorLSJT(
      filename,
      relative_cm_space,
      operator_labels,relative_cm_component_sectors,relative_cm_component_matrices,
      true  // verbose
    );

}

void ReadTestRelativeCM(const std::string& filename)
{
  basis::RelativeCMSpaceLSJT relative_cm_space;
  basis::RelativeCMOperatorParametersLSJT operator_parameters;
  std::array<basis::RelativeCMSectorsLSJT,3> relative_cm_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_cm_component_matrices;
  basis::ReadRelativeCMOperatorLSJT(
      filename,
      relative_cm_space,
      operator_parameters,relative_cm_component_sectors,relative_cm_component_matrices,
      true  // verbose
    );
}

void WriteTestTwoBody(const std::string& filename)
{
  ////////////////////////////////////////////////////////////////
  // defining operator
  ////////////////////////////////////////////////////////////////

  // since we set the operator to a "naive identity" operator, the
  // matrix elements are taken to be NAS

  std::cout << "Setup" << std::endl;

  // set up space
  int Nmax = 2;
  basis::TwoBodySpaceLSJT space(basis::Rank::kTwoBody,Nmax);

  // set up operator containers
  //
  // These are vectors to store information for T0=0/1/2 components.
  std::array<basis::TwoBodySectorsLSJT,3> component_sectors;
  std::array<basis::OperatorBlocks<double>,3> component_matrices;

  // populate operator containers
  int J0 = 0;
  int g0 = 0;
  int T0_min=0;
  int T0_max=0;
  for (int T0=T0_min; T0<=T0_max; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      component_sectors[T0] = basis::TwoBodySectorsLSJT(space,J0,T0,g0);
      std::cout << " T0 " << T0 << " size " << component_sectors[T0].size() << std::endl;

      // populate matrices
      if (T0==0)
        basis::SetOperatorToIdentity(component_sectors[T0],component_matrices[T0]);
      else
        basis::SetOperatorToZero(component_sectors[T0],component_matrices[T0]);
    }

  // set up stream for output
  std::ostringstream os;

  // write matrices
  for (int T0=T0_min; T0<=T0_max; ++T0)
    {
      basis::WriteTwoBodyOperatorComponentLSJT(
          os,
          T0,
          component_sectors[T0],component_matrices[T0],
          basis::NormalizationConversion::kNone
        );
    }

  // dump to file
  std::ofstream ofile(filename.c_str());
  ofile << os.str();
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  std::string relative_filename_text("lsjt_operator_test_relative_identity_Nmax02.dat");
  WriteTestRelative(relative_filename_text);
  ReadTestRelative(relative_filename_text);

  std::string relative_filename_binary("lsjt_operator_test_relative_identity_Nmax02.bin");
  WriteTestRelative(relative_filename_binary);
  ReadTestRelative(relative_filename_binary);

  std::string relative_cm_filename_text("lsjt_operator_test_relative_cm_identity_Nmax02.dat");
  WriteTestRelativeCM(relative_cm_filename_text);
  ReadTestRelativeCM(relative_cm_filename_text);

  std::string relative_cm_filename_binary("lsjt_operator_test_relative_cm_identity_Nmax02.bin");
  WriteTestRelativeCM(relative_cm_filename_binary);
  ReadTestRelativeCM(relative_cm_filename_binary);

  std::string two_body_filename("lsjt_operator_test_two_body_identity_nas_Nmax02.dat");
  WriteTestTwoBody(two_body_filename);

  // termination
  return 0;
}
