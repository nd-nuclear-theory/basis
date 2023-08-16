/****************************************************************

  lsjt_operator.cpp

  Mark A. Caprio, Patrick J. Fasano
  University of Notre Dame.

****************************************************************/

#include "lsjt_operator.h"

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <tuple>
#include <utility>
#include <vector>

#include "eigen3/Eigen/Core"

#include "basis.h"
#include "jt_operator.h"
#include "lsjt_scheme.h"
#include "many_body.h"
#include "operator.h"
#include "mcutils/parsing.h"
#include "mcutils/io.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // relative two-body operator
  ////////////////////////////////////////////////////////////////

  void ReadRelativeOperatorParametersLSJT(
      std::istream& is,
      basis::RelativeOperatorParametersLSJT& parameters,
      mcutils::IOMode io_mode
    )
  {

    std::string line;

    if (io_mode == mcutils::IOMode::kText)
    {
      // line 1: version -- but first gobble any comment lines
      while (std::getline(is,line), line[0]=='#') {};
      int version;
      std::stringstream(line) >> version;
      assert(version==1);

      // line 2: operator tensor properties
      std::getline(is,line);
      int symmetry_phase_mode_int;
      std::stringstream(line)
        >> parameters.J0 >> parameters.g0
        >> parameters.T0_min >> parameters.T0_max
        >> symmetry_phase_mode_int;
      assert (symmetry_phase_mode_int==0);
      parameters.symmetry_phase_mode = basis::SymmetryPhaseMode(symmetry_phase_mode_int);

      // line 3: relative basis truncation
      std::getline(is,line);
      std::stringstream(line) >> parameters.Nmax >> parameters.Jmax;
    }
    else // if (io_mode == mcutils::IOMode::kBinary)
    {
      // field 1: version
      mcutils::VerifyBinary<int32_t>(is, 1, "Invalid version", "version");

      // fields 2-6: operator tensor properties
      mcutils::ReadBinary<int32_t>(is, parameters.J0);
      mcutils::ReadBinary<int32_t>(is, parameters.g0);
      mcutils::ReadBinary<int32_t>(is, parameters.T0_min);
      mcutils::ReadBinary<int32_t>(is, parameters.T0_max);
      mcutils::VerifyBinary<int32_t>(is, 0, "Unsupported symmetry phase mode", "symmetry phase mode");
      parameters.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

      // fields 7-8: relative basis truncation
      mcutils::ReadBinary<int32_t>(is, parameters.Nmax);
      mcutils::ReadBinary<int32_t>(is, parameters.Jmax);
    }
  }

  void WriteRelativeOperatorParametersLSJT(
      std::ostream& os,
      const basis::RelativeOperatorParametersLSJT& parameters,
      mcutils::IOMode io_mode
    )
  {
    int version = 1;
    if (io_mode == mcutils::IOMode::kText)
    {
      os
        << "# RELATIVE LSJT" << std::endl
        << "#   version" << std::endl
        << "#   J0 g0 T0_min T0_max symmetry_phase_mode  [P0=(-)^g0]" << std::endl
        << "#   Nmax Jmax" << std::endl
        << "#   T0   N' L' S' J' T'   N L S J T   JT-RME" << std::endl
        << " " << version << std::endl
        << " " << parameters.J0 << " " << parameters.g0
        << " " << parameters.T0_min << " " << parameters.T0_max
        << " " << int(parameters.symmetry_phase_mode) << std::endl
        << " " << parameters.Nmax << " " << parameters.Jmax << std::endl;
    }
    else  // if (io_mode == mcutils::IOMode::kBinary)
    {
      // field 1: version
      mcutils::WriteBinary<int32_t>(os, version);

      // fields 2-6: operator tensor properties
      mcutils::WriteBinary<int32_t>(os, parameters.J0);
      mcutils::WriteBinary<int32_t>(os, parameters.g0);
      mcutils::WriteBinary<int32_t>(os, parameters.T0_min);
      mcutils::WriteBinary<int32_t>(os, parameters.T0_max);
      assert(parameters.symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian);
      mcutils::WriteBinary<int32_t>(os, static_cast<int32_t>(parameters.symmetry_phase_mode));

      // fields 7-8: relative basis truncation
      mcutils::WriteBinary<int32_t>(os, parameters.Nmax);
      mcutils::WriteBinary<int32_t>(os, parameters.Jmax);
    }
  }

  void ReadRelativeOperatorComponentLSJT(
      std::istream& is,
      int T0,
      const RelativeSectorsLSJT& sectors,
      OperatorBlocks<double>& matrices,
      mcutils::IOMode io_mode
    )
  {

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
        const typename RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
        const typename RelativeSectorsLSJT::SubspaceType& bra_subspace = sector.bra_subspace();
        const typename RelativeSectorsLSJT::SubspaceType& ket_subspace = sector.ket_subspace();

        // construct zero matrix for sector
        Eigen::MatrixXd sector_matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());

        std::vector<double> buffer;
        if (io_mode == mcutils::IOMode::kBinary)
          {
                  // calculate number of matrix elements in sector
            std::size_t sector_entries = 0;
            if (sector.IsDiagonal())
              // diagonal sector
              {
                std::size_t dimension = ket_subspace.size();
                sector_entries = dimension*(dimension+1)/2;
              }
            else  // if (sector.IsUpperTriangle())
              // upper triangle sector (but not diagonal)
              {
                std::size_t bra_dimension = bra_subspace.size();
                std::size_t ket_dimension = ket_subspace.size();
                sector_entries = bra_dimension*ket_dimension;
              }

            // read entire sector to temporary buffer
            buffer.resize(sector_entries, 0.);
            mcutils::ReadBinary<double>(is, buffer.data(), sector_entries);

          }

        std::size_t i = 0;
        // retrieve matrix elements
        for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
          for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
            {

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              double matrix_element;
              if (io_mode == mcutils::IOMode::kText)
              {
                // define states
                const basis::RelativeStateLSJT bra(bra_subspace,bra_index);
                const basis::RelativeStateLSJT ket(ket_subspace,ket_index);

                // read numbers from input line
                int
                  input_T0,
                  bra_N,
                  bra_L,
                  bra_S,
                  bra_J,
                  bra_T,
                  ket_N,
                  ket_L,
                  ket_S,
                  ket_J,
                  ket_T;
                is
                  >> input_T0
                  >> bra_N
                  >> bra_L
                  >> bra_S
                  >> bra_J
                  >> bra_T
                  >> ket_N
                  >> ket_L
                  >> ket_S
                  >> ket_J
                  >> ket_T
                  >> matrix_element;

                // validate labels
                bool expected_labels = true
                  && (input_T0==T0)
                  && (bra_N==bra.N())
                  && (bra_L==bra.L())
                  && (bra_S==bra.S())
                  && (bra_J==bra.J())
                  && (bra_T==bra.T())
                  && (ket_N==ket.N())
                  && (ket_L==ket.L())
                  && (ket_S==ket.S())
                  && (ket_J==ket.J())
                  && (ket_T==ket.T());
                assert(expected_labels);
              }
              else
              {
                matrix_element = buffer[i++];
              }

              // save matrix element
              sector_matrix(bra_index,ket_index) = matrix_element;
            }

        // store matrix for sector
        matrices.push_back(sector_matrix);

      }
  }

  void WriteRelativeOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const RelativeSectorsLSJT& sectors,
      const OperatorBlocks<double>& matrices,
      mcutils::IOMode io_mode
    )
  {

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
        const typename RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
        const typename RelativeSectorsLSJT::SubspaceType& bra_subspace = sector.bra_subspace();
        const typename RelativeSectorsLSJT::SubspaceType& ket_subspace = sector.ket_subspace();

        // verify that sector is canonical
        //
        // This is a check that the caller's sector construction
        // followed the specification that only "upper triangle"
        // sectors are stored.
        assert(sector.bra_subspace_index()<=sector.ket_subspace_index());

        // temporary buffer for binary write
        std::vector<double> buffer;
        std::size_t sector_entries = 0;
        if (io_mode == mcutils::IOMode::kBinary)
          {
            // calculate number of matrix elements in sector
            if (sector.IsDiagonal())
              // diagonal sector
              {
                const std::size_t& dimension = ket_subspace.size();
                sector_entries = dimension*(dimension+1)/2;
              }
            else  // if (sector.IsUpperTriangle())
              // upper triangle sector (but not diagonal)
              {
                const std::size_t& bra_dimension = bra_subspace.size();
                const std::size_t& ket_dimension = ket_subspace.size();
                sector_entries = bra_dimension*ket_dimension;
              }

            // allocate buffer
            buffer.resize(sector_entries, 0.);
          }

        // iterate over matrix elements
        std::size_t i = 0;
        for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
          for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
            {

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // define states
              const basis::RelativeStateLSJT bra(bra_subspace,bra_index);
              const basis::RelativeStateLSJT ket(ket_subspace,ket_index);

              // extract matrix element
              const double matrix_element = matrices[sector_index](bra_index,ket_index);

              if (io_mode == mcutils::IOMode::kText)
                {
                  // generate output line
                  const int width = 3;
                  const int precision = 14;  // less than 16 to provide some roundoff and avoid ugliness on doubles
                  os << std::setprecision(precision);
                  os
                    << " " << std::setw(width) << T0
                    << " " << "  "
                    << " " << std::setw(width) << bra.N()
                    << " " << std::setw(width) << bra.L()
                    << " " << std::setw(width) << bra.S()
                    << " " << std::setw(width) << bra.J()
                    << " " << std::setw(width) << bra.T()
                    << " " << "  "
                    << " " << std::setw(width) << ket.N()
                    << " " << std::setw(width) << ket.L()
                    << " " << std::setw(width) << ket.S()
                    << " " << std::setw(width) << ket.J()
                    << " " << std::setw(width) << ket.T()
                    << " " << "  "
                    << " " << std::showpoint << std::scientific << matrix_element
                    << std::endl;
                }
              else  // if (io_mode == mcutils::IOMode::kText)
                {
                  buffer[i++] = matrix_element;
                }
            }

        // write temporary buffer to file
        if (io_mode == mcutils::IOMode::kBinary)
          mcutils::WriteBinary<double>(os,buffer.data(),sector_entries);
      }
  }

  void ReadRelativeOperatorLSJT(
      const std::string& relative_filename,
      basis::RelativeSpaceLSJT& relative_space,
      basis::RelativeOperatorParametersLSJT& operator_parameters,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      bool verbose
    )
  // FUTURE: change to line-based input with std::getline and parsing checks for easier debugging
  // FUTURE: check file status on open
  {

    // deduce I/O mode
    mcutils::IOMode io_mode = mcutils::DeducedIOMode(relative_filename);

    // open stream for reading
    if (verbose)
      {
        std::cout
          << "Reading relative operator file..." << std::endl
          << "  Filename: " << relative_filename << std::endl
          << "  Mode: " << ((io_mode == mcutils::IOMode::kText) ? "text" : "binary")
          << std::endl;
      }
    std::ios_base::openmode mode_argument = std::ios_base::in;
    if (io_mode==mcutils::IOMode::kBinary)
      mode_argument |= std::ios_base::binary;
    std::ifstream is(relative_filename.c_str(), mode_argument);

    // read header parameters
    basis::ReadRelativeOperatorParametersLSJT(is,operator_parameters,io_mode);
    if (verbose)
      {
        std::cout
          << "  Operator properties:"
          << " J0 " << operator_parameters.J0
          << " g0 " << operator_parameters.g0
          << " T0_min " << operator_parameters.T0_min
          << " T0_max " << operator_parameters.T0_max
          << " symmetry " << int(operator_parameters.symmetry_phase_mode)
          << std::endl
          << "  Truncation:"
          << " Nmax " << operator_parameters.Nmax
          << " Jmax " << operator_parameters.Jmax
          << std::endl;
      }

    // set up relative space
    relative_space = basis::RelativeSpaceLSJT(
        operator_parameters.Nmax,operator_parameters.Jmax
      );

    // populate sectors and matrices
    for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
      // for each isospin component
      {
        // enumerate sectors
        relative_component_sectors[T0]
          = basis::RelativeSectorsLSJT(relative_space,operator_parameters.J0,T0,operator_parameters.g0);

        // read matrices
        basis::ReadRelativeOperatorComponentLSJT(
            is,
            T0,
            relative_component_sectors[T0],relative_component_matrices[T0],
            io_mode
          );
      }

    // write diagnostics
    if (verbose)
      {
        std::cout << "  Matrix elements:";
        for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
          std::cout << " " << basis::UpperTriangularEntries(relative_component_sectors[T0]);
        std::cout << std::endl;
        std::cout << "  Allocated:";
        for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
          std::cout << " " << basis::AllocatedEntries(relative_component_matrices[T0]);
        std::cout << std::endl;
      }
  }

  void WriteRelativeOperatorLSJT(
      const std::string& relative_filename,
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::OperatorLabelsJT& operator_labels,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      bool verbose
    )
  // FUTURE: check file status on open
  {

    // deduce I/O mode
    mcutils::IOMode io_mode = mcutils::DeducedIOMode(relative_filename);

    // open stream for writing
    if (verbose)
      {
        std::cout
          << "Writing relative operator file..." << std::endl
          << "  Filename: " << relative_filename << std::endl;
      }
    std::ios_base::openmode mode_argument = std::ios_base::out;
    if (io_mode==mcutils::IOMode::kBinary)
      mode_argument |= std::ios_base::binary;
    std::ofstream os(relative_filename.c_str(), mode_argument);

    // set up full operator file parameters
    int Nmax = relative_space.Nmax();
    int Jmax = relative_space.Jmax();
    RelativeOperatorParametersLSJT operator_parameters(operator_labels,Nmax,Jmax);
    if (verbose)
      {
        std::cout
          << "  Operator properties:"
          << " J0 " << operator_parameters.J0
          << " g0 " << operator_parameters.g0
          << " T0_min " << operator_parameters.T0_min
          << " T0_max " << operator_parameters.T0_max
          << " symmetry " << int(operator_parameters.symmetry_phase_mode)
          << std::endl
          << "  Truncation:"
          << " Nmax " << operator_parameters.Nmax
          << " Jmax " << operator_parameters.Jmax
          << std::endl;
      }

    // write header parameters
    basis::WriteRelativeOperatorParametersLSJT(os,operator_parameters,io_mode);

    // write matrices
    for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
      {
        basis::WriteRelativeOperatorComponentLSJT(
            os,
            T0,
            relative_component_sectors[T0],relative_component_matrices[T0],
            io_mode
          );
      }

  }

  ////////////////////////////////////////////////////////////////
  // relative LSJT operator construction
  ////////////////////////////////////////////////////////////////

  void ConstructZeroOperatorRelativeLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices
    )
  {

    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate sectors
        relative_component_sectors[T0]
          = basis::RelativeSectorsLSJT(relative_space,operator_labels.J0,T0,operator_labels.g0);

        // populate matrices
        relative_component_matrices[T0].resize(relative_component_sectors[T0].size());
        basis::SetOperatorToZero(relative_component_sectors[T0],relative_component_matrices[T0]);
      }
  }

  void ConstructIdentityOperatorRelativeLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices
    )
  {

    // validate operator labels
    assert(operator_labels.J0==0);
    assert(operator_labels.g0==0);
    assert(operator_labels.symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian);

    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate sectors
        relative_component_sectors[T0]
          = basis::RelativeSectorsLSJT(relative_space,operator_labels.J0,T0,operator_labels.g0);

        // populate matrices
        relative_component_matrices[T0].resize(relative_component_sectors[T0].size());
        if (T0==0)
          // identity matrices in
          basis::SetOperatorToIdentity(relative_component_sectors[T0],relative_component_matrices[T0]);
        else
          basis::SetOperatorToZero(relative_component_sectors[T0],relative_component_matrices[T0]);
      }
  }

  ////////////////////////////////////////////////////////////////
  // relative-cm LSJT operator manipulation
  ////////////////////////////////////////////////////////////////

  inline
  void RecastLabelsRelativeCMLSJTToRelativeCMLSJTN(
      const RelativeCMSubspaceLSJTLabels& relative_cm_lsjt_subspace_labels,
      const RelativeCMStateLSJTLabels& relative_cm_lsjt_state_labels,
      RelativeCMSubspaceLSJTNLabels& relative_cm_lsjtn_subspace_labels,
      RelativeCMStateLSJTNLabels& relative_cm_lsjtn_state_labels
    )
  // Recast labels for (subspace,state) from RelativeCMLSJT scheme to
  // RelativeCMLSJTN scheme.
  //
  // This direction of conversion (from target labels to source
  // labels) is as needed for looking up the source matrix element.
  //
  // This switches N from being the most-significant state label (that
  // is, least-rapidly-varying in the lexicographic ordering scheme),
  // to being the least-significant subspace label (that is,
  // most-rapidly-varying):
  //
  //   (L,S,J,T,g) : ([N],N1,l1,N2,l2) -> (L,S,J,T,g,N) : (N1,l1,N2,l2)
  //
  // As a state label, N is actually an implicit label, not stored,
  // but effectively the first label for ordering purposes.  Its value
  // may be recovered as N1+N2.
  //
  // Arguments:
  //   relative_cm_lsjt_subspace_labels (...) : source subspace labels
  //   relative_cm_lsjt_state_labels (...) : source state labels
  //   relative_cm_lsjtn_subspace_labels (...,output) : target subspace labels
  //   relative_cm_lsjtn_state_labels (...,output) : target state labels
  {
    // extract labels
    int L, S, J, T, g;
    std::tie(L,S,J,T,g) = relative_cm_lsjt_subspace_labels;
    int Nr, lr, Nc, lc;
    std::tie(Nr,lr,Nc,lc) = relative_cm_lsjt_state_labels;

    // repackage labels
    int N = Nr+Nc;
    relative_cm_lsjtn_subspace_labels = RelativeCMSubspaceLSJTNLabels(L,S,J,T,g,N);
    relative_cm_lsjtn_state_labels = relative_cm_lsjt_state_labels;
  }


  inline
  void RecastLabelsRelativeCMLSJTNToRelativeCMLSJT(
      const RelativeCMSubspaceLSJTNLabels& relative_cm_lsjtn_subspace_labels,
      const RelativeCMStateLSJTNLabels& relative_cm_lsjtn_state_labels,
      RelativeCMSubspaceLSJTLabels& relative_cm_lsjt_subspace_labels,
      RelativeCMStateLSJTLabels& relative_cm_lsjt_state_labels
    )
  // Recast labels for (subspace,state) from RelativeCMLSJTN scheme to
  // RelativeCMLSJT scheme.
  //
  // This direction of conversion (from target labels to source
  // labels) is as needed for looking up the source matrix element.
  //
  // This switches N from being the least-significant state label (that
  // is, most-rapidly-varying in the lexicographic ordering scheme),
  // to being the most-significant subspace label (that is,
  // least-rapidly-varying):
  //
  //   (L,S,J,T,g,N) : (N1,l1,N2,l2) -> (L,S,J,T,g) : ([N],N1,l1,N2,l2)
  //
  // As a state label, N is actually an implicit label, not stored,
  // but effectively the last label for ordering purposes.  Its value
  // may be recovered as Nr+Nc.
  //
  // Arguments:
  //   relative_cm_lsjtn_subspace_labels (...,output) : source subspace labels
  //   relative_cm_lsjtn_state_labels (...,output) : source state labels
  //   relative_cm_lsjt_subspace_labels (...) : target subspace labels
  //   relative_cm_lsjt_state_labels (...) : target state labels
  {
    // extract labels
    int L, S, J, T, g, N;
    std::tie(L,S,J,T,g,N) = relative_cm_lsjtn_subspace_labels;
    int Nr, lr, Nc, lc;
    std::tie(Nr,lr,Nc,lc) = relative_cm_lsjtn_state_labels;

    // repackage labels
    relative_cm_lsjt_subspace_labels = RelativeCMSubspaceLSJTLabels(L,S,J,T,g);
    relative_cm_lsjt_state_labels = relative_cm_lsjtn_state_labels;
  }


  void GatherOperatorRelativeCMLSJTNToRelativeCMLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
      const std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& relative_cm_lsjtn_component_matrices,
      const basis::RelativeCMSpaceLSJT& relative_cm_lsjt_space,
      std::array<basis::RelativeCMSectorsLSJT,3>& relative_cm_lsjt_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_cm_lsjt_component_matrices
    )
  {
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate target sectors
        relative_cm_lsjt_component_sectors[T0]
          = basis::RelativeCMSectorsLSJT(relative_cm_lsjt_space,operator_labels.J0,T0,operator_labels.g0);

        // populate matrices
        relative_cm_lsjt_component_matrices[T0].resize(relative_cm_lsjt_component_sectors[T0].size());
        for (std::size_t sector_index=0; sector_index<relative_cm_lsjt_component_sectors[T0].size(); ++sector_index)
          {
            // retrieve target sector
            const basis::RelativeCMSectorsLSJT::SectorType& relative_cm_lsjt_sector
              = relative_cm_lsjt_component_sectors[T0].GetSector(sector_index);

            // initialize matrix
            Eigen::MatrixXd& relative_cm_lsjt_matrix = relative_cm_lsjt_component_matrices[T0][sector_index];
            relative_cm_lsjt_matrix = Eigen::MatrixXd::Zero(
                relative_cm_lsjt_sector.bra_subspace().size(),
                relative_cm_lsjt_sector.ket_subspace().size()
              );

            // populate matrix elements
            for (std::size_t bra_index = 0; bra_index < relative_cm_lsjt_sector.bra_subspace().size(); ++bra_index)
              for (std::size_t ket_index = 0; ket_index < relative_cm_lsjt_sector.ket_subspace().size(); ++ket_index)
                // for each target matrix element
                {

                  // ensure canonical matrix element if diagonal sector
                  if (relative_cm_lsjt_sector.IsDiagonal())
                    if (!(bra_index<=ket_index))
                      continue;

                  // retrieve target states
                  basis::RelativeCMStateLSJT relative_cm_lsjt_bra(relative_cm_lsjt_sector.bra_subspace(),bra_index);
                  basis::RelativeCMStateLSJT relative_cm_lsjt_ket(relative_cm_lsjt_sector.ket_subspace(),ket_index);

                  // extract source bra labels
                  RelativeCMSubspaceLSJTLabels relative_cm_lsjt_subspace_labels_bra
                    = relative_cm_lsjt_bra.subspace().labels();
                  RelativeCMStateLSJTLabels relative_cm_lsjt_state_labels_bra
                    = relative_cm_lsjt_bra.labels();
                  RelativeCMSubspaceLSJTNLabels relative_cm_lsjtn_subspace_labels_bra;
                  RelativeCMStateLSJTNLabels relative_cm_lsjtn_state_labels_bra;
                  RecastLabelsRelativeCMLSJTToRelativeCMLSJTN(
                      relative_cm_lsjt_subspace_labels_bra,
                      relative_cm_lsjt_state_labels_bra,
                      relative_cm_lsjtn_subspace_labels_bra,
                      relative_cm_lsjtn_state_labels_bra
                    );

                  // extract source bra indices
                  std::size_t relative_cm_lsjtn_subspace_index_bra
                    = relative_cm_lsjtn_space.LookUpSubspaceIndex(
                        relative_cm_lsjtn_subspace_labels_bra
                      );
                  std::size_t relative_cm_lsjtn_state_index_bra
                    = relative_cm_lsjtn_space.GetSubspace(relative_cm_lsjtn_subspace_index_bra).LookUpStateIndex(
                        relative_cm_lsjtn_state_labels_bra
                      );

                  // extract source ket labels
                  RelativeCMSubspaceLSJTLabels relative_cm_lsjt_subspace_labels_ket
                    = relative_cm_lsjt_ket.subspace().labels();
                  RelativeCMStateLSJTLabels relative_cm_lsjt_state_labels_ket
                    = relative_cm_lsjt_ket.labels();
                  RelativeCMSubspaceLSJTNLabels relative_cm_lsjtn_subspace_labels_ket;
                  RelativeCMStateLSJTNLabels relative_cm_lsjtn_state_labels_ket;
                  RecastLabelsRelativeCMLSJTToRelativeCMLSJTN(
                      relative_cm_lsjt_subspace_labels_ket,
                      relative_cm_lsjt_state_labels_ket,
                      relative_cm_lsjtn_subspace_labels_ket,
                      relative_cm_lsjtn_state_labels_ket
                    );

                  // extract source ket indices
                  std::size_t relative_cm_lsjtn_subspace_index_ket
                    = relative_cm_lsjtn_space.LookUpSubspaceIndex(
                        relative_cm_lsjtn_subspace_labels_ket
                      );
                  std::size_t relative_cm_lsjtn_state_index_ket
                    = relative_cm_lsjtn_space.GetSubspace(relative_cm_lsjtn_subspace_index_ket).LookUpStateIndex(
                        relative_cm_lsjtn_state_labels_ket
                      );

                  // // debugging
                  // const RelativeCMSubspaceLSJTN& relative_cm_lsjtn_subspace_bra
                  //   = relative_cm_lsjtn_space.GetSubspace(relative_cm_lsjtn_subspace_index_bra);
                  // const RelativeCMSubspaceLSJTN& relative_cm_lsjtn_subspace_ket
                  //   = relative_cm_lsjtn_space.GetSubspace(relative_cm_lsjtn_subspace_index_ket);
                  // std::cout << " pre-lookup "
                  //           << " " << relative_cm_lsjtn_subspace_index_bra
                  //           << " " << relative_cm_lsjtn_subspace_index_ket
                  //           << " " << ";"
                  //           << " " << relative_cm_lsjtn_state_index_bra
                  //           << " " << relative_cm_lsjtn_state_index_ket
                  //           << " " << ";"
                  //           << " " << relative_cm_lsjtn_subspace_bra.LabelStr()
                  //           << " " << relative_cm_lsjtn_subspace_ket.LabelStr()
                  //           << " " << relative_cm_lsjtn_subspace_bra.size()
                  //           << " " << relative_cm_lsjtn_subspace_ket.size()
                  //           << std::endl;

                  // Note on canonicalization of indices for lookup
                  // (or lack thereof)
                  //
                  // Matrix elements which are canonical by LSJT
                  // sector, and ordered by N within a LSJT sector,
                  // should also be canonical by (LSJT;N) sector.
                  // That is, the upper triangle of a matrix with
                  // basis states ordered by
                  //
                  //   (L,S,J,T,g) : ([N],N1,l1,N2,l2)
                  //
                  // or
                  //
                  //   (L,S,J,T,g,N) : (N1,l1,N2,l2)
                  //
                  // should be identical.
                  //
                  // So canonicalization would only be necessary if we
                  // were to fill in a *lower* triangle matrix element
                  // of a diagonal target sector.  Then we would have
                  // to ensure that we look up a canonical (upper
                  // triangular) LSJTN sector.

                  // look up matrix element
                  std::size_t relative_cm_lsjtn_sector_index
                    = relative_cm_lsjtn_component_sectors[T0].LookUpSectorIndex(
                        relative_cm_lsjtn_subspace_index_bra,
                        relative_cm_lsjtn_subspace_index_ket
                      );

                  const Eigen::MatrixXd& relative_cm_lsjtn_matrix
                    = relative_cm_lsjtn_component_matrices[T0][relative_cm_lsjtn_sector_index];
                  double relative_cm_lsjtn_matrix_element = relative_cm_lsjtn_matrix(
                      relative_cm_lsjtn_state_index_bra,relative_cm_lsjtn_state_index_ket
                    );

                  relative_cm_lsjt_matrix(bra_index,ket_index) = relative_cm_lsjtn_matrix_element;

                }
          }
      }

  }

  void ScatterOperatorRelativeCMLSJTToRelativeCMLSJTN(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceLSJT& relative_cm_lsjt_space,
      const std::array<basis::RelativeCMSectorsLSJT,3>& relative_cm_lsjt_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& relative_cm_lsjt_component_matrices,
      const basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
      std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_cm_lsjtn_component_matrices
    )
  {
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate target sectors
        relative_cm_lsjtn_component_sectors[T0]
          = basis::RelativeCMSectorsLSJTN(relative_cm_lsjtn_space,operator_labels.J0,T0,operator_labels.g0);

        // populate matrices
        relative_cm_lsjtn_component_matrices[T0].resize(relative_cm_lsjtn_component_sectors[T0].size());
        for (std::size_t sector_index=0; sector_index<relative_cm_lsjtn_component_sectors[T0].size(); ++sector_index)
          {
            // retrieve target sector
            const basis::RelativeCMSectorsLSJTN::SectorType& relative_cm_lsjtn_sector
              = relative_cm_lsjtn_component_sectors[T0].GetSector(sector_index);

            // initialize matrix
            Eigen::MatrixXd& relative_cm_lsjtn_matrix = relative_cm_lsjtn_component_matrices[T0][sector_index];
            relative_cm_lsjtn_matrix = Eigen::MatrixXd::Zero(
                relative_cm_lsjtn_sector.bra_subspace().size(),
                relative_cm_lsjtn_sector.ket_subspace().size()
              );

            // populate matrix elements
            for (std::size_t bra_index = 0; bra_index < relative_cm_lsjtn_sector.bra_subspace().size(); ++bra_index)
              {
                // retreive target bra state
                basis::RelativeCMStateLSJTN relative_cm_lsjtn_bra(relative_cm_lsjtn_sector.bra_subspace(),bra_index);
                // extract source bra labels
                RelativeCMSubspaceLSJTNLabels relative_cm_lsjtn_subspace_labels_bra
                  = relative_cm_lsjtn_bra.subspace().labels();
                RelativeCMStateLSJTNLabels relative_cm_lsjtn_state_labels_bra
                  = relative_cm_lsjtn_bra.labels();
                RelativeCMSubspaceLSJTLabels relative_cm_lsjt_subspace_labels_bra;
                RelativeCMStateLSJTLabels relative_cm_lsjt_state_labels_bra;
                RecastLabelsRelativeCMLSJTNToRelativeCMLSJT(
                    relative_cm_lsjtn_subspace_labels_bra,
                    relative_cm_lsjtn_state_labels_bra,
                    relative_cm_lsjt_subspace_labels_bra,
                    relative_cm_lsjt_state_labels_bra
                  );

                // extract source bra indices
                std::size_t relative_cm_lsjt_subspace_index_bra
                  = relative_cm_lsjt_space.LookUpSubspaceIndex(
                      relative_cm_lsjt_subspace_labels_bra
                    );
                std::size_t relative_cm_lsjt_state_index_bra
                  = relative_cm_lsjt_space.GetSubspace(relative_cm_lsjt_subspace_index_bra).LookUpStateIndex(
                      relative_cm_lsjt_state_labels_bra
                    );

                for (std::size_t ket_index = 0; ket_index < relative_cm_lsjtn_sector.ket_subspace().size(); ++ket_index)
                  // for each target matrix element
                  {

                    // ensure canonical matrix element if diagonal sector
                    if (relative_cm_lsjtn_sector.IsDiagonal())
                      if (!(bra_index<=ket_index))
                        continue;

                    // retrieve target ket state
                    basis::RelativeCMStateLSJTN relative_cm_lsjtn_ket(relative_cm_lsjtn_sector.ket_subspace(),ket_index);

                    // extract source ket labels
                    RelativeCMSubspaceLSJTNLabels relative_cm_lsjtn_subspace_labels_ket
                      = relative_cm_lsjtn_ket.subspace().labels();
                    RelativeCMStateLSJTNLabels relative_cm_lsjtn_state_labels_ket
                      = relative_cm_lsjtn_ket.labels();
                    RelativeCMSubspaceLSJTLabels relative_cm_lsjt_subspace_labels_ket;
                    RelativeCMStateLSJTLabels relative_cm_lsjt_state_labels_ket;
                    RecastLabelsRelativeCMLSJTNToRelativeCMLSJT(
                        relative_cm_lsjtn_subspace_labels_ket,
                        relative_cm_lsjtn_state_labels_ket,
                        relative_cm_lsjt_subspace_labels_ket,
                        relative_cm_lsjt_state_labels_ket
                      );

                    // extract source ket indices
                    std::size_t relative_cm_lsjt_subspace_index_ket
                      = relative_cm_lsjt_space.LookUpSubspaceIndex(
                          relative_cm_lsjt_subspace_labels_ket
                        );
                    std::size_t relative_cm_lsjt_state_index_ket
                      = relative_cm_lsjt_space.GetSubspace(relative_cm_lsjt_subspace_index_ket).LookUpStateIndex(
                          relative_cm_lsjt_state_labels_ket
                        );

                    // // debugging
                    // const RelativeCMSubspaceLSJT& relative_cm_lsjt_subspace_bra
                    //   = relative_cm_lsjt_space.GetSubspace(relative_cm_lsjt_subspace_index_bra);
                    // const RelativeCMSubspaceLSJT& relative_cm_lsjt_subspace_ket
                    //   = relative_cm_lsjt_space.GetSubspace(relative_cm_lsjt_subspace_index_ket);
                    // std::cout << " pre-lookup "
                    //           << " " << relative_cm_lsjt_subspace_index_bra
                    //           << " " << relative_cm_lsjt_subspace_index_ket
                    //           << " " << ";"
                    //           << " " << relative_cm_lsjt_state_index_bra
                    //           << " " << relative_cm_lsjt_state_index_ket
                    //           << " " << ";"
                    //           << " " << relative_cm_lsjt_subspace_bra.LabelStr()
                    //           << " " << relative_cm_lsjt_subspace_ket.LabelStr()
                    //           << " " << relative_cm_lsjt_subspace_bra.size()
                    //           << " " << relative_cm_lsjt_subspace_ket.size()
                    //           << std::endl;

                    // Note on canonicalization of indices for lookup
                    // (or lack thereof)
                    //
                    // Matrix elements which are canonical by LSJT
                    // sector, and ordered by N within a LSJT sector,
                    // should also be canonical by (LSJT;N) sector.
                    // That is, the upper triangle of a matrix with
                    // basis states ordered by
                    //
                    //   (L,S,J,T,g) : ([N],N1,l1,N2,l2)
                    //
                    // or
                    //
                    //   (L,S,J,T,g,N) : (N1,l1,N2,l2)
                    //
                    // should be identical.
                    //
                    // So canonicalization would only be necessary if we
                    // were to fill in a *lower* triangle matrix element
                    // of a diagonal target sector.  Then we would have
                    // to ensure that we look up a canonical (upper
                    // triangular) LSJT sector.

                    // look up matrix element
                    std::size_t relative_cm_lsjt_sector_index
                      = relative_cm_lsjt_component_sectors[T0].LookUpSectorIndex(
                          relative_cm_lsjt_subspace_index_bra,
                          relative_cm_lsjt_subspace_index_ket
                        );

                    const Eigen::MatrixXd& relative_cm_lsjt_matrix
                      = relative_cm_lsjt_component_matrices[T0][relative_cm_lsjt_sector_index];
                    double relative_cm_lsjt_matrix_element = relative_cm_lsjt_matrix(
                        relative_cm_lsjt_state_index_bra,relative_cm_lsjt_state_index_ket
                      );

                    relative_cm_lsjtn_matrix(bra_index,ket_index) = relative_cm_lsjt_matrix_element;

                  }
              }
          }
      }

  }

  ////////////////////////////////////////////////////////////////
  // two-body LSJT operator manipulation
  ////////////////////////////////////////////////////////////////

  inline
  void RecastLabelsTwoBodyLSJTToTwoBodyLSJTN(
      const TwoBodySubspaceLSJTLabels& two_body_lsjt_subspace_labels,
      const TwoBodyStateLSJTLabels& two_body_lsjt_state_labels,
      TwoBodySubspaceLSJTNLabels& two_body_lsjtn_subspace_labels,
      TwoBodyStateLSJTNLabels& two_body_lsjtn_state_labels
    )
  // Recast labels for (subspace,state) from TwoBodyLSJT scheme to
  // TwoBodyLSJTN scheme.
  //
  // This direction of conversion (from target labels to source
  // labels) is as needed for looking up the source matrix element.
  //
  // This switches N from being the most-significant state label (that
  // is, least-rapidly-varying in the lexicographic ordering scheme),
  // to being the least-significant subspace label (that is,
  // most-rapidly-varying):
  //
  //   (L,S,J,T,g) : ([N],N1,l1,N2,l2) -> (L,S,J,T,g,N) : (N1,l1,N2,l2)
  //
  // As a state label, N is actually an implicit label, not stored,
  // but effectively the first label for ordering purposes.  Its value
  // may be recovered as N1+N2.
  //
  // Arguments:
  //   two_body_lsjt_subspace_labels (...) : source subspace labels
  //   two_body_lsjt_state_labels (...) : source state labels
  //   two_body_lsjtn_subspace_labels (...,output) : target subspace labels
  //   two_body_lsjtn_state_labels (...,output) : target state labels
  {
    // extract labels
    int L, S, J, T, g;
    std::tie(L,S,J,T,g) = two_body_lsjt_subspace_labels;
    int N1, l1, N2, l2;
    std::tie(N1,l1,N2,l2) = two_body_lsjt_state_labels;

    // repackage labels
    int N = N1+N2;
    two_body_lsjtn_subspace_labels = TwoBodySubspaceLSJTNLabels(L,S,J,T,g,N);
    two_body_lsjtn_state_labels = two_body_lsjt_state_labels;
  }


  void GatherOperatorTwoBodyLSJTNToTwoBodyLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceLSJTN& two_body_lsjtn_space,
      const std::array<basis::TwoBodySectorsLSJTN,3>& two_body_lsjtn_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& two_body_lsjtn_component_matrices,
      const basis::TwoBodySpaceLSJT& two_body_lsjt_space,
      std::array<basis::TwoBodySectorsLSJT,3>& two_body_lsjt_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& two_body_lsjt_component_matrices
    )
  {
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate target sectors
        two_body_lsjt_component_sectors[T0]
          = basis::TwoBodySectorsLSJT(two_body_lsjt_space,operator_labels.J0,T0,operator_labels.g0);

        // populate matrices
        two_body_lsjt_component_matrices[T0].resize(two_body_lsjt_component_sectors[T0].size());
        for (std::size_t sector_index=0; sector_index<two_body_lsjt_component_sectors[T0].size(); ++sector_index)
          {
            // retrieve target sector
            const basis::TwoBodySectorsLSJT::SectorType& two_body_lsjt_sector
              = two_body_lsjt_component_sectors[T0].GetSector(sector_index);

            // initialize matrix
            Eigen::MatrixXd& two_body_lsjt_matrix = two_body_lsjt_component_matrices[T0][sector_index];
            two_body_lsjt_matrix = Eigen::MatrixXd::Zero(
                two_body_lsjt_sector.bra_subspace().size(),
                two_body_lsjt_sector.ket_subspace().size()
              );

            // populate matrix elements
            for (std::size_t bra_index = 0; bra_index < two_body_lsjt_sector.bra_subspace().size(); ++bra_index)
              for (std::size_t ket_index = 0; ket_index < two_body_lsjt_sector.ket_subspace().size(); ++ket_index)
                // for each target matrix element
                {

                  // ensure canonical matrix element if diagonal sector
                  if (two_body_lsjt_sector.IsDiagonal())
                    if (!(bra_index<=ket_index))
                      continue;

                  // retrieve target states
                  basis::TwoBodyStateLSJT two_body_lsjt_bra(two_body_lsjt_sector.bra_subspace(),bra_index);
                  basis::TwoBodyStateLSJT two_body_lsjt_ket(two_body_lsjt_sector.ket_subspace(),ket_index);

                  // extract source bra labels
                  TwoBodySubspaceLSJTLabels two_body_lsjt_subspace_labels_bra
                    = two_body_lsjt_bra.subspace().labels();
                  TwoBodyStateLSJTLabels two_body_lsjt_state_labels_bra
                    = two_body_lsjt_bra.labels();
                  TwoBodySubspaceLSJTNLabels two_body_lsjtn_subspace_labels_bra;
                  TwoBodyStateLSJTNLabels two_body_lsjtn_state_labels_bra;
                  RecastLabelsTwoBodyLSJTToTwoBodyLSJTN(
                      two_body_lsjt_subspace_labels_bra,
                      two_body_lsjt_state_labels_bra,
                      two_body_lsjtn_subspace_labels_bra,
                      two_body_lsjtn_state_labels_bra
                    );

                  // extract source bra indices
                  std::size_t two_body_lsjtn_subspace_index_bra
                    = two_body_lsjtn_space.LookUpSubspaceIndex(
                        two_body_lsjtn_subspace_labels_bra
                      );
                  std::size_t two_body_lsjtn_state_index_bra
                    = two_body_lsjtn_space.GetSubspace(two_body_lsjtn_subspace_index_bra).LookUpStateIndex(
                        two_body_lsjtn_state_labels_bra
                      );

                  // extract source ket labels
                  TwoBodySubspaceLSJTLabels two_body_lsjt_subspace_labels_ket
                    = two_body_lsjt_ket.subspace().labels();
                  TwoBodyStateLSJTLabels two_body_lsjt_state_labels_ket
                    = two_body_lsjt_ket.labels();
                  TwoBodySubspaceLSJTNLabels two_body_lsjtn_subspace_labels_ket;
                  TwoBodyStateLSJTNLabels two_body_lsjtn_state_labels_ket;
                  RecastLabelsTwoBodyLSJTToTwoBodyLSJTN(
                      two_body_lsjt_subspace_labels_ket,
                      two_body_lsjt_state_labels_ket,
                      two_body_lsjtn_subspace_labels_ket,
                      two_body_lsjtn_state_labels_ket
                    );

                  // extract source ket indices
                  std::size_t two_body_lsjtn_subspace_index_ket
                    = two_body_lsjtn_space.LookUpSubspaceIndex(
                        two_body_lsjtn_subspace_labels_ket
                      );
                  std::size_t two_body_lsjtn_state_index_ket
                    = two_body_lsjtn_space.GetSubspace(two_body_lsjtn_subspace_index_ket).LookUpStateIndex(
                        two_body_lsjtn_state_labels_ket
                      );

                  // Note on canonicalization of indices for lookup
                  // (or lack thereof)
                  //
                  // Matrix elements which are canonical by LSJT
                  // sector, and ordered by N within a LSJT sector,
                  // should also be canonical by (LSJT;N) sector.
                  // That is, the upper triangle of a matrix with
                  // basis states ordered by
                  //
                  //   (L,S,J,T,g) : ([N],N1,l1,N2,l2)
                  //
                  // or
                  //
                  //   (L,S,J,T,g,N) : (N1,l1,N2,l2)
                  //
                  // should be identical.
                  //
                  // So canonicalization would only be necessary if we
                  // were to fill in a *lower* triangle matrix element
                  // of a diagonal target sector.  Then we would have
                  // to ensure that we look up a canonical (upper
                  // triangular) LSJTN sector.

                  // look up matrix element
                  std::size_t two_body_lsjtn_sector_index
                    = two_body_lsjtn_component_sectors[T0].LookUpSectorIndex(
                        two_body_lsjtn_subspace_index_bra,
                        two_body_lsjtn_subspace_index_ket
                      );

                  const Eigen::MatrixXd& two_body_lsjtn_matrix
                    = two_body_lsjtn_component_matrices[T0][two_body_lsjtn_sector_index];
                  double two_body_lsjtn_matrix_element = two_body_lsjtn_matrix(
                      two_body_lsjtn_state_index_bra,two_body_lsjtn_state_index_ket
                    );

                  two_body_lsjt_matrix(bra_index,ket_index) = two_body_lsjtn_matrix_element;

                }
          }
      }

  }

  ////////////////////////////////////////////////////////////////
  // relative-cm LSJT operator I/O
  ////////////////////////////////////////////////////////////////

  void ReadRelativeCMOperatorParametersLSJT(
      std::istream& is,
      basis::RelativeCMOperatorParametersLSJT& parameters,
      mcutils::IOMode io_mode
    )
  {

    std::string line;

    if (io_mode == mcutils::IOMode::kText)
    {
      // line 1: version -- but first gobble any comment lines
      while (std::getline(is,line), line[0]=='#') {};
      int version;
      std::stringstream(line) >> version;
      assert(version==1);

      // line 2: operator tensor properties
      std::getline(is,line);
      int symmetry_phase_mode_int;
      std::stringstream(line)
        >> parameters.J0 >> parameters.g0
        >> parameters.T0_min >> parameters.T0_max
        >> symmetry_phase_mode_int;
      assert (symmetry_phase_mode_int==0);
      parameters.symmetry_phase_mode = basis::SymmetryPhaseMode(symmetry_phase_mode_int);

      // line 3: relative basis truncation
      std::getline(is,line);
      std::stringstream(line) >> parameters.Nmax;
    }
    else // if (io_mode == mcutils::IOMode::kBinary)
    {
      // field 1: version
      mcutils::VerifyBinary<int32_t>(is, 1, "Invalid version", "version");

      // fields 2-6: operator tensor properties
      int symmetry_phase_mode_int;
      mcutils::ReadBinary<int32_t>(is, parameters.J0);
      mcutils::ReadBinary<int32_t>(is, parameters.g0);
      mcutils::ReadBinary<int32_t>(is, parameters.T0_min);
      mcutils::ReadBinary<int32_t>(is, parameters.T0_max);
      mcutils::ReadBinary<int32_t>(is, symmetry_phase_mode_int);
      assert(symmetry_phase_mode_int==0);
      parameters.symmetry_phase_mode = basis::SymmetryPhaseMode(symmetry_phase_mode_int);

      // field 7: relative basis truncation
      mcutils::ReadBinary<int32_t>(is, parameters.Nmax);
    }
  }

  void WriteRelativeCMOperatorParametersLSJT(
      std::ostream& os,
      const basis::RelativeCMOperatorParametersLSJT& parameters,
      mcutils::IOMode io_mode
    )
  {
    int version = 1;
    if (io_mode == mcutils::IOMode::kText)
    {
      os
        << "# RELATIVE-CM LSJT" << std::endl
        << "#   version" << std::endl
        << "#   J0 g0 T0_min T0_max symmetry_phase_mode  [P0=(-)^g0]" << std::endl
        << "#   Nmax" << std::endl
        << "#   T0   Nr' lr' Nc' lc' L' S' J' T' g'   Nr lr Nc lc L S J T g   JT-RME" << std::endl
        << " " << version << std::endl
        << " " << parameters.J0 << " " << parameters.g0
        << " " << parameters.T0_min << " " << parameters.T0_max
        << " " << int(parameters.symmetry_phase_mode) << std::endl
        << " " << parameters.Nmax << std::endl;
    }
    else  // if (io_mode == mcutils::IOMode::kBinary)
    {
      // field 1: version
      mcutils::WriteBinary<int32_t>(os, version);

      // fields 2-6: operator tensor properties
      mcutils::WriteBinary<int32_t>(os, parameters.J0);
      mcutils::WriteBinary<int32_t>(os, parameters.g0);
      mcutils::WriteBinary<int32_t>(os, parameters.T0_min);
      mcutils::WriteBinary<int32_t>(os, parameters.T0_max);
      assert(parameters.symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian);
      mcutils::WriteBinary<int32_t>(os, static_cast<int32_t>(parameters.symmetry_phase_mode));

      // fields 7: relative basis truncation
      mcutils::WriteBinary<int32_t>(os, parameters.Nmax);
    }
  }

  void ReadRelativeCMOperatorComponentLSJT(
      std::istream& is,
      int T0,
      const basis::RelativeCMSectorsLSJT& sectors,
      basis::OperatorBlocks<double>& matrices,
      mcutils::IOMode io_mode
    )
  {
    std::string line;
    int line_count = 0;

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {
        // extract sector
        const auto& sector = sectors.GetSector(sector_index);
        const auto& bra_subspace = sector.bra_subspace();
        const auto& ket_subspace = sector.ket_subspace();

        // verify that sector is canonical
        //
        // This is a check that the caller's sector construction
        // followed the specification that only "upper triangle"
        // sectors are stored.
        assert(sector.bra_subspace_index()<=sector.ket_subspace_index());

        // construct zero matrix for sector
        basis::OperatorBlock<double> sector_matrix =
          basis::OperatorBlock<double>::Zero(bra_subspace.size(),ket_subspace.size());

        std::vector<double> buffer;
        if (io_mode == mcutils::IOMode::kBinary)
          {
                  // calculate number of matrix elements in sector
            std::size_t sector_entries = 0;
            if (sector.IsDiagonal())
              // diagonal sector
              {
                std::size_t dimension = ket_subspace.size();
                sector_entries = dimension*(dimension+1)/2;
              }
            else  // if (sector.IsUpperTriangle())
              // upper triangle sector (but not diagonal)
              {
                std::size_t bra_dimension = bra_subspace.size();
                std::size_t ket_dimension = ket_subspace.size();
                sector_entries = bra_dimension*ket_dimension;
              }

            // read entire sector to temporary buffer
            buffer.resize(sector_entries, 0.);
            mcutils::ReadBinary<double>(is, buffer.data(), sector_entries);

          }

        std::size_t i = 0;
        // iterate over matrix elements
        for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
          for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
            {

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              double matrix_element;
              if (io_mode == mcutils::IOMode::kText)
              {
                // define states
                const basis::RelativeCMStateLSJT bra(bra_subspace,bra_index);
                const basis::RelativeCMStateLSJT ket(ket_subspace,ket_index);

                // define input variables
                int input_T0;
                int bra_Nr, bra_lr, bra_Nc, bra_lc, bra_L, bra_S, bra_J, bra_T, bra_g;
                int ket_Nr, ket_lr, ket_Nc, ket_lc, ket_L, ket_S, ket_J, ket_T, ket_g;

                // read input line
                mcutils::GetLine(is, line, line_count);
                std::istringstream line_stream(line);
                line_stream
                  >> input_T0
                  >> bra_Nr >> bra_lr >> bra_Nc >> bra_lc
                  >> bra_L  >> bra_S  >> bra_J  >> bra_T  >> bra_g
                  >> ket_Nr >> ket_lr >> ket_Nc >> ket_lc
                  >> ket_L  >> ket_S  >> ket_J  >> ket_T  >> ket_g
                  >> matrix_element;
                mcutils::ParsingCheck(line_stream, line_count, line);

                // validate labels
                bool expected_labels = true
                  && (input_T0==T0)
                  && (bra_Nr==bra.Nr())
                  && (bra_lr==bra.lr())
                  && (bra_Nc==bra.Nc())
                  && (bra_lc==bra.lc())
                  && (bra_L==bra.L())
                  && (bra_S==bra.S())
                  && (bra_J==bra.J())
                  && (bra_T==bra.T())
                  && (ket_Nr==ket.Nr())
                  && (ket_lr==ket.lr())
                  && (ket_Nc==ket.Nc())
                  && (ket_lc==ket.lc())
                  && (ket_L==ket.L())
                  && (ket_S==ket.S())
                  && (ket_J==ket.J())
                  && (ket_T==ket.T());
                assert(expected_labels);
              }
              else
              {
                matrix_element = buffer[i++];
              }

              // save matrix element
              sector_matrix(bra_index,ket_index) = matrix_element;
            }

          // store matrix for sector
          matrices.push_back(std::move(sector_matrix));
      }
  }

  void WriteRelativeCMOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const basis::RelativeCMSectorsLSJT& sectors,
      const basis::OperatorBlocks<double>& matrices,
      mcutils::IOMode io_mode
    )
  {

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
        const typename RelativeCMSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
        const typename RelativeCMSectorsLSJT::SubspaceType& bra_subspace = sector.bra_subspace();
        const typename RelativeCMSectorsLSJT::SubspaceType& ket_subspace = sector.ket_subspace();

        // verify that sector is canonical
        //
        // This is a check that the caller's sector construction
        // followed the specification that only "upper triangle"
        // sectors are stored.
        assert(sector.bra_subspace_index()<=sector.ket_subspace_index());

        // temporary buffer for binary write
        std::vector<double> buffer;
        std::size_t sector_entries = 0;
        if (io_mode == mcutils::IOMode::kBinary)
          {
            // calculate number of matrix elements in sector
            if (sector.IsDiagonal())
              // diagonal sector
              {
                const std::size_t& dimension = ket_subspace.size();
                sector_entries = dimension*(dimension+1)/2;
              }
            else  // if (sector.IsUpperTriangle())
              // upper triangle sector (but not diagonal)
              {
                const std::size_t& bra_dimension = bra_subspace.size();
                const std::size_t& ket_dimension = ket_subspace.size();
                sector_entries = bra_dimension*ket_dimension;
              }

            // allocate buffer
            buffer.resize(sector_entries, 0.);
          }

        // iterate over matrix elements
        std::size_t i = 0;
        for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
          for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
            {

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // define states
              const basis::RelativeCMStateLSJT bra(bra_subspace,bra_index);
              const basis::RelativeCMStateLSJT ket(ket_subspace,ket_index);

              // extract matrix element
              const double matrix_element = matrices[sector_index](bra_index,ket_index);

              if (io_mode == mcutils::IOMode::kText)
                {
                  // generate output line
                  const int width = 3;
                  const int precision = 8;  // for approximately single precision output
                  os << std::setprecision(precision);
                  os
                    << " " << std::setw(width) << T0
                    << " " << "  "
                    << " " << std::setw(width) << bra.Nr()
                    << " " << std::setw(width) << bra.lr()
                    << " " << std::setw(width) << bra.Nc()
                    << " " << std::setw(width) << bra.lc()
                    << " " << std::setw(width) << bra.L()
                    << " " << std::setw(width) << bra.S()
                    << " " << std::setw(width) << bra.J()
                    << " " << std::setw(width) << bra.T()
                    << " " << std::setw(width) << bra.g()
                    << " " << "    "
                    << " " << std::setw(width) << ket.Nr()
                    << " " << std::setw(width) << ket.lr()
                    << " " << std::setw(width) << ket.Nc()
                    << " " << std::setw(width) << ket.lc()
                    << " " << std::setw(width) << ket.L()
                    << " " << std::setw(width) << ket.S()
                    << " " << std::setw(width) << ket.J()
                    << " " << std::setw(width) << ket.T()
                    << " " << std::setw(width) << ket.g()
                    << " " << "    "
                    << " " << std::showpoint << std::scientific << matrix_element
                    << std::endl;
                }
              else  // if (io_mode == mcutils::IOMode::kText)
                {
                  buffer[i++] = matrix_element;
                }

            }

        // write temporary buffer to file
        if (io_mode == mcutils::IOMode::kBinary)
          mcutils::WriteBinary<double>(os,buffer.data(),sector_entries);
      }
  }


  void ReadRelativeCMOperatorLSJT(
      const std::string& filename,
      basis::RelativeCMSpaceLSJT& space,
      basis::RelativeCMOperatorParametersLSJT& operator_parameters,
      std::array<basis::RelativeCMSectorsLSJT,3>& component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& component_matrices,
      bool verbose
    )
  // FUTURE: change to line-based input with std::getline and parsing checks for easier debugging
  // FUTURE: check file status on open
  {

    // deduce I/O mode
    mcutils::IOMode io_mode = mcutils::DeducedIOMode(filename);

    // open stream for reading
    if (verbose)
      {
        std::cout
          << "Reading relative-cm operator file..." << std::endl
          << "  Filename: " << filename << std::endl
          << "  Mode: " << ((io_mode == mcutils::IOMode::kText) ? "text" : "binary")
          << std::endl;
      }
    std::ios_base::openmode mode_argument = std::ios_base::in;
    if (io_mode==mcutils::IOMode::kBinary)
      mode_argument |= std::ios_base::binary;
    std::ifstream is(filename.c_str(), mode_argument);

    // read header parameters
    basis::ReadRelativeCMOperatorParametersLSJT(is, operator_parameters, io_mode);
    if (verbose)
      {
        std::cout
          << "  Operator properties:"
          << " J0 " << operator_parameters.J0
          << " g0 " << operator_parameters.g0
          << " T0_min " << operator_parameters.T0_min
          << " T0_max " << operator_parameters.T0_max
          << " symmetry " << int(operator_parameters.symmetry_phase_mode)
          << std::endl
          << "  Truncation:"
          << " Nmax " << operator_parameters.Nmax
          << std::endl;
      }

    // set up relative space
    space = basis::RelativeCMSpaceLSJT(
        operator_parameters.Nmax
      );

    // populate sectors and matrices
    for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
      // for each isospin component
      {
        // enumerate sectors
        component_sectors[T0]
          = basis::RelativeCMSectorsLSJT(space,operator_parameters.J0,T0,operator_parameters.g0);

        // read matrices
        basis::ReadRelativeCMOperatorComponentLSJT(
            is,
            T0,
            component_sectors[T0], component_matrices[T0],
            io_mode
          );
      }

    // write diagnostics
    if (verbose)
      {
        std::cout << "  Matrix elements:";
        for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
          std::cout << " " << basis::UpperTriangularEntries(component_sectors[T0]);
        std::cout << std::endl;
        std::cout << "  Allocated:";
        for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
          std::cout << " " << basis::AllocatedEntries(component_matrices[T0]);
        std::cout << std::endl;
      }
  }

  void WriteRelativeCMOperatorLSJT(
      const std::string& filename,
      const basis::RelativeCMSpaceLSJT& space,
      const basis::OperatorLabelsJT& operator_labels,
      const std::array<basis::RelativeCMSectorsLSJT,3>& component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& component_matrices,
      bool verbose
    )
  // FUTURE: check file status on open
  {

    // deduce I/O mode
    mcutils::IOMode io_mode = mcutils::DeducedIOMode(filename);

    // open stream for writing
    if (verbose)
      {
        std::cout
          << "Writing relative-cm operator file..." << std::endl
          << "  Filename: " << filename << std::endl
          << "  Mode: " << ((io_mode == mcutils::IOMode::kText) ? "text" : "binary")
          << std::endl;
      }
    std::ios_base::openmode mode_argument = std::ios_base::out;
    if (io_mode==mcutils::IOMode::kBinary)
      mode_argument |= std::ios_base::binary;
    std::ofstream os(filename.c_str(), mode_argument);

    // set up full operator file parameters
    int Nmax = space.Nmax();
    RelativeCMOperatorParametersLSJT operator_parameters(operator_labels,Nmax);
    if (verbose)
      {
        std::cout
          << "  Operator properties:"
          << " J0 " << operator_parameters.J0
          << " g0 " << operator_parameters.g0
          << " T0_min " << operator_parameters.T0_min
          << " T0_max " << operator_parameters.T0_max
          << " symmetry " << int(operator_parameters.symmetry_phase_mode)
          << std::endl
          << "  Truncation:"
          << " Nmax " << operator_parameters.Nmax
          << std::endl;
      }

    // write header parameters
    basis::WriteRelativeCMOperatorParametersLSJT(os,operator_parameters,io_mode);

    // write matrices
    for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
      {
        basis::WriteRelativeCMOperatorComponentLSJT(
            os,
            T0,
            component_sectors[T0], component_matrices[T0],
            io_mode
          );
      }

  }

  ////////////////////////////////////////////////////////////////
  // relative LSJT operator construction
  ////////////////////////////////////////////////////////////////

  void ConstructZeroOperatorRelativeCMLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceLSJT& relative_cm_space,
      std::array<basis::RelativeCMSectorsLSJT,3>& relative_cm_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_cm_component_matrices
    )
  {

    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate sectors
        relative_cm_component_sectors[T0]
          = basis::RelativeCMSectorsLSJT(relative_cm_space,operator_labels.J0,T0,operator_labels.g0);

        // populate matrices
        relative_cm_component_matrices[T0].resize(relative_cm_component_sectors[T0].size());
        basis::SetOperatorToZero(relative_cm_component_sectors[T0],relative_cm_component_matrices[T0]);
      }
  }

  void ConstructIdentityOperatorRelativeCMLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceLSJT& relative_cm_space,
      std::array<basis::RelativeCMSectorsLSJT,3>& relative_cm_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_cm_component_matrices
    )
  {

    // validate operator labels
    assert(operator_labels.J0==0);
    assert(operator_labels.g0==0);
    assert(operator_labels.symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian);

    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate sectors
        relative_cm_component_sectors[T0]
          = basis::RelativeCMSectorsLSJT(relative_cm_space,operator_labels.J0,T0,operator_labels.g0);

        // populate matrices
        relative_cm_component_matrices[T0].resize(relative_cm_component_sectors[T0].size());
        if (T0==0)
          // identity matrices in
          basis::SetOperatorToIdentity(relative_cm_component_sectors[T0],relative_cm_component_matrices[T0]);
        else
          basis::SetOperatorToZero(relative_cm_component_sectors[T0],relative_cm_component_matrices[T0]);
      }
  }

  ////////////////////////////////////////////////////////////////
  // two-body LSJT operator output
  ////////////////////////////////////////////////////////////////

  void WriteTwoBodyOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const TwoBodySectorsLSJT& sectors,
      const OperatorBlocks<double>& matrices,
      basis::NormalizationConversion conversion_mode
    )
  {

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
        const typename TwoBodySectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
        const typename TwoBodySectorsLSJT::SubspaceType& bra_subspace = sector.bra_subspace();
        const typename TwoBodySectorsLSJT::SubspaceType& ket_subspace = sector.ket_subspace();

        // verify that sector is canonical
        //
        // This is a check that the caller's sector construction
        // followed the specification that only "upper triangle"
        // sectors are stored.
        assert(sector.bra_subspace_index()<=sector.ket_subspace_index());

        // iterate over matrix elements
        for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
          for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
            {

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // define states
              const basis::TwoBodyStateLSJT bra(bra_subspace,bra_index);
              const basis::TwoBodyStateLSJT ket(ket_subspace,ket_index);

              // determine matrix element normalization factor
              double conversion_factor = 1.;
              if (conversion_mode == basis::NormalizationConversion::kASToNAS)
                {
                  if ((bra.N1()==bra.N2())&&(bra.l1()==bra.l2()))
                    conversion_factor *= (1/sqrt(2.));
                  if ((ket.N1()==ket.N2())&&(ket.l1()==ket.l2()))
                    conversion_factor *= (1/sqrt(2.));
                }
              else if (conversion_mode == basis::NormalizationConversion::kNASToAS)
                {
                  if ((bra.N1()==bra.N2())&&(bra.l1()==bra.l2()))
                    conversion_factor *= sqrt(2.);
                  if ((ket.N1()==ket.N2())&&(ket.l1()==ket.l2()))
                    conversion_factor *= sqrt(2.);
                }

              // extract matrix element
              const double matrix_element = conversion_factor*matrices[sector_index](bra_index,ket_index);

              // generate output line
              const int width = 3;
              const int precision = 8;  // for approximately single precision output
              os << std::setprecision(precision);
              os
                << " " << std::setw(width) << T0
                << " " << "  "
                << " " << std::setw(width) << bra.N1()
                << " " << std::setw(width) << bra.l1()
                << " " << std::setw(width) << bra.N2()
                << " " << std::setw(width) << bra.l2()
                << " " << std::setw(width) << bra.L()
                << " " << std::setw(width) << bra.S()
                << " " << std::setw(width) << bra.J()
                << " " << std::setw(width) << bra.T()
                << " " << std::setw(width) << bra.g()
                << " " << "    "
                << " " << std::setw(width) << ket.N1()
                << " " << std::setw(width) << ket.l1()
                << " " << std::setw(width) << ket.N2()
                << " " << std::setw(width) << ket.l2()
                << " " << std::setw(width) << ket.L()
                << " " << std::setw(width) << ket.S()
                << " " << std::setw(width) << ket.J()
                << " " << std::setw(width) << ket.T()
                << " " << std::setw(width) << ket.g()
                << " " << "    "
                << " " << std::showpoint << std::scientific << matrix_element
                << std::endl;

            }

      };
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
