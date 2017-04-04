/****************************************************************

  lsjt_operator.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/

#include "lsjt_operator.h"

#include <fstream>

namespace basis {

  ////////////////////////////////////////////////////////////////
  // relative two-body operator
  ////////////////////////////////////////////////////////////////

  void WriteRelativeOperatorParametersLSJT(
      std::ostream& os,
      const basis::RelativeOperatorParametersLSJT& parameters
    )
  {
    int version = 1;
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

  void ReadRelativeOperatorParametersLSJT(
      std::istream& is,
      basis::RelativeOperatorParametersLSJT& parameters
    )
  {

    std::string line;

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

  void WriteRelativeOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const RelativeSectorsLSJT& sectors,
      const MatrixVector& matrices
    )
  {

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
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

	// iterate over matrix elements
	for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
	  for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
	    {

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // define states
	      const basis::RelativeStateLSJT bra(bra_subspace,bra_index);
	      const basis::RelativeStateLSJT ket(ket_subspace,ket_index);

              // extract matrix element factor
      	      const double matrix_element = matrices[sector_index](bra_index,ket_index);
	    
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

      };
  }

  void ReadRelativeOperatorComponentLSJT(
      std::istream& is,
      int T0,
      const RelativeSectorsLSJT& sectors,
      MatrixVector& matrices
    )
  {

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const typename RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
	const typename RelativeSectorsLSJT::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename RelativeSectorsLSJT::SubspaceType& ket_subspace = sector.ket_subspace();

        // construct zero matrix for sector
	Eigen::MatrixXd sector_matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());

	// retrieve matrix elements
	for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
	  for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
	    {

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

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
	      double matrix_element;
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

              // save matrix element
	      sector_matrix(bra_index,ket_index) = matrix_element;
	    }

        // store matrix for sector
	matrices.push_back(sector_matrix);

      }
  }

  void ReadRelativeOperatorLSJT(
      const std::string& relative_filename,
      basis::RelativeSpaceLSJT& relative_space,
      basis::OperatorLabelsJT& operator_labels,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      bool verbose
    )
  // FUTURE: change to line-based input with std::getline and parsing checks for easier debugging
  // FUTURE: check file status on open
  {

    // open stream for reading
    if (verbose)
      {
        std::cout
          << "Reading relative operator file..." << std::endl
          << "  Filename: " << relative_filename << std::endl;
      }
    std::ifstream is(relative_filename.c_str());

    // read header parameters
    basis::RelativeOperatorParametersLSJT operator_parameters;
    basis::ReadRelativeOperatorParametersLSJT(is,operator_parameters);
    operator_labels = static_cast<basis::OperatorLabelsJT>(operator_parameters);
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
          = basis::RelativeSectorsLSJT(relative_space,operator_labels.J0,T0,operator_labels.g0);

        // read matrices
        basis::ReadRelativeOperatorComponentLSJT(
            is,
            T0,
            relative_component_sectors[T0],relative_component_matrices[T0]
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
      const std::array<basis::MatrixVector,3>& relative_component_matrices,
      bool verbose
    )
  // FUTURE: check file status on open
  {

    // open stream for writing
    if (verbose)
      {
        std::cout
          << "Writing relative operator file..." << std::endl
          << "  Filename: " << relative_filename << std::endl;
      }
    std::ofstream os(relative_filename.c_str());

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
    basis::WriteRelativeOperatorParametersLSJT(os,operator_parameters);

    // write matrices
    for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
      {
        basis::WriteRelativeOperatorComponentLSJT(
            os,
            T0,
            relative_component_sectors[T0],relative_component_matrices[T0]
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
      std::array<basis::MatrixVector,3>& relative_component_matrices
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
        basis::SetOperatorToZero(relative_component_sectors[T0],relative_component_matrices[T0]);
      }
  }

  void ConstructIdentityOperatorRelativeLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices
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


  void GatherOperatorRelativeCMLSJTNToRelativeCMLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
      const std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
      const std::array<basis::MatrixVector,3>& relative_cm_lsjtn_component_matrices,
      const basis::RelativeCMSpaceLSJT& relative_cm_lsjt_space,
      std::array<basis::RelativeCMSectorsLSJT,3>& relative_cm_lsjt_component_sectors,
      std::array<basis::MatrixVector,3>& relative_cm_lsjt_component_matrices
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
        for (int sector_index=0; sector_index<relative_cm_lsjt_component_sectors[T0].size(); ++sector_index)
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
            for (int bra_index = 0; bra_index < relative_cm_lsjt_sector.bra_subspace().size(); ++bra_index)
              for (int ket_index = 0; ket_index < relative_cm_lsjt_sector.ket_subspace().size(); ++ket_index)
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
                  int relative_cm_lsjtn_subspace_index_bra
                    = relative_cm_lsjtn_space.LookUpSubspaceIndex(
                        relative_cm_lsjtn_subspace_labels_bra
                      );
                  int relative_cm_lsjtn_state_index_bra
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
                  int relative_cm_lsjtn_subspace_index_ket
                    = relative_cm_lsjtn_space.LookUpSubspaceIndex(
                        relative_cm_lsjtn_subspace_labels_ket
                      );
                  int relative_cm_lsjtn_state_index_ket
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
                  int relative_cm_lsjtn_sector_index
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
      const std::array<basis::MatrixVector,3>& two_body_lsjtn_component_matrices,
      const basis::TwoBodySpaceLSJT& two_body_lsjt_space,
      std::array<basis::TwoBodySectorsLSJT,3>& two_body_lsjt_component_sectors,
      std::array<basis::MatrixVector,3>& two_body_lsjt_component_matrices
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
        for (int sector_index=0; sector_index<two_body_lsjt_component_sectors[T0].size(); ++sector_index)
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
            for (int bra_index = 0; bra_index < two_body_lsjt_sector.bra_subspace().size(); ++bra_index)
              for (int ket_index = 0; ket_index < two_body_lsjt_sector.ket_subspace().size(); ++ket_index)
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
                  int two_body_lsjtn_subspace_index_bra
                    = two_body_lsjtn_space.LookUpSubspaceIndex(
                        two_body_lsjtn_subspace_labels_bra
                      );
                  int two_body_lsjtn_state_index_bra
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
                  int two_body_lsjtn_subspace_index_ket
                    = two_body_lsjtn_space.LookUpSubspaceIndex(
                        two_body_lsjtn_subspace_labels_ket
                      );
                  int two_body_lsjtn_state_index_ket
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
                  int two_body_lsjtn_sector_index
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
  // relative-cm LSJT operator output
  ////////////////////////////////////////////////////////////////

  void WriteRelativeCMOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const basis::RelativeCMSectorsLSJT& sectors,
      const basis::MatrixVector& matrices
    )
  {

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
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

	// iterate over matrix elements
	for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
	  for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
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

      };
  }


  ////////////////////////////////////////////////////////////////
  // two-body LSJT operator output
  ////////////////////////////////////////////////////////////////

  void WriteTwoBodyOperatorComponentLSJT(
      std::ostream& os,
      int T0,
      const TwoBodySectorsLSJT& sectors,
      const MatrixVector& matrices,
      basis::NormalizationConversion conversion_mode
    )
  {

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
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
	for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
	  for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
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
