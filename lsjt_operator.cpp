/****************************************************************

  lsjt_operator.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/

#include "lsjt_operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // relative two-body operator
  ////////////////////////////////////////////////////////////////

  void WriteRelativeOperatorParametersLSJT(
      std::ostream& os,
      const basis::RelativeOperatorParametersLSJT& parameters
    )
  {
    assert(parameters.version==1);
    os 
      << "# RELATIVE LSJT" << std::endl
      << "#   version" << std::endl
      << "#   J0 g0 T0_min T0_max symmetry_phase_mode  [P0=(-)^g0]" << std::endl
      << "#   Nmax Jmax" << std::endl
      << "#   T0   N' L' S' J' T'   N L S J T   JT-RME" << std::endl
      << " " << parameters.version << std::endl
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

    // line 1: version -- and gobble any comment lines
    while (std::getline(is,line), line[0]=='#') {};
    std::stringstream(line) >> parameters.version;
    assert(parameters.version==1);
    
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
              const int precision = 16;
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

  ////////////////////////////////////////////////////////////////
  // operator matrix element canonicalizaton
  ////////////////////////////////////////////////////////////////

  void CanonicalizeIndicesRelativeLSJT(
      const basis::RelativeSpaceLSJT& space,
      int& bra_subspace_index, int& ket_subspace_index,
      int& bra_state_index, int& ket_state_index,
      double& canonicalization_factor,
      int J0, int T0, int g0,
      basis::SymmetryPhaseMode symmetry_phase_mode
    )
  // DEPRECATED
  {

    // canonicalize indices
    bool swapped_subspaces, swapped_states;
    basis::CanonicalizeIndices(
        bra_subspace_index, ket_subspace_index,
        swapped_subspaces,
        bra_state_index, ket_state_index,
        swapped_states
      );

    // calculate canonicalization factor
    //
    // Beware that the indices now describe the "new" bra and ket
    // *after* any swap, so one must take care in matching up bra and
    // ket labels to those in any formula describing the symmetry.

    // check that case is covered
    //
    // Phase definitions are currently only provided for
    // Hamiltonian-like operators.
    assert(
        (symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian)
        && (J0==0) && (g0==0)
      );
    
    // case: Hamiltonian-like operator
    //
    // Recall symmetry relation:
    //
    //     <a,J,T,g || A_{T0} || a',J,T',g>
    //       = (-)^(T'-T)*Hat(T')/Hat(T)
    //         * <a',J,T',g || A_{T0} || a,J,T,g>
    
    canonicalization_factor = 1.;
    if (swapped_subspaces)
      {

        // retrieve sector labels (*after* swap, i.e., canonical m.e. on RHS)
        const basis::RelativeSubspaceLSJT& bra_subspace = space.GetSubspace(
            bra_subspace_index
          );
        const basis::RelativeSubspaceLSJT& ket_subspace = space.GetSubspace(
            ket_subspace_index
          );
        int Tp = bra_subspace.T();
        int T = ket_subspace.T();

        canonicalization_factor *= ParitySign(Tp-T)*Hat(Tp)/Hat(T);
      }
  }

  ////////////////////////////////////////////////////////////////
  // relative LSJT operator construction
  ////////////////////////////////////////////////////////////////

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
  // two-body LSJT operator manipulation
  ////////////////////////////////////////////////////////////////

  inline
  void RecastLabelsTwoBodyLSJTToTwoBodyNLSJT(
      const TwoBodySubspaceLSJTLabels& two_body_lsjt_subspace_labels,
      const TwoBodyStateLSJTLabels& two_body_lsjt_state_labels,
      TwoBodySubspaceNLSJTLabels& two_body_nlsjt_subspace_labels,
      TwoBodyStateNLSJTLabels& two_body_nlsjt_state_labels
    )
  // Recast labels for (subspace,state) from TwoBodyLSJT scheme to
  // TwoBodyNLSJT scheme.
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
  //   two_body_nlsjt_subspace_labels (...,output) : target subspace labels
  //   two_body_nlsjt_state_labels (...,output) : target state labels
  {
    // extract labels
    int L, S, J, T, g;
    std::tie(L,S,J,T,g) = two_body_lsjt_subspace_labels;
    int N1, l1, N2, l2;
    std::tie(N1,l1,N2,l2) = two_body_lsjt_state_labels;

    // repackage labels
    int N = N1+N2;
    two_body_nlsjt_subspace_labels = TwoBodySubspaceNLSJTLabels(L,S,J,T,g,N);
    two_body_nlsjt_state_labels = two_body_lsjt_state_labels;
  }

  void GatherOperatorTwoBodyNLSJTToTwoBodyLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceNLSJT& two_body_nlsjt_space,
      const std::array<basis::TwoBodySectorsNLSJT,3>& two_body_nlsjt_component_sectors,
      const std::array<basis::MatrixVector,3>& two_body_nlsjt_component_matrices,
      const basis::TwoBodySpaceLSJT& two_body_lsjt_space,
      std::array<basis::TwoBodySectorsLSJT,3>& two_body_lsjt_component_sectors,
      std::array<basis::MatrixVector,3>& two_body_lsjt_component_matrices
    )
  {
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    // for each isospin component
    {

      // enumerate sectors
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
                
                // // look up source subspace indices
                // int two_body_nlsjt_bra_subspace_index = two_body_nlsjt_space.LookUpSubspaceIndex(
                //     basis::TwoBodySubspaceNLSJTLabels(
                //         two_body_lsjt_bra.L(),
                //         two_body_lsjt_bra.S(),
                //         two_body_lsjt_bra.J(),
                //         two_body_lsjt_bra.T(),
                //         two_body_lsjt_bra.g(),
                //         two_body_lsjt_bra.N()
                //       )
                //   );
                // int two_body_nlsjt_ket_subspace_index = two_body_nlsjt_space.LookUpSubspaceIndex(
                //     basis::TwoBodySubspaceNLSJTLabels(
                //         two_body_lsjt_ket.L(),
                //         two_body_lsjt_ket.S(),
                //         two_body_lsjt_ket.J(),
                //         two_body_lsjt_ket.T(),
                //         two_body_lsjt_ket.g(),
                //         two_body_lsjt_ket.N()
                //       )
                //   );
                // 
                // // look up source matrix element indices
                // const basis::TwoBodySubspaceNLSJT& two_body_nlsjt_bra_subspace
                //   = two_body_nlsjt_space.GetSubspace(two_body_nlsjt_bra_subspace_index);
                // const basis::TwoBodySubspaceNLSJT& two_body_nlsjt_ket_subspace
                //   = two_body_nlsjt_space.GetSubspace(two_body_nlsjt_ket_subspace_index);
                // int two_body_nlsjt_bra_index = two_body_nlsjt_bra_subspace.LookUpStateIndex(
                //     basis::TwoBodyStateNLSJT::StateLabelsType(
                //         two_body_lsjt_bra.N1(),
                //         two_body_lsjt_bra.l1(),
                //         two_body_lsjt_bra.N2(),
                //         two_body_lsjt_bra.l2()
                //       )
                //   );
                // int two_body_nlsjt_ket_index = two_body_nlsjt_ket_subspace.LookUpStateIndex(
                //     basis::TwoBodyStateNLSJT::StateLabelsType(
                //         two_body_lsjt_ket.N1(),
                //         two_body_lsjt_ket.l1(),
                //         two_body_lsjt_ket.N2(),
                //         two_body_lsjt_ket.l2()
                //       )
                //   );

                // extract source bra labels
                TwoBodySubspaceLSJTLabels two_body_lsjt_subspace_labels_bra
                  = two_body_lsjt_bra.subspace().labels();
                TwoBodyStateLSJTLabels two_body_lsjt_state_labels_bra
                  = two_body_lsjt_bra.labels();
                TwoBodySubspaceNLSJTLabels two_body_nlsjt_subspace_labels_bra;
                TwoBodyStateNLSJTLabels two_body_nlsjt_state_labels_bra;
                RecastLabelsTwoBodyLSJTToTwoBodyNLSJT(
                    two_body_lsjt_subspace_labels_bra,
                    two_body_lsjt_state_labels_bra,
                    two_body_nlsjt_subspace_labels_bra,
                    two_body_nlsjt_state_labels_bra
                  );

                // extract source bra indices
                int two_body_nlsjt_subspace_index_bra
                  = two_body_nlsjt_space.LookUpSubspaceIndex(
                      two_body_nlsjt_subspace_labels_bra
                    );
                int two_body_nlsjt_state_index_bra
                  = two_body_nlsjt_space.GetSubspace(two_body_nlsjt_subspace_index_bra).LookUpStateIndex(
                      two_body_nlsjt_state_labels_bra
                    );

                // extract source ket labels
                TwoBodySubspaceLSJTLabels two_body_lsjt_subspace_labels_ket
                  = two_body_lsjt_ket.subspace().labels();
                TwoBodyStateLSJTLabels two_body_lsjt_state_labels_ket
                  = two_body_lsjt_ket.labels();
                TwoBodySubspaceNLSJTLabels two_body_nlsjt_subspace_labels_ket;
                TwoBodyStateNLSJTLabels two_body_nlsjt_state_labels_ket;
                RecastLabelsTwoBodyLSJTToTwoBodyNLSJT(
                    two_body_lsjt_subspace_labels_ket,
                    two_body_lsjt_state_labels_ket,
                    two_body_nlsjt_subspace_labels_ket,
                    two_body_nlsjt_state_labels_ket
                  );

                // extract source ket indices
                int two_body_nlsjt_subspace_index_ket
                  = two_body_nlsjt_space.LookUpSubspaceIndex(
                      two_body_nlsjt_subspace_labels_ket
                    );
                int two_body_nlsjt_state_index_ket
                  = two_body_nlsjt_space.GetSubspace(two_body_nlsjt_subspace_index_ket).LookUpStateIndex(
                      two_body_nlsjt_state_labels_ket
                    );

                // // canonicalize indices for matrix element lookup
                // //
                // // Matrix elements which are canonical by
                // // LSJT sector, and ordered by N within a LSJT sector,
                // // should also be canonical by (LSJT;N) sector.  That
                // // is, the upper triangle of a matrix with basis
                // // states ordered by
                // //
                // //   (L,S,J,T,g) : ([N],N1,l1,N2,l2)
                // //
                // // or
                // //
                // //   (L,S,J,T,g,N) : (N1,l1,N2,l2)
                // //
                // // should be identical.
                // //
                // // So canonicalization should only be necessary when
                // // we are filling in a *lower* triangle matrix element
                // // of a diagonal target sector.  Then we must ensure
                // // that we look up a canonical (upper triangular)
                // // NLSJT sector.
                // //
                // // If the matrix element happens to fall within a
                // // diagonal NLSJT subblocks, canonicalization of the
                // // state indices should not actually be necessary,
                // // since the Moshinsky transformation machinery up
                // // until this point has actually been populating the
                // // full (square) matrices for the diagonal sectors.
                // //
                // // Note that no canonicalization factor is needed.
                // // Since N is the least significant subspace label
                // // in the lexicographic ordering, canonical
                // // swaps will never entail swapping subspace LSJT
                // // labels, just the N labels.
                // //
                // // Maybe a saner alternative is just to look up upper
                // // triangular matrix elements, then brute force
                // // symmetrize all LSJT diagonal sectors.
                // 
                // bool swapped_subspaces, swapped_states;
                // basis::CanonicalizeIndices(
                //     two_body_nlsjt_subspace_index_bra,two_body_nlsjt_subspace_index_ket,
                //     swapped_subspaces,
                //     two_body_nlsjt_state_index_bra,two_body_nlsjt_state_index_ket,
                //     swapped_states
                //   );

                // look up matrix element
                int two_body_nlsjt_sector_index
                  = two_body_nlsjt_component_sectors[T0].LookUpSectorIndex(
                      two_body_nlsjt_subspace_index_bra,
                      two_body_nlsjt_subspace_index_ket
                    );

                const Eigen::MatrixXd& two_body_nlsjt_matrix
                  = two_body_nlsjt_component_matrices[T0][two_body_nlsjt_sector_index];
                double two_body_nlsjt_matrix_element = two_body_nlsjt_matrix(
                    two_body_nlsjt_state_index_bra,two_body_nlsjt_state_index_ket
                  );

                two_body_lsjt_matrix(bra_index,ket_index) = two_body_nlsjt_matrix_element;

              }
        }
    }

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

              // extract matrix element factor
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
              const double matrix_element = conversion_factor*matrices[sector_index](bra_index,ket_index);

              // generate output line
              const int width = 3;
              const int precision = 16;
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
