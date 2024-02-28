/****************************************************************
  oscillator_orbital_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <iomanip>
#include <iostream>
#include <limits>

#include "am/halfint.h"

#include "operator.h"
#include "oscillator_orbital.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void TestOscillatorOrbitalSubspace()
{
  // subspace construction
  //
  // Normally subspaces are always constructed as part of a space.  We would not
  // construct a standalone subspace.  But, for testing purposes, it is
  // convenient to do so, since the code for subspace and state can be written
  // and tested independently of the code for space and/or sectors.

  std::cout << "Subspace construction" << std::endl;
  std::cout << std::endl;

  const basis::OscillatorOrbitalSubspace subspace(2, 4);  // l=2, Nmax=4
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();

  // state construction
  //
  // Let us give the state type a more thorough workout.  Note that the state
  // type was already used in the implementation of the subspace's DebugStr()
  // above.

  std::cout << std::endl;
  std::cout << "State construction" << std::endl;
  std::cout << std::endl;

  // illustrate construction by index vs. state labels

  // construct by index -- invalid cases

  if (false)
    {
      // Index -1 should be trapped and lead to assertion failure.
      //
      // Note: Negative index cast to size_t will be large positive index, by
      // two-complement.

      std::size_t test_index = (std::size_t)(-1);
      std::cout << "Constructing with index " << test_index << std::endl;
      const basis::OscillatorOrbitalState state_from_index(subspace, test_index);
    }

  if (false)
    {
      // But a very large negative integer could wrap and become a valid index again!
      //
      // std::size_t test_index = (std::size_t)(-18446744073709551615);

      std::size_t max_size = std::numeric_limits<std::size_t>().max();
      std::size_t test_index = -max_size;  // should be equivalent to +1
      std::cout << "Constructing with index " << test_index << std::endl;
      const basis::OscillatorOrbitalState state_from_index(subspace, test_index);
    }

  // construct by index -- valid case
  std::size_t test_index = 1;
  std::cout << "Constructing with index " << test_index << std::endl;
  const basis::OscillatorOrbitalState state_from_index(subspace, test_index);
  std::cout << "  label string " << state_from_index.LabelStr() << std::endl;
  std::cout << "  quantum numbers"
            << " l " << state_from_index.l()
            << " g " << state_from_index.g()
            << " n " << state_from_index.n()
            << " N " << state_from_index.N()
            << std::endl;

  // construct by state labels
  int n = 0;
  std::cout << "Constructing with n " << n << std::endl;
  const basis::OscillatorOrbitalState state_from_labels(subspace, basis::OscillatorOrbitalState::StateLabelsType(n));
  std::cout << "  label string " << state_from_labels.LabelStr() << std::endl;
  std::cout << "  quantum numbers"
            << " l " << state_from_labels.l()
            << " g " << state_from_labels.g()
            << " n " << state_from_labels.n()
            << " N " << state_from_labels.N()
            << std::endl;

  // try out accessors
  std::cout << "Try out accessors" << std::endl;
  std::cout << "subspace " << state_from_labels.subspace().LabelStr() << std::endl;
  std::cout << "labels " << std::get<0>(state_from_labels.labels()) << std::endl;
  std::cout << "index " << state_from_labels.index() << std::endl;

  // text equality operator
  std::cout << "  equality test " << (state_from_index == state_from_index) << std::endl;
  std::cout << "  equality test " << (state_from_index == state_from_labels) << std::endl;

  // iterate over states within subspace
  //
  // here we also try out the "state factory" member function of the subspace
  std::cout << std::endl;
  std::cout << "Iterate over states in subspace" << std::endl;
  std::cout << std::endl;

  // iterate by state index
  std::cout << "by index" << std::endl;
  for (std::size_t state_index=0; state_index<subspace.size(); ++state_index)
    {
      // const basis::OscillatorOrbitalState state(subspace,state_index);
      const basis::OscillatorOrbitalState state = subspace.GetState(state_index);
      std::cout << "index " << state.index() << " N " << state.N() << std::endl;
    };

  // iterate by labels
  std::cout << "by labels" << std::endl;
  int n_max = (subspace.Nmax()-subspace.l())/2;
  for (int n=0; n<=n_max; ++n)
    {
      const basis::OscillatorOrbitalState state = subspace.GetState(basis::OscillatorOrbitalState::StateLabelsType(n));
      std::cout << "index " << state.index() << " N " << state.N() << std::endl;
    };

  std::cout << std::endl;

}


void TestOscillatorOrbitalSpace()
{

  std::cout << "Space construction" << std::endl;
  std::cout << std::endl;

  // construct space
  int Nmax = 4;
  basis::OscillatorOrbitalSpace space(Nmax);

  // print diagnostics
  std::cout << space.DebugStr();

  std::cout << std::endl;

  // try out accessors
  std::cout << "Try out accessors" << std::endl;
  std::cout << "Nmax " << space.Nmax() << std::endl;
  std::cout << "ContainsSubspace " << space.ContainsSubspace(basis::OscillatorOrbitalSubspace::SubspaceLabelsType(5)) << std::endl;
  std::cout << "LookUpSubspaceIndex " << space.LookUpSubspaceIndex(basis::OscillatorOrbitalSubspace::SubspaceLabelsType(2)) << std::endl;
  std::cout << "LookUpSubspace " << space.LookUpSubspace(basis::OscillatorOrbitalSubspace::SubspaceLabelsType(2)).LabelStr() << std::endl;
  std::cout << "size " << space.size() << std::endl;
  std::cout << "dimension " << space.dimension() << std::endl;

  const basis::OscillatorOrbitalSubspace& subspace = space.LookUpSubspace(basis::OscillatorOrbitalSubspace::SubspaceLabelsType(2));
  std::cout << "LookUpSubspace with a reference variable " << subspace.LabelStr() << std::endl;
}

void TestOscillatorOrbitalSectors()
{

  std::cout << "Sectors construction" << std::endl;
  std::cout << std::endl;

  // construct space
  int Nmax = 4;
  basis::OscillatorOrbitalSpace space(Nmax);

  // construct sectors -- L0=0, g0=0 (Hamiltonian-like)
  basis::OscillatorOrbitalSectors hamiltonian_sectors(space, 0, 0);
  std::cout << "hamiltonian_sectors" << std::endl
            << hamiltonian_sectors.DebugStr()
            << std::endl;

  // find indices for a matrix element from labels
  std::cout << "Find indices for a matrix element from labels" << std::endl;
  std::size_t bra_subspace_index=space.LookUpSubspaceIndex(basis::OscillatorOrbitalSubspace::SubspaceLabelsType(2));
  std::size_t ket_subspace_index=space.LookUpSubspaceIndex(basis::OscillatorOrbitalSubspace::SubspaceLabelsType(2));
  std::size_t sector_index=hamiltonian_sectors.LookUpSectorIndex(bra_subspace_index,ket_subspace_index);
  std::cout << "sector index " << sector_index << std::endl;
  basis::OscillatorOrbitalSectors::SectorType sector=hamiltonian_sectors.GetSector(sector_index);
  // const basis::OscillatorOrbitalSectors::SectorType::BraSubspaceType& bra_subspace=sector.bra_subspace();
  // const basis::OscillatorOrbitalSectors::SectorType::KetSubspaceType& ket_subspace=sector.ket_subspace();
  const basis::OscillatorOrbitalSubspace& bra_subspace=sector.bra_subspace();
  const basis::OscillatorOrbitalSubspace& ket_subspace=sector.ket_subspace();
  // auto& bra_subspace=sector.bra_subspace();
  // auto& ket_subspace=sector.ket_subspace();
  std::cout << "bra_subspace " << bra_subspace.LabelStr() << std::endl;
  std::cout << "ket_subspace " << ket_subspace.LabelStr() << std::endl;
  std::cout << "bra_subspace " << std::endl << bra_subspace.DebugStr() << std::endl;
  std::cout << "ket_subspace " << std::endl << ket_subspace.DebugStr() << std::endl;
  std::size_t bra_state_index=bra_subspace.LookUpStateIndex(basis::OscillatorOrbitalState::StateLabelsType(0));
  std::size_t ket_state_index=ket_subspace.LookUpStateIndex(basis::OscillatorOrbitalState::StateLabelsType(1));
  std::cout << "bra state index " << bra_state_index << std::endl;
  std::cout << "ket state index " << ket_state_index << std::endl;

  // construct sectors -- L0=1, g0=1 (E1-like)
  basis::OscillatorOrbitalSectors e1_sectors(space, 1, 1);
  std::cout << "e1_sectors" << std::endl
            << e1_sectors.DebugStr()
            << std::endl;

  // construct sectors -- L0=1, g0=1 (E1-like) -- but including noncanonical ("lower-triangle") sectors
  basis::OscillatorOrbitalSectors e1_sectors_both_ways(
      space, 1, 1,
      basis::SectorDirection::kBoth
    );
  std::cout << "e1_sectors_both_ways" << std::endl
            << e1_sectors_both_ways.DebugStr()
            << std::endl;


  // construct sectors -- L0=2, g0=0 (E2-like)
  basis::OscillatorOrbitalSectors e2_sectors(space, 2, 0);
  std::cout << "e2_sectors" << std::endl
            << e2_sectors.DebugStr()
            << std::endl;

  // construct sectors -- L0=2, g0=1 (M2-like)
  basis::OscillatorOrbitalSectors m2_sectors(space, 2, 1);
  std::cout << "n2_sectors" << std::endl
            << m2_sectors.DebugStr()
            << std::endl;

}

////////////////////////////////////////////////////////////////
// ladder operator implementation and test
////////////////////////////////////////////////////////////////

// enumerated type for ladder operator
enum class LadderOperatorType {
  kRaising,
  kLowering,
};

// ladder operator implementation

void LadderOperator(
    LadderOperatorType operator_type,
    const basis::OscillatorOrbitalSpace& space,
    const basis::OscillatorOrbitalSectors& sectors,
    basis::OperatorBlocks<double>& matrices
  )
// Generate reduced matrix elements of harmonic oscillator ladder operator
// $c^\dagger$ or $\tilde{c}$.
//
// TODO: Resolve $\tilde{c}$ vs. plain $c$ and radial convention.
//
// RMEs are in Rose convention.
//
// Arguments:
//   operator_type (LadderOperatorType): radial operator type
//   space (basis::OscillatorOrbitalSpace): space
//   sectors (basis::OscillatorOrbitalSectors): sectors
//   matrices (basis::OperatorBlocks): operator matrices
{
  // validate operator quantum numbers
  const int L0 = sectors.L0();
  const int g0 = sectors.g0();
  assert((L0==1) && (g0==1));

  // initialize output matrices
  basis::SetOperatorToZero(sectors, matrices);

  // loop over sectors
  for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    // make aliases for sector and block
    const basis::OscillatorOrbitalSectors::SectorType& sector =
        sectors.GetSector(sector_index);
    basis::OperatorBlock<double>& sector_matrix = matrices[sector_index];

    // get sizes
    const basis::OscillatorOrbitalSubspace& bra_subspace = sector.bra_subspace();
    const basis::OscillatorOrbitalSubspace& ket_subspace = sector.ket_subspace();

    // angular momentum factor
    const int bra_l = bra_subspace.l();
    const int ket_l = ket_subspace.l();
    const int delta_l = bra_l - ket_l;
    double angular_matrix_element=0;

    if (delta_l==+1)
      angular_matrix_element = std::sqrt(float(ket_l+1)/(2*ket_l+3));
    else if (delta_l==-1)
      angular_matrix_element = std::sqrt(float(ket_l+1)/(2*ket_l+3));
    else
      assert(0);

    // main loop
    //
    // Note: Here we are looping over all matrix elements.  However, the ladder
    // operators have only a single nonzero diagonal, so at most a single ket is
    // connected to each bra within the matrix.  It would potentially be more
    // efficient to efficiently simply iterate over this diagonal.  However, for
    // this example, we do naive double iteration over the full matrix, to avoid
    // having to look up the bra state index from the bra qn connected to the
    // given ket.
    // #pragma omp parallel for collapse(2)
    for (std::size_t bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
    {
      for (std::size_t ket_index = 0; ket_index < ket_subspace.size(); ++ket_index)
      {
        // get states
        const basis::OscillatorOrbitalState bra_state(sector.bra_subspace(), bra_index);
        const basis::OscillatorOrbitalState ket_state(sector.ket_subspace(), ket_index);

        // calculate radial matrix element
        //
        // TODO mac (03/24/23): confirm radial wave function convention (this is
        // sigma=-1) and c vs. c-tilde convention
        const int bra_N = bra_state.N();
        const int ket_N = ket_state.N();
        const int delta_N = bra_N - ket_N;
        double radial_matrix_element = 0;
        const int N = bra_N, l = bra_l;  // aliases for use in matrix element formula
        if (operator_type == LadderOperatorType::kRaising)
          {
            if (delta_N == +1)
              {
                if (delta_l == +1)
                  radial_matrix_element = std::sqrt(N + l + 3);
                else if (delta_l == -1)
                  radial_matrix_element = -std::sqrt(N - l + 2);
              }
          }
        else if (operator_type == LadderOperatorType::kLowering)
          {
            if (delta_N == -1)
              {
                if (delta_l == +1)
                  radial_matrix_element = -std::sqrt(N - l);
                else if (delta_l == -1)
                  radial_matrix_element = std::sqrt(N + l + 1);
              }
          }

        // save full matrix element
        sector_matrix(bra_index, ket_index) = radial_matrix_element * angular_matrix_element;
      }
    }
  }
}

void TestLadderOperators()
{

  std::cout << "Ladder operators" << std::endl;
  std::cout << std::endl;

  // construct space
  int Nmax = 4;
  basis::OscillatorOrbitalSpace space(Nmax);

  // construct raising operator -- L0=1, g0=1
  basis::OscillatorOrbitalSectors raising_operator_sectors(space, 1, 1);
  basis::OperatorBlocks<double> raising_operator_matrices;
  LadderOperator(LadderOperatorType::kRaising, space, raising_operator_sectors, raising_operator_matrices);

  // inspect raising operator
  std::cout << "raising" << std::endl;
  for (std::size_t sector_index = 0; sector_index < raising_operator_sectors.size(); ++sector_index)
    {
      // make aliases for sector and block
      const basis::OscillatorOrbitalSectors::SectorType& sector =
        raising_operator_sectors.GetSector(sector_index);

      // print sector info
      std::cout << "  sector " << sector_index
                << "  bra index " << sector.bra_subspace_index()
                << " labels " << sector.bra_subspace().LabelStr()
                << " dim " << sector.bra_subspace().size()
                << "  ket index " << sector.ket_subspace_index()
                << " labels " << sector.ket_subspace().LabelStr()
                << " dim " << sector.ket_subspace().size()
                << "  multiplicity index " << sector.multiplicity_index()
                << std::endl;

      // print matrix
      std::cout << raising_operator_matrices[sector_index] << std::endl;
    }

  // construct lowering operator -- L0=1, g0=1
  basis::OscillatorOrbitalSectors lowering_operator_sectors(space, 1, 1);
  basis::OperatorBlocks<double> lowering_operator_matrices;
  LadderOperator(LadderOperatorType::kLowering, space, lowering_operator_sectors, lowering_operator_matrices);

  // inspect lowering operator
  std::cout << "lowering" << std::endl;
  for (std::size_t sector_index = 0; sector_index < lowering_operator_sectors.size(); ++sector_index)
    {
      // make aliases for sector and block
      const basis::OscillatorOrbitalSectors::SectorType& sector =
        lowering_operator_sectors.GetSector(sector_index);

      // print sector info
      std::cout << "  sector " << sector_index
                << "  bra index " << sector.bra_subspace_index()
                << " labels " << sector.bra_subspace().LabelStr()
                << " dim " << sector.bra_subspace().size()
                << "  ket index " << sector.ket_subspace_index()
                << " labels " << sector.ket_subspace().LabelStr()
                << " dim " << sector.ket_subspace().size()
                << "  multiplicity index " << sector.multiplicity_index()
                << std::endl;

      // print matrix
      std::cout << lowering_operator_matrices[sector_index] << std::endl;
    }
}

////////////////////////////////////////////////////////////////
// working with operators: constructing a product operator
////////////////////////////////////////////////////////////////

// construct number operator as spherical scalar product of ladder operators
//
// Apply Racah's "one-system" reduction formula, for the reduced matrix
// elements of the spherical tensor coupled product of two spherical tensor
// operators operating on a single space (here recall we are using Rose's
// normalization convention for the reduced matrix elements):
//
// < n' l || A_L0 . B_L0 || n l >
//     = sum_(n'' l'') (-)^(l''+l) * hat(l'') / hat(l)
//         * < n' l || A || n'' l''> < n'' l'' || B || n l >
//
// where A_L0 . B_L0 = (-)^L0 * hat(L0) * (A_L0 x B_L0)_00.
// [e.g., Brink & Satchler]
//
// Recall that in Rose convention, for a scalar operator, the RME and the
// matrix element are identical, with no angular momentum factor required.
//
// For our present case, L0=1.
//
// The sum over intermediate n'' is accomplished through the matrix
// multiplication of blocks.  The sum over intermediate l'' is subject to the
// triangle inequality triangle(l,L0,l) as well as the parity constraint
// pi(l'',1,l).  Thus, we have l'' in {l-1,l+1}, restricted to nonnegative
// values, while l''=l is excluded by parity.
//
// Thus, for each "target" sector, we must accumulate the contributions from
// ladder operator sectors corresponding to intermediate angular momenta l''.
// There are two approaches to accumulating these contributions:
//
// Approach #1: Blind loop over factor operator sectors
//
//   Iterate over all sectors of operators A and B, and see which pairs
//   satisfy the selection rules for the product.
//
// Approach #2: Loop over quantum numbers of the intermediate subspace
//
//   Here we identify the quantum numbers for the intermediate subspace,
//   based on the selection rules, and then we need to look up the
//   corresponding sector.
//
// However, either way, we must beware that, by hermiticity, only the
// "canonical" sectors in (bra_subspace_index,ket_subspace_index) actually
// appear in our indexing.  So, in approach #1, we must consider the lower
// triangle sectors as well, or, in approach #2, we must canonicalize the
// order of quantum numbers and look up the transpose if appropriate.
//
// Implementing these two approaches illustrates different functionality of
// the basis package.  Here are the attempts we make:
//
//   Approach1A -- Naively looping over sectors in the source (raising and
//   lowering) operators fails, since only the upper triangle sectors show up in
//   that iteration.  But we need some sectors from both the lower and upper
//   triangles of both operators.  We keep this attempt as a pedagogical example
//   of what goes wrong!
//
//   Approach1B -- We could patch this up by adding a second level of iteration,
//   so that for each sector of the source operators, we also consider the
//   (unstored) transpose sectors.  This seems messy.
//
//   Approach2A -- Naively looping over the allowed quantum numbers for the
//   intermediate subspace and retrieving the appropriate sectors of the source
//   operators again runs into trouble that only the upper triangle sectors are
//   defined and stored.  We again keep this attempt as a pedagogical example of
//   what goes wrong!
//
//   Approach2B -- As we loop over the quantum numbers for the intermediate
//   subspace, and retrieve the appropriate sectors of the source operators, we
//   make sure to "canonicalize" the sectors, flipping bra and ket subspaces as
//   necessary, to ensure that we only attempt to retrieve upper triangle
//   sectors.
//
//   Approach3A -- Naively looping over *all* subspaces as candidate
//   intermediate subspaces, and checking to see if the corresponding sectors
//   are defined for the source operators.  (This is what is done in shell's
//   obme.cpp, where the "full" operator, not just the upper triangle, is
//   stored.)  This makes sense if if is okay to assume that allowed sectors are
//   "dense", or, in any case, looping over all subspaces as candidate
//   intermediate subspaces is not overwhelmingly inefficient.
//
//   Approach3B -- As we check if the relevant sectors exist, first canonicalize
//   the subspace pair...

void PopulateProductOperatorApproach1A(
    const basis::OscillatorOrbitalSectors& raising_operator_sectors, const basis::OperatorBlocks<double>& raising_operator_matrices,
    const basis::OscillatorOrbitalSectors& lowering_operator_sectors, const basis::OperatorBlocks<double>& lowering_operator_matrices,
    const basis::OscillatorOrbitalSectors& number_operator_sectors, basis::OperatorBlocks<double>& number_operator_matrices
  )
// Naive implementation of Approach #1.  Fails to recognize that, for both the
// raising an lower operators, only the sectors constituting the upper triangle
// show up when iterating over sectors, and only this upper triangle is
// explicitly stored.
//
// Iteration only covers upper triangle sectors for both raising and lowering
// operators (and none of these even contribute!):
//
// Here is some of the diagnostic output (with comments):
//
//     target sector: [ 0 ] [ 0 ]
//
//          Here we are generating the l=0 target sector.
//
//       raising operator sector: [ 0 ] [ 1 ]
//       lowering operator sector: [ 0 ] [ 1 ]
//
//          Problem: This raising operator sector <l=0|l=1> (upper triangle)
//          needs to be paired with the lowering operator sector <l=1|l=0>
//          (lower triangle).
//
//       raising operator sector: [ 0 ] [ 1 ]
//       lowering operator sector: [ 1 ] [ 2 ]
//
//           Okay: This pair of sectors is genuinely irrelevant, since we are
//           generating a scalar operator.
{

  // for each number (target) operator block, accumulate contributions
  for (
      std::size_t number_operator_sector_index = 0;
      number_operator_sector_index < number_operator_sectors.size();
      ++number_operator_sector_index
    )
    {
      // make aliases for sector and block
      const basis::OscillatorOrbitalSectors::SectorType& number_operator_sector
        = number_operator_sectors.GetSector(number_operator_sector_index);
      basis::OperatorBlock<double>& number_operator_sector_matrix
        = number_operator_matrices[number_operator_sector_index];

      // extract target sector quantum numbers
      //
      // Note that this code is specialized to a rank 0 (dot) product, where
      // therefore l_bra=l_ket.
      assert(number_operator_sector.bra_subspace().l()==number_operator_sector.bra_subspace().l());
      int l = number_operator_sector.bra_subspace().l();
      std::cout << "target sector: " << number_operator_sector.bra_subspace().LabelStr()
                << " " << number_operator_sector.bra_subspace().LabelStr() << std::endl;

      // accumulate contributions

      // for each raising sector
      for (
          std::size_t raising_operator_sector_index = 0;
          raising_operator_sector_index < raising_operator_sectors.size();
          ++raising_operator_sector_index
        )
        {
          // for each lowering sector
          for (
              std::size_t lowering_operator_sector_index = 0;
              lowering_operator_sector_index < lowering_operator_sectors.size();
              ++lowering_operator_sector_index
            )

            {
              // make aliases for sectors and blocks
              const basis::OscillatorOrbitalSectors::SectorType& raising_operator_sector
                = raising_operator_sectors.GetSector(raising_operator_sector_index);
              const basis::OperatorBlock<double>& raising_operator_sector_matrix
                = raising_operator_matrices[raising_operator_sector_index];
              const basis::OscillatorOrbitalSectors::SectorType& lowering_operator_sector
                = lowering_operator_sectors.GetSector(lowering_operator_sector_index);
              const basis::OperatorBlock<double>& lowering_operator_sector_matrix
                = lowering_operator_matrices[lowering_operator_sector_index];

              // diagnostic output
              std::cout
                << " raising operator sector: " << raising_operator_sector.bra_subspace().LabelStr() << " " << raising_operator_sector.ket_subspace().LabelStr()
                << std::endl
                << " lowering operator sector: " << lowering_operator_sector.bra_subspace().LabelStr() << " " << lowering_operator_sector.ket_subspace().LabelStr()
                << std::endl;

              // check if this pair of raising and lowering operator sectors contributes
              //
              //   - "outer" bra and ket match those of target sector (l',l),
              //     which in this example, for a scalar product operator, is
              //     simply (l,l)
              //
              //   - "inner" bra and ket match those of intermediate subspace (l'')
              //
              // A comparison of labels() is a tuple comparison on (l,), which,
              // for this basis scheme, is equivalent to a simple integer
              // comparison on l().  We use the generic notation labels() to
              // provide an easily generalizable idiom.
              bool target_subspaces_match = (
                  (number_operator_sector.bra_subspace().labels() == raising_operator_sector.bra_subspace().labels())
                  &&
                  (lowering_operator_sector.ket_subspace().labels() == number_operator_sector.ket_subspace().labels())
                );
              bool intermediate_subspaces_match = (
                  raising_operator_sector.ket_subspace().labels() == lowering_operator_sector.bra_subspace().labels()
                );
              if (!(target_subspaces_match && intermediate_subspaces_match))
                continue;

              // calculate prefactor
              //
              //   (-)^(l''+l) * hat(l'') / hat(l)
              int lpp = raising_operator_sector.ket_subspace().l();
              double prefactor = ParitySign(l+lpp) * Hat(lpp) / Hat(l);
              number_operator_sector_matrix += prefactor * raising_operator_sector_matrix * lowering_operator_sector_matrix;

              // Note: In naive Approach 1A, we never even get here, since we
              // never encountered a valid pair of raising and lowering operator
              // sectors!

              std::cout << "matrix element calculation "
                        << " l " << l
                        << " lpp " << lpp
                        << " prefactor " << prefactor
                        << std::endl
                        << "raising matrix " << std::endl
                        << raising_operator_sector_matrix
                        << std::endl
                        << "lowering matrix " << std::endl
                        << lowering_operator_sector_matrix
                        << std::endl
                        << "product matrix " << std::endl
                        << raising_operator_sector_matrix * lowering_operator_sector_matrix
                        << std::endl;
            }
        }
    }

}

void PopulateProductOperatorApproach3A(
    const basis::OscillatorOrbitalSpace& space,
    const basis::OscillatorOrbitalSectors& raising_operator_sectors, const basis::OperatorBlocks<double>& raising_operator_matrices,
    const basis::OscillatorOrbitalSectors& lowering_operator_sectors, const basis::OperatorBlocks<double>& lowering_operator_matrices,
    const basis::OscillatorOrbitalSectors& number_operator_sectors, basis::OperatorBlocks<double>& number_operator_matrices
  )
// Naive implementation of Approach #3.  Fails to recognize that, for both the
// raising an lower operators, only the sectors constituting the upper triangle
// show up when iterating over sectors, and only this upper triangle is
// explicitly stored.
{

  // for each number (target) operator block, accumulate contributions
  for (
      std::size_t number_operator_sector_index = 0;
      number_operator_sector_index < number_operator_sectors.size();
      ++number_operator_sector_index
    )
    {
      // make aliases for sector and block
      const basis::OscillatorOrbitalSectors::SectorType& number_operator_sector
        = number_operator_sectors.GetSector(number_operator_sector_index);
      basis::OperatorBlock<double>& number_operator_sector_matrix
        = number_operator_matrices[number_operator_sector_index];

      // extract target sector quantum numbers
      //
      // Note that this code is specialized to a rank 0 (dot) product, where
      // therefore l_bra=l_ket.
      assert(number_operator_sector.bra_subspace().l()==number_operator_sector.bra_subspace().l());
      int l = number_operator_sector.bra_subspace().l();
      std::cout << "target sector: " << number_operator_sector.bra_subspace().LabelStr()
                << " " << number_operator_sector.bra_subspace().LabelStr() << std::endl;

      // extract target sector subspace indices
      const std::size_t number_operator_bra_subspace_index = number_operator_sector.bra_subspace_index();
      const std::size_t number_operator_ket_subspace_index = number_operator_sector.ket_subspace_index();

      for (
          std::size_t inner_subspace_index = 0; inner_subspace_index < space.size();
           ++inner_subspace_index
        )
        {

          // short circuit: check if resulting sectors exist in source operators
          if (
              !(
                  raising_operator_sectors.ContainsSector(number_operator_bra_subspace_index, inner_subspace_index)
                  && lowering_operator_sectors.ContainsSector(inner_subspace_index, number_operator_ket_subspace_index)
                )
            )
            continue;

          // retrieve corresponding blocks
          const std::size_t raising_operator_sector_index
            = raising_operator_sectors.LookUpSectorIndex(number_operator_bra_subspace_index, inner_subspace_index);
          const basis::OperatorBlock<double>& raising_operator_sector_matrix
            = raising_operator_matrices[raising_operator_sector_index];
          const std::size_t lowering_operator_sector_index
            = lowering_operator_sectors.LookUpSectorIndex(inner_subspace_index, number_operator_ket_subspace_index);
          const basis::OperatorBlock<double>& lowering_operator_sector_matrix
            = lowering_operator_matrices[lowering_operator_sector_index];

//               // diagnostic output
//               std::cout
//                 << " raising operator sector: " << raising_operator_sector.bra_subspace().LabelStr() << " " << raising_operator_sector.ket_subspace().LabelStr()
//                 << std::endl
//                 << " lowering operator sector: " << lowering_operator_sector.bra_subspace().LabelStr() << " " << lowering_operator_sector.ket_subspace().LabelStr()
//                 << std::endl;
//
//               // check if this pair of raising and lowering operator sectors contributes
//               //
//               //   - "outer" bra and ket match those of target sector (l',l),
//               //     which in this example, for a scalar product operator, is
//               //     simply (l,l)
//               //
//               //   - "inner" bra and ket match those of intermediate subspace (l'')
//               //
//               // A comparison of labels() is a tuple comparison on (l,), which,
//               // for this basis scheme, is equivalent to a simple integer
//               // comparison on l().  We use the generic notation labels() to
//               // provide an easily generalizable idiom.
//               bool target_subspaces_match = (
//                   (number_operator_sector.bra_subspace().labels() == raising_operator_sector.bra_subspace().labels())
//                   &&
//                   (lowering_operator_sector.ket_subspace().labels() == number_operator_sector.ket_subspace().labels())
//                 );
//               bool intermediate_subspaces_match = (
//                   raising_operator_sector.ket_subspace().labels() == lowering_operator_sector.bra_subspace().labels()
//                 );
//               if (!(target_subspaces_match && intermediate_subspaces_match))
//                 continue;
//
//               // calculate prefactor
//               //
//               //   (-)^(l''+l) * hat(l'') / hat(l)
//               int lpp = raising_operator_sector.ket_subspace().l();
//               double prefactor = ParitySign(l+lpp) * Hat(lpp) / Hat(l);
//               number_operator_sector_matrix += prefactor * raising_operator_sector_matrix * lowering_operator_sector_matrix;
//
//               // Note: In naive Approach 1A, we never even get here, since we
//               // never encountered a valid pair of raising and lowering operator
//               // sectors!
//
//               std::cout << "matrix element calculation "
//                         << " l " << l
//                         << " lpp " << lpp
//                         << " prefactor " << prefactor
//                         << std::endl
//                         << "raising matrix " << std::endl
//                         << raising_operator_sector_matrix
//                         << std::endl
//                         << "lowering matrix " << std::endl
//                         << lowering_operator_sector_matrix
//                         << std::endl
//                         << "product matrix " << std::endl
//                         << raising_operator_sector_matrix * lowering_operator_sector_matrix
//                         << std::endl;
        }
    }
}

void TestProductOperator()
{

  std::cout << "Number operator as product operator" << std::endl;
  std::cout << std::endl;

  // construct space
  int Nmax = 4;
  basis::OscillatorOrbitalSpace space(Nmax);

  // construct raising operator -- L0=1, g0=1
  basis::OscillatorOrbitalSectors raising_operator_sectors(space, 1, 1);
  basis::OperatorBlocks<double> raising_operator_matrices;
  LadderOperator(LadderOperatorType::kRaising, space, raising_operator_sectors, raising_operator_matrices);

  // construct lowering operator -- L0=1, g0=1
  basis::OscillatorOrbitalSectors lowering_operator_sectors(space, 1, 1);
  basis::OperatorBlocks<double> lowering_operator_matrices;
  LadderOperator(LadderOperatorType::kLowering, space, lowering_operator_sectors, lowering_operator_matrices);

  // construct (empty) product operator -- L0=0, g0=0
  basis::OscillatorOrbitalSectors number_operator_sectors(space, 0, 0);
  basis::OperatorBlocks<double> number_operator_matrices;

  // construct zeroed-out number operator
  basis::SetOperatorToZero(number_operator_sectors, number_operator_matrices);

  PopulateProductOperatorApproach1A(
      raising_operator_sectors, raising_operator_matrices,
      lowering_operator_sectors, lowering_operator_matrices,
      number_operator_sectors, number_operator_matrices
    );

  // inspect resulting product (number?) operator
  std::cout << "number" << std::endl;
  for (std::size_t sector_index = 0; sector_index < number_operator_sectors.size(); ++sector_index)
    {
      // make aliases for sector and block
      const basis::OscillatorOrbitalSectors::SectorType& sector =
        number_operator_sectors.GetSector(sector_index);

      // print sector info
      std::cout << "  sector " << sector_index
                << "  bra index " << sector.bra_subspace_index()
                << " labels " << sector.bra_subspace().LabelStr()
                << " dim " << sector.bra_subspace().size()
                << "  ket index " << sector.ket_subspace_index()
                << " labels " << sector.ket_subspace().LabelStr()
                << " dim " << sector.ket_subspace().size()
                << "  multiplicity index " << sector.multiplicity_index()
                << std::endl;

      // print matrix
      std::cout << number_operator_matrices[sector_index] << std::endl;
    }


}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  std::cout << "Oscillator orbital basis" << std::endl;
  std::cout << std::endl;

  // TestOscillatorOrbitalSubspace();
  // TestOscillatorOrbitalSpace();
  TestOscillatorOrbitalSectors();
  // TestLadderOperators();
  // TestProductOperator();

  // termination
  return 0;
}
