/****************************************************************
  jjjpn_scheme_general.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>

#include "am/am.h"

#include "jjjpnorb_scheme.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle orbitals
  ////////////////////////////////////////////////////////////////

  OrbitalSubspacePN::OrbitalSubspacePN(OrbitalSpeciesPN orbital_species, int Nmax)
  {

    // set values
    labels_ = SubspaceLabelsType(orbital_species);
    weight_max_ = double(Nmax);
    Nmax_ = Nmax;

    // iterate over total oscillator quanta
    for (int N = 0; N <= Nmax; ++N)
      // iterate over j within shell
      for (HalfInt j = HalfInt(1,2); j <= N+HalfInt(1,2); ++j) 
        {
          // recover derived quantum numbers (n,l) from (N,j)
          int l = (TwiceValue(j)-1)/2 + (N+(TwiceValue(j)-1)/2)%2;
          int n = (N-l)/2;
            
          // save state
          PushStateLabels(StateLabelsType(n,l,j)); 
            
          // save oscillator quantum number as weight
          weights_.push_back(double(N));
        }
  }

  std::string OrbitalSubspacePN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << int(orbital_species())
       << " " << "]";

    return os.str();
  }

  std::string OrbitalSubspacePN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        OrbitalStatePN state(*this,state_index);

        os
	  << " " << "index"
	  << " " << std::setw(width) << state_index
	  << " " << "nlj"
	  << " " << std::setw(width) << state.n()
	  << " " << std::setw(width) << state.l()
	  << " " << std::setw(width+2) << state.j().Str()
	  << " " << "weight"
          << " " << state.weight()
          << std::endl;
      }

    return os.str();

  }


  OrbitalSpacePN::OrbitalSpacePN(int Nmax)
  {

    // save truncation
    weight_max_ = double(Nmax);
    Nmax_ = Nmax;

    // iterate over species
    for (OrbitalSpeciesPN orbital_species : {OrbitalSpeciesPN::kP,OrbitalSpeciesPN::kN})
      {
        OrbitalSubspacePN subspace(orbital_species,Nmax);
        PushSubspace(subspace);
      }
  }

  std::string OrbitalSpacePN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);
        
	os
	  << " " << "index"
	  << " " << std::setw(width) << subspace_index
	  << " " << "species"
	  << " " << std::setw(width) << int(subspace.orbital_species()) 
	  << " " << "dim"
	  << " " << std::setw(width) << subspace.size()
	  << " " << std::endl;
      }

    return os.str();

  }


  std::string OrbitalSpacePN::OrbitalDefinitionStr() const
  {

    std::ostringstream os;

    // header comments
    os << "# MFDn SPorbital file" << std::endl;
    os << "#   version" << std::endl;
    os << "#   norb_p norb_n" << std::endl;
    os << "#   index n l 2*j species weight" << std::endl;

    // header line 1: version
    int version = 15055;
    os << version << std::endl;

    // header line 2: dimensions
    os << GetSubspace(0).size() << " " << GetSubspace(1).size() << std::endl;
        
    // data
    const int width = 3;
    const int precision = 8;
    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        // retrieve subspace
	const SubspaceType& subspace = GetSubspace(subspace_index);

        // iterate over states
        for (int state_index=0; state_index<subspace.size(); ++state_index)
          {
            OrbitalStatePN state(subspace,state_index);
            
            os
              << " " << std::setw(width) << state.index()+1  // 1-based
              << " " << std::setw(width) << state.n()
              << " " << std::setw(width) << state.l()
              << " " << std::setw(width) << TwiceValue(state.j())
              << " " << std::setw(width) << int(state.orbital_species())+1 // 1-based
              << " " << std::fixed << std::setw(width+1+precision)
              << std::setprecision(precision) << state.weight()
              << std::endl;
          }
      }
    
    return os.str();

  }

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJpn scheme with general orbitals
  ////////////////////////////////////////////////////////////////

  // notational definitions for two-body state species -- static array definitions
  const std::array<int,3> kTwoBodySpeciesPNCodeTz({+1,-1,0});
  const std::array<int,3> kTwoBodySpeciesPNCodeDecimal({11,22,12});
  const std::array<const char*,3> kTwoBodySpeciesPNCodeChar({"pp","nn","pn"});


  TwoBodySubspaceJJJPN::TwoBodySubspaceJJJPN(
      const OrbitalSpacePN& orbital_space,
      TwoBodySpeciesPN two_body_species, int J, int g,
      const WeightMax& weight_max
    )
  {

    // std::cout << "Subspace construction"
    //           << " " << int(two_body_species)
    //           << " " << J
    //           << " " << g
    //           << std::endl;

    // set values
    labels_ = SubspaceLabelsType(two_body_species,J,g);
    weight_max_ = weight_max;

    // identify orbital subspaces
    OrbitalSpeciesPN orbital_species1 = 
      (two_body_species==TwoBodySpeciesPN::kPP) || (two_body_species==TwoBodySpeciesPN::kPN)
      ? OrbitalSpeciesPN::kP : OrbitalSpeciesPN::kN;
    OrbitalSpeciesPN orbital_species2 = 
      (two_body_species==TwoBodySpeciesPN::kPP) 
      ? OrbitalSpeciesPN::kP : OrbitalSpeciesPN::kN;
    orbital_subspace1_ptr_ = &(orbital_space.GetSubspace(int(orbital_species1)));
    orbital_subspace2_ptr_ = &(orbital_space.GetSubspace(int(orbital_species2)));

    // set up indexing
    // iterate over orbitals for particle 1
    for (int index1 = 0; index1 < orbital_subspace1().size(); ++index1)
      // iterate over orbitals for particle 2
      for (int index2 = 0; index2 < orbital_subspace2().size(); ++index2)
        {

          // retrieve orbitals
          OrbitalStatePN orbital1(orbital_subspace1(),index1);
          OrbitalStatePN orbital2(orbital_subspace2(),index2);

          // for like-particle subspaces
          if ((two_body_species==TwoBodySpeciesPN::kPP) || (two_body_species==TwoBodySpeciesPN::kNN))
            {
              // impose canonical ordering on single-particle states
              if (!(index1<=index2))
                continue;

              // impose antisymmetry
              if ((index1==index2)&&(!(J%2==0)))
                continue;
            }

          // impose triangularity
          if (!(am::AllowedTriangle(orbital1.j(),orbital2.j(),J)))
            continue;

          // impose parity
          if (!((orbital1.g()+orbital2.g()+g)%2==0))
            continue;

          // impose weight cutoffs
          double w1 = orbital1.weight();
          double w2 = orbital2.weight();
          double w = w1 + w2;
          double w1_max = weight_max.one_body[int(orbital_species1)];
          double w2_max = weight_max.one_body[int(orbital_species2)];
          double w_max = weight_max.two_body[int(two_body_species)];
          bool weight_allowed = ((w1<=w1_max)&&(w2<=w2_max)&&(w<=w_max));
          // std::cout << " " << w1
          //           << " " << w2
          //           << " " << w
          //           << " " << w1_max
          //           << " " << w2_max
          //           << " " << w_max
          //           << " " << weight_allowed 
          //           << std::endl;
          if (!weight_allowed)
            continue;

          // keep surviving states
          PushStateLabels(StateLabelsType(index1,index2)); 
        }
  }
  
  std::string TwoBodySubspaceJJJPN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << kTwoBodySpeciesPNCodeTz[int(two_body_species())]
       << " " << std::setw(width) << J() 
       << " " << std::setw(width) << g()
       << " " << "]";

    return os.str();
  }

  std::string TwoBodySubspaceJJJPN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {

        // retrieve two-body state
        TwoBodyStateJJJPN state(*this,state_index);

        // retrieve orbitals
        OrbitalStatePN orbital1(state.orbital_subspace1(),state.index1());
        OrbitalStatePN orbital2(state.orbital_subspace2(),state.index2());

        os
	  << " " << "index"
	  << " " << std::setw(width) << state_index
	  << " " << "index1 index2"
	  << " " << std::setw(width) << state.index1()
	  << " " << std::setw(width) << state.index2()
	  << " " << "nlj1 nlj2"
	  << " " << std::setw(width) << orbital1.n()
	  << " " << std::setw(width) << orbital1.l()
	  << " " << std::setw(width+2) << orbital1.j().Str()
	  << " " << std::setw(width) << orbital2.n()
	  << " " << std::setw(width) << orbital2.l()
	  << " " << std::setw(width+2) << orbital2.j().Str()
          << std::endl;
      }

    return os.str();

  }

  std::string TwoBodyStateJJJPN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << kTwoBodySpeciesPNCodeTz[int(two_body_species())]
       << " " << std::setw(width) << J() 
       << " " << std::setw(width) << g()
       << " " << ";"
       << " " << std::setw(width) << index1()
       << " " << std::setw(width) << index2()
       << " " << "("
       << " " << std::setw(width) << GetOrbital1().N()
       << " " << std::setw(width) << GetOrbital1().j()
       << " " << std::setw(width) << GetOrbital2().N()
       << " " << std::setw(width) << GetOrbital2().j()
       << " " << ")"
       << " " << "]";

    return os.str();
  }


  TwoBodySpaceJJJPN::TwoBodySpaceJJJPN(
      const OrbitalSpacePN& orbital_space,
      const WeightMax& weight_max
    )
  {

    // save truncation
    weight_max_ = weight_max;

    // find putative Jmax from maximal j among orbitals
    HalfInt jmax=HalfInt(1,2);
    for (int subspace_index=0; subspace_index<orbital_space.size(); ++subspace_index)
      {
        // retrieve subspace
	const OrbitalSubspacePN& subspace = orbital_space.GetSubspace(subspace_index);

        // iterate over states
        for (int state_index=0; state_index<subspace.size(); ++state_index)
          {
            OrbitalStatePN state(subspace,state_index);
            jmax = std::max(jmax,state.j());
          }        
      }
    int Jmax = TwiceValue(jmax);

    // iterate over two-body species
    for (TwoBodySpeciesPN two_body_species : {TwoBodySpeciesPN::kPP,TwoBodySpeciesPN::kNN,TwoBodySpeciesPN::kPN})
      // iterate over J
      for (int J=0; J<=Jmax; ++J)
      {
        // iterate over g
        for (int g=0; g<=1; ++g)
          {
            TwoBodySubspaceJJJPN subspace(
                orbital_space,
                two_body_species,J,g,
                weight_max
              );

            if (subspace.size()!=0)
              PushSubspace(subspace);
          }
      }
  }


  std::string TwoBodySpaceJJJPN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);
	os
	  << " " << "index"
	  << " " << std::setw(width) << subspace_index
	  << " " << "sJg"
          << " " << std::setw(width) << int(subspace.two_body_species())
	  << " " << std::setw(width) << subspace.J() 
	  << " " << std::setw(width) << subspace.g() 
	  << " " << "dim"
	  << " " << std::setw(width) << subspace.size()
	  << " " << std::endl;
      }

    return os.str();

  }


  TwoBodySectorsJJJPN::TwoBodySectorsJJJPN(
      const TwoBodySpaceJJJPN& space,
      int J0, int g0,
      basis::SectorDirection sector_direction
    )
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
	{

          // enforce canonical ordering
          if (
              (sector_direction == basis::SectorDirection::kCanonical)
              && !(bra_subspace_index<=ket_subspace_index)
            )
            continue;

          // retrieve subspaces
          const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
          const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);

          // enforce particle conservation
          //
          // FUTURE: to change to fixed delta Tz condition
          if (!(bra_subspace.two_body_species()==ket_subspace.two_body_species()))
              continue;

          // verify angular momentum, isosopin, and parity selection rules
          bool allowed = true;
          allowed &= am::AllowedTriangle(ket_subspace.J(),J0,bra_subspace.J());
          allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);

          // push sector
	  if (allowed)
            PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
        }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
