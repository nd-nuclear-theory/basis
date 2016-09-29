/****************************************************************
  nlj_orbital.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>

#include "am/am.h"

#include "nlj_orbital.h"

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
  ////////////////////////////////////////////////////////////////
} // namespace
