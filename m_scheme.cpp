/****************************************************************
  m_scheme.cpp

  Patrick J. Fasano
  University of Notre Dame
  European Centre for Theoretical Studies
    in Nuclear Physics and Related Areas

****************************************************************/


#include <cstddef>
#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>
#include <set>
#include <algorithm>

#include "am/am.h"
#include "mcutils/parsing.h"
#include "basis/proton_neutron.h"

#include "basis/m_scheme.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle states
  ////////////////////////////////////////////////////////////////

  /**
   * Construct a subspace with a particular species from an OrbitalSubspacePN.
   *
   * @param[in] orbital_subspace orbital subspace
   */
  SingleParticleSubspacePN::SingleParticleSubspacePN(const OrbitalSubspacePN& orbital_subspace)
  {

    // set values
    labels_ = SubspaceLabelsType(orbital_subspace.labels());
    weight_max_ = orbital_subspace.weight_max();

    // iterate over all states, enumerating m-substates
    for (std::size_t i=0; i < orbital_subspace.size(); ++i) {
      OrbitalStatePN orbital(orbital_subspace, i);
      for (HalfInt m = -orbital.j(); m <= orbital.j(); ++m) {
        PushStateLabels(StateLabelsType(orbital.n(),orbital.l(),orbital.j(),m));
        weights_.push_back(orbital.weight());
      }
    }

    // check if oscillator-like and set Nmax
    is_oscillator_like_ = orbital_subspace.is_oscillator_like();
    if (is_oscillator_like()) {
      Nmax_ = orbital_subspace.Nmax();
    } else {
      Nmax_ = -1;
    }


  }

  /**
   * Generate a string representation of the subspace labels.
   * @return subspace labels as a string
   */
  std::string SingleParticleSubspacePN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << int(orbital_species())
       << " " << "]";

    return os.str();
  }

  /**
   * Generate a string representation, useful for debugging.
   * @return debug string
   */
  std::string SingleParticleSubspacePN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    std::string oscillator_like_indicator = (is_oscillator_like() ? "true" : "false");
    os << " weight_max " << weight_max()
       << " Nmax " << Nmax()
       << " (oscillator-like: " << oscillator_like_indicator << ")"
       << std::endl;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        SingleParticleStatePN state(*this,state_index);

        os
          << " " << "index"
          << " " << std::setw(width) << state_index
          << " " << "nljm"
          << " " << std::setw(width) << state.n()
          << " " << std::setw(width) << state.l()
          << " " << std::setw(width+2) << state.j().Str()
          << " " << std::setw(width+3) << state.m().Str()
          << " " << "weight"
          << " " << state.weight()
          << std::endl;
      }

    return os.str();

  }

  /**
   * Generate a string representation of the orbital labels.
   * @return orbital labels as a string
   */
  std::string SingleParticleStatePN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << int(orbital_species())
       << " " << std::setw(width) << index()
       << " :"
       << " " << std::setw(width) << n()
       << " " << std::setw(width) << l()
       << " " << std::setw(width) << j()
       << " " << std::setw(width) << m()
       << " :"
       << " " << std::setw(width) << weight()
       << " " << "]";

    return os.str();
  }

  /**
   * Construct a space with species subspaces from an OrbitalSpacePN.
   *
   * @param[in] states vector of orbitals
   */
  SingleParticleSpacePN::SingleParticleSpacePN(const OrbitalSpacePN& orbital_space)
  {
    weight_max_ = orbital_space.weight_max();

    // construct subspaces
    for (std::size_t i=0; i<orbital_space.size(); ++i)
    {
      const OrbitalSubspacePN& orbital_subspace = orbital_space.GetSubspace(i);
      PushSubspace(SingleParticleSubspacePN(orbital_subspace));
    }

    // check if oscillator-like and set Nmax
    is_oscillator_like_ = orbital_space.is_oscillator_like();
    if (is_oscillator_like()) {
      Nmax_ = orbital_space.Nmax();
    } else {
      Nmax_ = -1;
    }

  }

  /**
   * Generate a string representation, useful for debugging.
   * @return debug string
   */
  std::string SingleParticleSpacePN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    std::string oscillator_like_indicator = (is_oscillator_like() ? "true" : "false");
    os << " weight_max " << weight_max()
       << " Nmax " << Nmax()
       << " (oscillator-like: " << oscillator_like_indicator << ")"
       << std::endl;

    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
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

}  // namespace basis
