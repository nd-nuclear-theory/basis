/****************************************************************
  @file proton_neutron.h

  Provides convenience definitions for proton-neutron bookkeeping.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 7/1/17 (mac): Extracted from jjjpn_scheme.h.

****************************************************************/

#ifndef BASIS_PROTON_NEUTRON_H_
#define BASIS_PROTON_NEUTRON_H_

#include <array>

#include "am/halfint.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle orbital label types and coding schemes
  ////////////////////////////////////////////////////////////////

  // enumerated type for orbital species
  //
  // Note: Follows same sequence as MFDn, but MFDn uses 1-based
  // numbering.

  enum class OrbitalSpeciesPN {kP=0,kN=1};

  // notational definitions for orbital species
  //
  // Use of this array requires conversion of the OrbitalSpeciesPN to int.
  //
  // Example:
  //   basis::OrbitalSpeciesPN orbital_species;
  //   ...
  //   os << basis::kOrbitalSpeciesPNCodeTz[int(orbital_species)];
  //
  // Note: Tz uses "up quark is positive" convention.
  extern const std::array<HalfInt, 2> kOrbitalSpeciesPNCodeTz;
  extern const std::array<int, 2> kOrbitalSpeciesPNCodeDecimal;
  extern const std::array<const char*, 2> kOrbitalSpeciesPNCodeChar;


  ////////////////////////////////////////////////////////////////
  // two-body label types and coding schemes
  ////////////////////////////////////////////////////////////////

  // enumerated type for two-body state species
  //
  // Note: Ordering of pp/nn/pn labels follows the same sequence as
  // used in MFDn.  However, note that MFDn uses 1-based numbering.

  enum class TwoBodySpeciesPN {kPP=0,kNN=1,kPN=2};

  // notational definitions for two-body state species
  //
  // Use of these arrays requires conversion of the TwoBodySpeciesPN to int.
  //
  // Example:
  //   basis::TwoBodySpeciesPN two_body_species;
  //   ...
  //   os << basis::kTwoBodySpeciesPNCodeTz[int(two_body_species)];
  extern const std::array<int,3> kTwoBodySpeciesPNCodeTz;  // {+1,-1,0} -- "up quark is positive" convention
  extern const std::array<int,3> kTwoBodySpeciesPNCodeDecimal;  // {11,22,12} -- used by MFDn
  extern const std::array<const char*,3> kTwoBodySpeciesPNCodeChar;  // {"pp","nn","pn"}

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
