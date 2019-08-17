/****************************************************************
  @file proton_neutron.h

  Provides convenience definitions for proton-neutron bookkeeping.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 07/01/17 (mac): Extracted from jjjpn_scheme.h.
  + 10/19/17 (mac): Define OperatorTypePN.
  + 02/27/19 (pjf): Add reverse lookups for OrbitalSpeciesPN.
  + 08/14/19 (pjf): Add reverse lookup for OperatorTypePN.

****************************************************************/

#ifndef BASIS_PROTON_NEUTRON_H_
#define BASIS_PROTON_NEUTRON_H_

#include <array>
#include <string>
#include <unordered_map>

#include "am/halfint.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle orbital label types and coding schemes
  ////////////////////////////////////////////////////////////////

  // enumerated type for orbital species
  //
  // Note: Follows same sequence as MFDn, but MFDn uses 1-based
  // numbering.

  enum class OrbitalSpeciesPN : int {kP=0,kN=1};

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

  // notational reverse definitions for orbital species
  //
  // Example:
  //   std::string orbital_species_code = "p";
  //   ...
  //   os << basis::kTzOrbitalSpeciesPNCode[orbital_species_code];
  //
  // Note: Tz uses "up quark is positive" convention.
  extern const std::unordered_map<HalfInt, OrbitalSpeciesPN> kTzCodeOrbitalSpeciesPN;
  extern const std::unordered_map<int, OrbitalSpeciesPN> kIndexCodeOrbitalSpeciesPN;
  extern const std::unordered_map<int, OrbitalSpeciesPN> kDecimalCodeOrbitalSpeciesPN;
  extern const std::unordered_map<std::string, OrbitalSpeciesPN> kCharCodeOrbitalSpeciesPN;

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
  // proton-neutron operator types
  ////////////////////////////////////////////////////////////////

  // enumerated type for proton-neutron operator types

  enum class OperatorTypePN {kP=0,kN=1,kTotal=2};

  // notational definitions for proton-neutron operator types
  //
  // Use of these arrays requires conversion of the OperatorTypePN to int.
  extern const std::array<const char*,3> kOperatorTypePNCodeChar;  // {"p","n","total"}

  // notational reverse definitions for proton-neutron operator types
  //
  // Example:
  //   std::string operator_species_code = "p";
  //   ...
  //   os << kCharCodeOperatorTypePN[operator_species_code];
  extern const std::unordered_map<std::string,OperatorTypePN> kCharCodeOperatorTypePN;


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
