/************************************************************//**
  @file jjjpn_operator.h

  Defines functions for I/O and manipulation of two-body operator
  matrices in jjJpn coupling scheme, based on general single-particle
  orbital sets.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 7/19/16 (mac): Created, adapting code from jjjt_operator (jjjpnorb_operator).
  + 10/15/16 (mac): Update to use new sector method IsUpperTriangle().
  + 10/25/16 (mac): Rename to jjjpn_operator.
  + 9/27/17 (mac): Fix AS/NAS conversion factors in WriteTwoBodyOperatorJJJPN.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.
  + 03/07/23 (mac): Add canonicalization routines (refactored from xpn2h2).

****************************************************************/

#ifndef BASIS_JJJPN_OPERATOR_H_
#define BASIS_JJJPN_OPERATOR_H_

#include <iosfwd>

#include "jjjpn_scheme.h"
#include "many_body.h"
#include "operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body jjJpn state lookup (with canonicalization) 
  ////////////////////////////////////////////////////////////////

  std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,double> CanonicalizeTwoBodyOrbitalIndicesJJJPN(
      const OrbitalSpacePN& space,
      int J,
      std::size_t subspace_index1, std::size_t subspace_index2,
      std::size_t state_index1, std::size_t state_index2
    );
  // Convert subspace (proton-neutron) and state (orbital) indices for a JJJPN
  // two-particle state to canonical indices for state look-up.
  //
  // Canonicalization of orbitals is lexicographical by (subspace,state), so (1)
  // np state will be swapped to pn, then (2) like-species orbitals will be
  // swapped to numerical order.
  //
  // The canonicalization factor for swapping orbitals is -(-)^(J-j1-j2).  See,
  // e.g., equation (C3) of csbasis
  // [http://dx.doi.org/10.1103/PhysRevC.86.034312].
  //
  // Arguments:
  //
  //   space (basis::OrbitalSpacePN): space, for retrieving
  //     orbital quantum numbers to calculate canonicalization factor
  //   J (int): two-particle state J, to calculate canonicalization factor
  //   subspace_index1, subspace_index2 (std::size_t):
  //     orbital subspace indices, possibly to be swapped (if np state given)
  //     (may be implicitly cast as int from basis::OrbitalSpeciesPN)
  //   state_index1, state_index2 (std::size_t):
  //     naive orbital indices, possibly to be swapped if sector
  //     is diagonal sector
  //
  // Returns:
  //   canonicalized indices and swap flag and phase as:
  //
  //        subspace_index1, subspace_index2,
  //        state_index1, state_index2,
  //        canonicalization_factor

  ////////////////////////////////////////////////////////////////
  // two-body jjJpn RME lookup (with canonicalization) 
  ////////////////////////////////////////////////////////////////
  
  std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,double> CanonicalizeIndicesJJJPN(
      const TwoBodySpaceJJJPN& space,
      int J0, int g0,
      std::size_t subspace_index_bra, std::size_t subspace_index_ket,
      std::size_t state_index_bra, std::size_t state_index_ket
    );
  // Convert subspace and state indices for a matrix element to
  // canonical ("upper triangle") indices.
  //
  // It is assumed that the operator is a standard hermitian self-adjoint
  // spherical tensor operator.
  //
  // This is a customized wrapper for basis::CanonicalizeIndices (see
  // operator.h), for use with JJJPN operators.
  //
  // Arguments:
  //   space (basis::TwoBodySpaceJJJPN): space, for retrieving
  //     subspace quantum numbers to calculate canonicalization factor
  //   J0, g0 (int): operator tensorial properties
  //   bra_subspace_index, ket_subspace_index (std::size_t):
  //     naive sector bra and ket subspace indices, possibly to be swapped
  //   bra_state_index, ket_state_index (std::size_t):
  //     naive bra and ket state indices, possibly to be swapped if sector
  //     is diagonal sector
  //
  // Returns:
  //   canonicalized indices and swap flag and phase as:
  //
  //        subspace_index_bra, subspace_index_ket,
  //        state_index_bra, state_index_ket,
  //        swapped_subspaces,
  //        canonicalization_factor

  
  ////////////////////////////////////////////////////////////////
  // two-body jjJpn operator output
  ////////////////////////////////////////////////////////////////

  // Note that the primary intention of this output for two-body
  // operators in JJJPN scheme is for diagnostic purposes.  It is
  // based on the preliminary specification for the text version of
  // the MFDn Version 15 h2 format, i.e., h2v0.
  //
  // Data lines are of the form:
  //
  //   i1' i2' J' g' Tz'   i1 i2 J g Tz   JT-RME
  //
  // The single particle orbital indices (i1,i2,...) may be written
  // either 0-based (native) or 1-based (e.g., MFDn convention), as
  // controlled by the parameter indexing_base.
  //
  // The isospin projection Tz is defined under the convention that
  // protons (or up quarks) are positive, i.e., "pp" -> +1, "pn" -> 0,
  // and "nn" -> -1.
  //
  // Although the parity grade label g is redundant (it can be deduced
  // from the l values of the orbitals), it is included to make the
  // sector structure more easily apparent to a human reader.
  //
  // Iteration follows the usual scheme within the basis module:
  // sectors are lexicographic by (bra,ket) subspace indices, then
  // matrix elements within a sector are lexicographic by (bra,ket)
  // state indices.
  //
  // Reminder: One should be sure to document whether one is writing
  // AS or NAS matrix elements!

  void WriteTwoBodyOperatorJJJPN(
      std::ostream& os,
      const basis::TwoBodySectorsJJJPN& sectors,
      const basis::OperatorBlocks<double>& matrices,
      basis::NormalizationConversion conversion_mode,
      int indexing_base
    );
  // Write two-body operator in JJJPN scheme.
  //
  // Side effect: The floating point precision attribute of the output
  // stream is modified.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   sector (basis::TwoBodySectorsJJJPN) : sectors defining operator
  //   matrices (basis::OperatorBlocks<double>) : matrices defining operator
  //   conversion (basis::NormalizationConversion) : specifies any
  //     conversion between AS and NAS for output
  //   indexing_base (int) : use 0-based or 1-based indexing

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
