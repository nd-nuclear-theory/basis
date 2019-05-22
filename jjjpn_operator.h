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

****************************************************************/

#ifndef BASIS_JJJPN_OPERATOR_H_
#define BASIS_JJJPN_OPERATOR_H_

#include "basis/jjjpn_scheme.h"
#include "basis/operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body jjJpn operator output
  ////////////////////////////////////////////////////////////////

  // Note that the primary intention of this output for two-body
  // operators in JJJPN scheme is for diagnostic purposes.  It is
  // based on the preliminary specification for the text version of
  // the MFDn Version 15 h2 format.
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
