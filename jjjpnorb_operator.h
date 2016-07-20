/****************************************************************
  jjjpnorb_operator.h

  Defines functions for I/O and manipulation of two-body operator
  matrices in jjJpn coupling scheme, based on general single-particle
  orbital sets.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/19/16 (mac): Created, adapting code from jjjt_operator.

****************************************************************/

#ifndef JJJPNORB_OPERATOR_H_
#define JJJPNORB_OPERATOR_H_

#include "basis/jjjpnorb_scheme.h"
#include "basis/operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body jjJpn operator output
  ////////////////////////////////////////////////////////////////

  // Note that the primary intention of this output for two-body
  // operators in JJJPN scheme is for diagnostic purposes.
  //
  // Other more specialized formats are used by MFDn.
  //
  // Data lines are of the form:
  //
  //   i1' i2' J' g' Tz'   i1 i2 J g Tz   JT-RME
  //
  // The single particle orbital indices (i1,i2,...) may be written
  // either 0-based (native) or 1-based (e.g., MFDn convention).
  //
  // The isospin projection Tz is defined under the convention that
  // protons (or up quarks) are positive.
  //
  // Although the g label is redundant (it can be deduced from the
  // orbitals), it is included to make the sector structure more
  // easily apparent to a human reader.
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
      const basis::MatrixVector& matrices,
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
  //   matrices (basis::MatrixVector) : matrices defining operator
  //   conversion (basis::NormalizationConversion) : specifies any 
  //     conversion between AS and NAS for output
  //   indexing_base (int) : use 0-based or 1-based indexing

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
