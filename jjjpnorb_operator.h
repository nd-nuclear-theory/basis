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

  // Note that the primary intention of the output for two-body
  // operators in JJJT scheme is for diagnostic purposes.  
  //
  // Data lines are of the form:
  //
  //   T0  N1' l1' j1' N2' l2' J' T' g'  N1 l1 j1 N2 l2 j2 J T g  JT-RME
  //
  // Although the g label is redundant (it can be deduced from l1 and
  // l2), it is included to make the sector structure more easily
  // apparent to a human reader.
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
  // Write single isospin component of a two-body operator in JJJT
  // scheme.
  //
  // Side effect: The floating point precision attribute of the output
  // stream is modified.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   T0 (int) : isospin for this isospin component
  //   sector (basis::TwoBodySectorsJJJT) : sectors defining operator
  //   matrices (basis::MatrixVector) : matrices defining operator
  //   conversion (basis::NormalizationConversion) : specifies any 
  //     conversion between AS and NAS for output




  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
