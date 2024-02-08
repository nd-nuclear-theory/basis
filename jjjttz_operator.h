/************************************************************//**
  @file jjjt_operator.h

  Defines functions for I/O and manipulation of two-body operator
  matrices in jjJT coupling scheme.  Written for use in Moshinsky
  transformation.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 07/13/16 (mac): Created, adapting code from lsjt_operator.
  + 11/04/16 (mac): Remove dependency on Eigen/Core.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.

****************************************************************/

#ifndef BASIS_JJJTTz_OPERATOR_H_
#define BASIS_JJJTTz_OPERATOR_H_

#include <array>
#include <iosfwd>


#include "jt_operator.h"
#include "jjjttz_scheme.h"
#include "many_body.h"
#include "operator.h"

namespace basis {

  void WriteTwoBodyOperatorJJJTTz(
      std::ostream& os,
      const basis::TwoBodySectorsJJJTTz& sectors,
      const basis::OperatorBlocks<double>& matrices,
      basis::NormalizationConversion conversion_mode
    );
  // Write two-body operator in JJJTTz
  // scheme.
  //
  // Side effect: The floating point precision attribute of the output
  // stream is modified.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   sector (basis::TwoBodySectorsJJJTTz) : sectors defining operator
  //   matrices (basis::OperatorBlocks<double>) : matrices defining operator
  //   conversion (basis::NormalizationConversion) : specifies any
  //     conversion between AS and NAS for output

  void SetTwoBodyOperatorMatrixElementJJJTTz(
      const basis::TwoBodySpaceJJJTTz& space, // TODO: add accessor to space in basis.h and remove it here
      const TwoBodySubspaceJJJTTzLabels& subspace_labels_bra,
      const TwoBodySubspaceJJJTTzLabels& subspace_labels_ket,
      const TwoBodyStateJJJTTzLabels& state_labels_bra,
      const TwoBodyStateJJJTTzLabels& state_labels_ket,
      const TwoBodySectorsJJJTTz& sectors,
      OperatorBlocks<double>& matrices,
      double matrix_element
    );

  double GetTwoBodyOperatorMatrixElementJJJTTz(
      const basis::TwoBodySpaceJJJTTz& space, // TODO: add accessor to space in basis.h and remove it here
      const TwoBodySubspaceJJJTTzLabels& subspace_labels_bra,
      const TwoBodySubspaceJJJTTzLabels& subspace_labels_ket,
      const TwoBodyStateJJJTTzLabels& state_labels_bra,
      const TwoBodyStateJJJTTzLabels& state_labels_ket,
      const TwoBodySectorsJJJTTz& sectors,
      const OperatorBlocks<double>& matrices
    );




  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
