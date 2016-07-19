/****************************************************************
  jjjpnorb_operator.cpp
                                 
  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "jjjpnorb_operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body jjJpn operator output
  ////////////////////////////////////////////////////////////////

  void WriteTwoBodyOperatorJJJPN(
      std::ostream& os,
      const basis::TwoBodySectorsJJJPN& sectors,
      const basis::MatrixVector& matrices,
      basis::NormalizationConversion conversion_mode,
      int indexing_base
    );


//   void OutMFDnH2Stream::WriteSectorText (const TwoBodyMatrixNljTzJP& matrix)
//   {
//     // recover sector properties
//     const TwoSpeciesStateType state_type = sector_.GetStateType();
//     const int J = sector_.GetJ();
//     const int g = sector_.GetGrade();
//     const int dimension = matrix.GetTwoBodyBasis().GetDimension(state_type,J,g);
// 		
//     // output floating point format setup
//     const int field_width = 3;
//     const int float_precision = 7; // precision 7 for 8 digits total
//     *stream_ << std::scientific << std::setprecision(float_precision);
// 		
//     // for canonical pairs of states (in lexicographical order)
//     for (int k1 = 0; k1 < dimension; ++k1)
//       for (int k2 = k1; k2 < dimension; ++k2)
// 	{
// 	  TwoBodyStateNlj s1 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1);
// 	  TwoBodyStateNlj s2 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2);
// 				
// 	  int twice_J = 2*J;
// 				
// 	  *stream_ << " " << std::setw(field_width) << s1.a1.GetIndex1() 
// 		   << " " << std::setw(field_width) << s1.a2.GetIndex1()
// 		   << " " << std::setw(field_width) << s2.a1.GetIndex1() 
// 		   << " " << std::setw(field_width) << s2.a2.GetIndex1()
// 		   << " " << std::setw(field_width) << twice_J 
// 		   << " " << std::setw(field_width) << TwoSpeciesStateTypeToOnsager(state_type) 
// 		   << " " << std::setw(16) << matrix.GetMatrixElementNAS(state_type, s1, s2) 
// 		   << std::endl;
// 	}
//   }
// 




  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
