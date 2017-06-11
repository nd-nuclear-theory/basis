/************************************************************//**
  @file multibasis_test.cpp

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "basis/multibasis.h"

#include "am/am.h"


namespace basis {

  ////////////////////////////////////////////////////////////////
  // relative states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // As for RelativeLSJT (in lsjt_scheme.h) but with...
  //
  // multiplicity labels: (M_J,M_T)
  //
  //   => multiplicity (2*J+1)*(2*T+1)
  //
  ////////////////////////////////////////////////////////////////

  // labels

  typedef std::tuple<int,int,int,int,int> RelativeSubspaceLSJTLabels;
  typedef std::tuple<int> RelativeStateLSJTLabels;

  // subspace

  class RelativeMultiSubspaceLSJT
    : public BaseMultiSubspace<RelativeSubspaceLSJTLabels,RelativeStateLSJTLabels>
    {

      public:

      // constructor

      RelativeMultiSubspaceLSJT() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      RelativeMultiSubspaceLSJT(int L, int S, int J, int T, int g, int Nmax);
      // Set up indexing.

      // accessors

      int L() const {return std::get<0>(labels_);}
      int S() const {return std::get<1>(labels_);}
      int J() const {return std::get<2>(labels_);}
      int T() const {return std::get<3>(labels_);}
      int g() const {return std::get<4>(labels_);}
      int Nmax() const {return Nmax_;}

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      //validation
      bool ValidLabels() const;

      // truncation
      int Nmax_;

    };

  // state

  class RelativeMultiStateLSJT
    : public BaseMultiState<RelativeMultiSubspaceLSJT>
  {

    public:

    // pass-through constructors

    RelativeMultiStateLSJT(const SubspaceType& subspace, int index)
      // Construct state by index.
      : BaseMultiState<RelativeMultiSubspaceLSJT>(subspace,index) {}

    RelativeMultiStateLSJT(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseMultiState<RelativeMultiSubspaceLSJT>(subspace,state_labels) {}

    // pass-through accessors
    int L() const {return subspace().L();}
    int S() const {return subspace().S();}
    int J() const {return subspace().J();}
    int T() const {return subspace().T();}
    int g() const {return subspace().g();}

    // state label accessors
    int N() const {return std::get<0>(labels());}

  };

  // space

  class RelativeMultiSpaceLSJT
    : public BaseMultiSpace<RelativeMultiSubspaceLSJT>
  {

    public:

    // constructor

    RelativeMultiSpaceLSJT() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    RelativeMultiSpaceLSJT(int Nmax, int Jmax);
    // Enumerate all relative LSJT subspaces of given dimension up to
    // a given relative oscillator cutoff and relative angular
    // momentum cutoff.
    //
    // The relative angular momentum cutoff is included in recognition
    // of the common practice of truncation by highest partial wave in
    // the representation of relative interations.
    //
    // Arguments:
    //   Nmax (int) : relative oscillator truncation on included subspaces
    //   Jmax (int) : relative angular momentum truncation on included
    //     subspaces (Jmax<=Nmax+1)

    // accessors
    int Nmax() const {return Nmax_;}
    int Jmax() const {return Jmax_;}

    // diagnostic string
    std::string DebugStr() const;

    private:
    // truncation
    int Nmax_, Jmax_;

  };


  ////////////////////////////////////////////////////////////////
  // implementation
  ////////////////////////////////////////////////////////////////

  RelativeMultiSubspaceLSJT::RelativeMultiSubspaceLSJT(int L, int S, int J, int T, int g, int Nmax)
  {

    // set values
    labels_ = SubspaceLabelsType(L,S,J,T,g);
    Nmax_ = Nmax;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing

    // manual version w/o lookup
    //	  dimension_ = (Nmax() - g()) / 2 + 1;

    // iterate over total oscillator quanta
    for (int N = L; N <= Nmax; N +=2)
      PushStateLabels(StateLabelsType(N),(2*J+1)*(2*T+1));

  }

  bool RelativeMultiSubspaceLSJT::ValidLabels() const
  {

    bool valid = true;

    // triangularity
    valid &= am::AllowedTriangle(L(),S(),J());
    // parity
    valid &= (g() == (L()%2));
    // antisymmetry
    valid &= ((L()+S()+T())%2 == 1);
    // truncation
    // valid &= ((Nmax()%2)==g());

    return valid;
  }

  std::string RelativeMultiSubspaceLSJT::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << L() 
       << " " << std::setw(width) << S() 
       << " " << std::setw(width) << J() 
       << " " << std::setw(width) << T() 
       << " " << std::setw(width) << g()
       << " " << "]";

    return os.str();
  }

  std::string RelativeMultiSubspaceLSJT::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        RelativeMultiStateLSJT state(*this,state_index);

        os
	  << " " << "index"
	  << " " << std::setw(width) << state_index
	  << " " << "N"
	  << " " << std::setw(width) << state.N()
	  << " " << "multiplicity"
	  << " " << std::setw(width) << state.multiplicity()
          << std::endl;
      }

    os << "dimensions: size " << size() << " full_dimension " << full_dimension() << std::endl;

    return os.str();

  }

  RelativeMultiSpaceLSJT::RelativeMultiSpaceLSJT(int Nmax, int Jmax)
    : Nmax_(Nmax), Jmax_(Jmax)
  {

    // iterate over L
    for (int L=0; L<=Nmax; ++L)
      {
        // set g
	int g = L%2;

	// iterate over S
	for (int S=0; S<=1; ++S)
	  {
            // set T
	    int T = (L+S+1)%2;

	    // iterate over J
            int J_limit = std::min(L+S,Jmax);
	    for (int J=abs(L-S); J<=J_limit; ++J)
	      {
		// downshift Nmax to match parity of subspace
		// required to pass label validity tests
		// int Nmax_subspace = Nmax - (Nmax-g)%2;

		RelativeMultiSubspaceLSJT subspace(L,S,J,T,g,Nmax);
		assert(subspace.size()!=0);
		PushSubspace(subspace);
	      }
	  }
      }
  }

  std::string RelativeMultiSpaceLSJT::DebugStr() const
  {
    std::ostringstream os;

    const int width = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);
	os
	  << " " << "index"
	  << " " << std::setw(width) << subspace_index
	  << " " << " (L,S,J,T,g) "
	  << " " << subspace.LabelStr()
	  << " " << "Nmax"
	  << " " << std::setw(width) << subspace.Nmax()
	  << " " << "dim"
	  << " " << std::setw(width) << subspace.size()
	  << " " << std::setw(width) << subspace.full_dimension()
	  << " " << std::endl;
      }

    return os.str();
  }



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

void Test()
{
  std::cout << "Relative space" << std::endl;
  int N_max = 2;
  int J_max = 3;

  basis::RelativeMultiSpaceLSJT space(N_max,J_max);
  std::cout << space.DebugStr() << std::endl;
  std::cout
    << "Dimension " << space.Dimension()
    << "FullDimension " << space.FullDimension() << std::endl;
  std::cout << std::endl;

  for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
    {
      const basis::RelativeMultiSubspaceLSJT& subspace = space.GetSubspace(subspace_index);
      std::cout << subspace.LabelStr() << std::endl;
      std::cout << subspace.DebugStr() << std::endl;
    }    
}

int main(int argc, char **argv)
{

  Test();

  // termination
  return EXIT_SUCCESS;
}
