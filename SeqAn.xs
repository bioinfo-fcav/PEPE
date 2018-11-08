#include <iostream>
#include <seqan/align.h>

using namespace seqan;

namespace seqan {
	// We have to create a new specialization of the ScoringMatrix_ class
	// for the DNA alphabet. For this, we first create a new tag.
	struct UserDefinedMatrix {};
	// Then, we specialize the class ScoringMatrix_ for the Dna5 alphabet.
	template <>
		struct ScoringMatrixData_<int, Dna5, UserDefinedMatrix>
		{
			enum
			{
				VALUE_SIZE = ValueSize<Dna5>::VALUE,
				TAB_SIZE = VALUE_SIZE * VALUE_SIZE
			};
			static inline int const * getData()
			{
				static int const _data[TAB_SIZE] =
				{
					 2, -1, -1, -1, 1,
					-1,  2, -1, -1, 1,
					-1, -1,  2, -1, 1,
					-1, -1, -1,  2, 1,
					 1,  1,  1,  1, 2 
				};
				return _data;
			}
		};


	//seqan::String<char> align(Dna5String seqV, Dna5String seqH)
	std::string align(Dna5String seqV, Dna5String seqH)
	{
		typedef int TValue;
		typedef Score<TValue, ScoreMatrix<Dna5, UserDefinedMatrix> > TScoringScheme;
		int const gapOpenScore = -2;
		int const gapExtendScore = -1;
		TScoringScheme userScoringSchemeDna(gapExtendScore, gapOpenScore);
		
		Align<Dna5String> align;
		resize(rows(align), 2);
		assignSource(row(align, 0), seqH);
		assignSource(row(align, 1), seqV);
		
		//AlignConfig<TTop, TLeft, TRight, TDown>
		AlignConfig<true, false, false, true> alignConfig;

		int result = globalAlignment(align, userScoringSchemeDna, alignConfig);

		//std::cout << "Score: " << result << "\n";
		//std::cout << "The resulting alignment is\n"
		//	<< row(align,0) << "\n" << row(align,1) << "\n";
		std::stringstream buffer;
		buffer << result << ";" << row(align,0) << ";" << row(align,1);
		return buffer.str();
	}

} // namespace seqan

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef __cplusplus
}
#endif


class MyClass {
public:
	MyClass() {

	}

	char * Align2Seq(char * seqV, char * seqH ) {
		std::string str = align(seqV, seqH);

		char* chr = strdup(str.c_str());
		return chr;
	}
	
	~MyClass() { 
			//std::cout << "Destruction is a way of life for me." << std::endl; 
		   }
};

MODULE = PEPE::SeqAn		PACKAGE = PEPE::SeqAn

MyClass *
MyClass::new();

char *
MyClass::Align2Seq(char * seqV, char * seqH)

void
MyClass::DESTROY()
