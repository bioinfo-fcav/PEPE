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

	
	std::string align(Dna5String seqV, Dna5String seqH, int alnType, int match, int mismatch, int gapExtendScore, int gapOpenScore)
	{
		typedef int TValue;
		typedef Score<TValue, ScoreMatrix<Dna5, UserDefinedMatrix> > TScoringScheme;
		//TScoringScheme userScoringSchemeDna(gapExtendScore, gapOpenScore);
		Score<int, Simple> userScoringSchemeDna(match, mismatch, gapExtendScore, gapOpenScore);

		Align<Dna5String> align;
		resize(rows(align), 2);
		assignSource(row(align, 0), seqH);
		assignSource(row(align, 1), seqV);
		
		int result;
		switch (alnType) {
			case 1:
				result = localAlignment(align, userScoringSchemeDna);
				break;
			case 2:
				result = globalAlignment(align, userScoringSchemeDna, AlignConfig<false, false, false, false>());
				// ordinary global alignment
				break;
			case 3:
				result = globalAlignment(align, userScoringSchemeDna, AlignConfig<true, false, false, true>());
				//semiglobal alignment, free begin and end gaps in second/vertical sequence
				break;
			case 4:
				result = globalAlignment(align, userScoringSchemeDna, AlignConfig<false, true, true, false>());
				//semiglobal alignment, free begin and end gaps in first/horizontal sequence
				break;
			case 5:
				result = globalAlignment(align, userScoringSchemeDna, AlignConfig<false, true, false, true>());
				//overlap alignment with second/vertical sequence overhanging to the left of first/horizontal
				break;
			case 6:
				result = globalAlignment(align, userScoringSchemeDna, AlignConfig<true, false, true, false>());
				//overlap alignment with first/horizontal sequence overhanging to the left of second/vertical
				break;
			case 7:
				result = globalAlignment(align, userScoringSchemeDna, AlignConfig<false, true, false, false>());
				//free begin gaps in second/vertical sequence only
				break;
			case 8:
				result = globalAlignment(align, userScoringSchemeDna, AlignConfig<false, false, true, false>());
				//free end gaps in second/vertical sequence only
				break;
			default:
				result = globalAlignment(align, userScoringSchemeDna, AlignConfig<true, false, false, true>());
				//semiglobal alignment, free begin and end gaps in second/vertical sequence
				break;
		}
		//std::cout << "Score: " << result << "\n";
		//std::cout << "The resulting alignment is\n"
		//	<< row(align,0) << "\n" << row(align,1) << "\n";
		std::stringstream buffer;
		buffer << result << ";" << row(align,0) << ";" << row(align,1) << ";" << (beginPosition(row(align,1))+1) << ";" << (endPosition(row(align,1))+1);
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

	char * Align2Seq(char * seqV, char * seqH, int alnType=0, int match=2, int mismatch=-1, int gapExtendScore=-1, int gapOpenScore=-2) {
		std::string str = align(seqV, seqH, alnType, match, mismatch, gapExtendScore, gapOpenScore);

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
MyClass::Align2Seq(char * seqV, char * seqH, int alnType=0, int match=2, int mismatch=-1, int gapExtendScore=-1, int gapOpenScore=-2)

void
MyClass::DESTROY()
