/*
 *   SeqAlignAlgorithm.h
 *
 */

#ifndef FLEXBAR_SEQALIGNALGORITHM_H
#define FLEXBAR_SEQALIGNALGORITHM_H


template <typename TSeqStr>
class SeqAlignAlgorithm {

private:
	
	typedef typename seqan::Value<TSeqStr>::Type        TChar;
	typedef typename seqan::Row<flexbar::TAlign>::Type  TRow;
	typedef typename seqan::Iterator<TRow>::Type        TRowIterator;
	
	typedef seqan::Score<int, seqan::Simple>                              TScoreSimple;
	typedef seqan::Score<int, seqan::ScoreMatrix<TChar, seqan::Default> > TScoreMatrix;
	
	TScoreSimple m_score;
	TScoreMatrix m_scoreMatrix;
	
	const bool m_randTag;
	const flexbar::LogAlign m_log;
	const flexbar::TrimEnd m_trimEnd;
	
public:
	
	SeqAlignAlgorithm(const Options &o, const int match, const int mismatch, const int gapCost, const flexbar::TrimEnd trimEnd):
			m_randTag(o.randTag),
			m_log(o.logAlign),
			m_trimEnd(trimEnd){
		
		using namespace seqan;
		
		m_score = Score<int, Simple>(match, mismatch, gapCost);
		
		m_scoreMatrix = TScoreMatrix(gapCost);
		
		for(unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i){
			for(unsigned j = 0; j < ValueSize<TChar>::VALUE; ++j){
				
				if(i == j || TChar(j) == 'N')
					 setScore(m_scoreMatrix, TChar(i), TChar(j), match);
				else setScore(m_scoreMatrix, TChar(i), TChar(j), mismatch);
			}
		}
		
		// printScoreMatrix(m_scoreMatrix);
	};
	
	
	virtual ~SeqAlignAlgorithm(){
	};
	
	
	void alignGlobal(const TSeqStr &querySeq, const TSeqStr &readSeq, int &gapsR, int &gapsA, int &mismatches, int &startPos, int &endPos, int &startPosA, int &endPosA, int &startPosS, int &endPosS, int &alScore, std::stringstream &alString, TSeqStr &tagSeq, flexbar::TAlignments &alignments, const flexbar::ComputeCycle cycle, unsigned int &aIdx){
		
		using namespace std;
		using namespace seqan;
		using namespace flexbar;
		
		
		// TAlign align2;
		// TAlign &align = align2;
		//
		// if(m_randTag){
		//
		// 	resize(rows(align), 2);
		// 	assignSource(row(align, 0), readSeq);
		// 	assignSource(row(align, 1), querySeq);
		//
		// 	if(m_trimEnd == RIGHT || m_trimEnd == RIGHT_TAIL){
		//
		// 		AlignConfig<true, false, true, true> ac;
		// 		alScore = globalAlignment(align, m_scoreMatrix, ac);
		// 	}
		// 	else if(m_trimEnd == LEFT || m_trimEnd == LEFT_TAIL){
		//
		// 		AlignConfig<true, true, false, true> ac;
		// 		alScore = globalAlignment(align, m_scoreMatrix, ac);
		// 	}
		// 	else{
		// 		AlignConfig<true, true, true, true> ac;
		// 		alScore = globalAlignment(align, m_scoreMatrix, ac);
		// 	}
		// }
		// else{
			if(cycle == COMPUTE){
				
				if(m_trimEnd == RIGHT || m_trimEnd == RIGHT_TAIL){
					
					AlignConfig<true, false, true, true> ac;
					alignments.second = globalAlignment(alignments.first, m_score, ac);
				}
				else if(m_trimEnd == LEFT || m_trimEnd == LEFT_TAIL){
					
					AlignConfig<true, true, false, true> ac;
					alignments.second = globalAlignment(alignments.first, m_score, ac);
				}
				else{
					AlignConfig<true, true, true, true> ac;
					alignments.second = globalAlignment(alignments.first, m_score, ac);
				}
			}
			
			TAlign &align = value(alignments.first,  aIdx);
			alScore       = value(alignments.second, aIdx);
		// }
		
		// cout << "Score: " << alScore << endl;
		// cout << "Align: " << align << endl;
		
		
		TRow &row1 = row(align, 0);
		TRow &row2 = row(align, 1);
		
		startPosS = toViewPosition(row1, 0);
		startPosA = toViewPosition(row2, 0);
		endPosS   = toViewPosition(row1, length(source(row1)));
		endPosA   = toViewPosition(row2, length(source(row2)));
		
		if(startPosA > startPosS) startPos = startPosA;
		else                      startPos = startPosS;
		
		if(endPosA > endPosS) endPos = endPosS;
		else                  endPos = endPosA;
		
		
		// cout << "\n\n" << startPosS << endl << startPosA << endl << endPosS << endl << endPosA;
		// cout << align << endl << alScore << endl;
		
		if(m_log != flexbar::NONE) alString << align;
		
		
		TRowIterator it1 = begin(row1);
		TRowIterator it2 = begin(row2);
		
		int aliPos = 0;
		gapsR      = 0;
		gapsA      = 0;
		mismatches = 0;
		
		for(; it1 != end(row1); ++it1){
			
			if(startPos <= aliPos && aliPos < endPos){
				     if(isGap(it1))                   ++gapsR;
				else if(isGap(it2))                   ++gapsA;
				else if(*it1 != *it2 && *it2 != 'N')  ++mismatches;
				else if(m_randTag    && *it2 == 'N')  append(tagSeq, (TChar) *it1);
			}
			++aliPos;
			++it2;
		}
		
		// cout << "\n\n" << gapsR << endl << gapsA << endl << mismatches << endl;
		
		++aIdx;
	}
	
	
	void printScoreMatrix(TScoreMatrix &scoreMatrix){
		
		using namespace std;
		using namespace seqan;
		
		cout << endl;
		for(unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i)
			cout << "\t" << TChar(i);
		cout << endl;
		
		for(unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i){
			cout << TChar(i);
			for(unsigned j = 0; j < ValueSize<TChar>::VALUE; ++j)
				cout << "\t" << score(scoreMatrix, TChar(i), TChar(j));
			cout << endl;
		}
	}
};


#endif