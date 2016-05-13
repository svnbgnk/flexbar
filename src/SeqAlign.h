/*
 *   SeqAlign.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_SEQALIGN_H
#define FLEXBAR_SEQALIGN_H


template <typename TSeqStr, typename TString, class TAlgorithm>
class SeqAlign {

private:
	
	const flexbar::TrimEnd m_trimEnd;
	const flexbar::LogAlign m_log;
	const flexbar::FileFormat m_format;
	
	const bool m_isBarcoding, m_writeTag, m_randTag, m_strictRegion;
	const int m_minLength, m_minOverlap, m_tailLength;
	const float m_threshold;
	
	const unsigned int m_bundleSize;
	
	tbb::atomic<unsigned long> m_nPreShortReads, m_modified;
	
	tbb::concurrent_vector<flexbar::TBar> *m_queries;
	tbb::concurrent_vector<unsigned long> *m_rmOverlaps;
	
	std::ostream *m_out;
	TAlgorithm *algo;
	
public:
	
	SeqAlign(tbb::concurrent_vector<flexbar::TBar> *queries, const Options &o, int minOverlap, float threshold, const int tailLength, const int match, const int mismatch, const int gapCost, const flexbar::TrimEnd end, const bool isBarcoding):
			
			m_minOverlap(minOverlap),
			m_threshold(threshold),
			m_tailLength(tailLength),
			m_trimEnd(end),
			m_isBarcoding(isBarcoding),
			m_randTag(o.randTag),
			m_minLength(o.min_readLen),
			m_log(o.logAlign),
			m_format(o.format),
			m_writeTag(o.useRemovalTag),
			m_strictRegion(! o.relaxRegion),
			m_bundleSize(o.bundleSize),
			m_out(o.out){
		
		m_queries = queries;
		
		m_nPreShortReads = 0;
		m_modified       = 0;
		
		algo = new TAlgorithm(o, match, mismatch, gapCost, m_trimEnd);
		
		m_rmOverlaps = new tbb::concurrent_vector<unsigned long>(flexbar::MAX_READLENGTH + 1, 0);
	};
	
	
	virtual ~SeqAlign(){
		delete algo;
		delete m_rmOverlaps;
	};
	
	
	int alignSeqRead(void* item, const bool performRemoval, flexbar::TAlignments &alignments, flexbar::ComputeCycle cycle, unsigned int &aIdx){
		
		using namespace std;
		using namespace flexbar;
		
		using seqan::prefix;
		using seqan::suffix;
		using seqan::infix;
		
		
		SeqRead<TSeqStr, TString> &seqRead = *static_cast< SeqRead<TSeqStr, TString>* >(item);
		
		TSeqStr &readSeq = seqRead.seq;
		int readLength   = length(readSeq);
		
		if(! m_isBarcoding && readLength < m_minLength){
			
			if(cycle != PRECYCLE) ++m_nPreShortReads;
			return 0;
		}
		
		
		if(cycle == PRECYCLE){
			
			unsigned int nQueries = m_queries->size();
			
			if(aIdx == 0){
				reserve(alignments.first, m_bundleSize * nQueries);
			}
			
			for(unsigned int i = 0; i < nQueries; ++i){
				
				TAlign align;
				resize(rows(align), 2);
				assignSource(row(align, 0), readSeq);
				assignSource(row(align, 1), m_queries->at(i).first->seq);
				
				appendValue(alignments.first, align);
				
				++aIdx;
			}
			return 0;
		}
		
		
		int fmismatches, fgapsR, fgapsA, foverlapLength, fqueryLength, ftailLength;
		int fstartPos, fstartPosA, fstartPosS, fendPos, fendPosS, fendPosA;
		
		int qIndex   = -1;
		int scoreMax = -1000000;
		
		float fallowedErrors;
		
		TSeqStr finalRandTag;
		TString finalAlStr;
		
		TString &readTag = seqRead.tag;
		TSeqStr sequence = readSeq;
		
		
		// align each query sequence and keep track of best one
		for(unsigned int i = 0; i < m_queries->size(); ++i){
			
			if(i > 0) cycle = RESULTS;
			
			TSeqStr &query = m_queries->at(i).first->seq;
			
			int queryLength = length(query);
			int tailLength  = (m_tailLength > 0) ? m_tailLength : queryLength;
			
			
			if(m_trimEnd == LEFT_TAIL || m_trimEnd == RIGHT_TAIL){
				
				if(tailLength < readLength){
					
					if(m_trimEnd == LEFT_TAIL){
						sequence = prefix(readSeq, tailLength);
					}else{
						sequence = suffix(readSeq, readLength - tailLength);
					}
				}
			}
			
			
			int startPos = 0, endPos = 0, startPosA = 0, endPosA = 0, startPosS = 0, endPosS = 0;
			int alScore = 0, mismatches = 0, gapsR = 0, gapsA = 0;
			
			TSeqStr randTag;
			stringstream alString;
			
			// align query to read sequence
			algo->alignGlobal(query, sequence, gapsR, gapsA, mismatches, startPos, endPos, startPosA, endPosA,
			                  startPosS, endPosS, alScore, alString, randTag, alignments, cycle, aIdx);
			
			
			int overlapLength = endPos - startPos;
			
			float allowedErrors = m_threshold * overlapLength / 10.0f;
			float madeErrors    = static_cast<float>(mismatches + gapsR + gapsA);
			
			int minOverlapValue = (m_isBarcoding && m_minOverlap == 0) ? queryLength : m_minOverlap;
			
			
			bool validAli = true;
			
			if(((m_trimEnd == RIGHT_TAIL || m_trimEnd == RIGHT) && startPosA < startPosS && m_strictRegion) ||
			   ((m_trimEnd == LEFT_TAIL  || m_trimEnd == LEFT)  && endPosA   > endPosS   && m_strictRegion) ||
			     overlapLength < 1){
				
				validAli = false;
			}
			
			
			// check if alignment is valid and score is max as well as if number of errors and overlap length are allowed
			if(validAli && alScore > scoreMax && madeErrors <= allowedErrors && overlapLength >= minOverlapValue){
				
				qIndex      = i;
				scoreMax    = alScore;
				fstartPos   = startPos;
				fstartPosA  = startPosA;
				fstartPosS  = startPosS;
				fendPos     = endPos;
				fendPosA    = endPosA;
				fendPosS    = endPosS;
				fgapsR      = gapsR;
				fgapsA      = gapsA;
				
				ftailLength    = tailLength;
				foverlapLength = overlapLength;
				fqueryLength   = queryLength;
				
				if(m_randTag) finalRandTag = randTag;
				
				if(m_log != NONE){
					fmismatches    = mismatches;
					finalAlStr     = alString.str();
					fallowedErrors = allowedErrors;
				}
			}
		}
		
		
		// valid alignment
		if(qIndex >= 0){
			
			TrimEnd trimEnd = m_trimEnd;
			
			// cut read according to best alignment
			if(performRemoval){
				
				if(trimEnd == ANY){
					
					if(fstartPosA <= fstartPosS && fendPosS <= fendPosA){
						seqRead.seq = "";
						if(m_format == FASTQ) seqRead.qual = "";
					}
					else if(fstartPosA - fstartPosS >= fendPosS - fendPosA){
						trimEnd = RIGHT;
					}
					else{
						trimEnd = LEFT;
					}
				}
				
				switch(trimEnd){
					
					int rCutPos;
					
					case LEFT_TAIL:
						sequence = readSeq;
					
					case LEFT:
						rCutPos = fendPos;
						
						// translate alignment end pos to read idx
						if(fstartPosS > 0) rCutPos -= fstartPosS;
						
						// adjust to inner read gaps
						rCutPos -= fgapsR;
						
						if(rCutPos > readLength) rCutPos = readLength;
						
						erase(sequence, 0, rCutPos);
						seqRead.seq = sequence;
						
						if(m_format == FASTQ){
							erase(seqRead.qual, 0, rCutPos);
						}
						break;
					
					case RIGHT_TAIL:
						sequence  = readSeq;
						// adjust cut pos to original read length
						fstartPos += readLength - ftailLength;
					
					case RIGHT:
						rCutPos = fstartPos;
						
						// skipped restriction
						if(rCutPos < 0) rCutPos = 0;
						
						erase(sequence, rCutPos, readLength);
						seqRead.seq = sequence;
						
						if(m_format == FASTQ){
							erase(seqRead.qual, rCutPos, readLength);
						}
						break;
						
	                case ANY:;
				}
				
				++m_modified;
				
				
				// count for each query number of removals
				m_queries->at(qIndex).second.first++;
				
				if(foverlapLength == fqueryLength){
					m_queries->at(qIndex).second.second++;
				}
				
				if(m_writeTag){
					append(seqRead.tag, "_Flexbar_removal");
					
					if(! m_isBarcoding){
						append(seqRead.tag, "_");
						append(seqRead.tag, m_queries->at(qIndex).first->tag);
					}
				}
				
				// store overlap occurrences for min, max, mean and median
				if(foverlapLength <= MAX_READLENGTH) m_rmOverlaps->at(foverlapLength)++;
				else cerr << "\nCompile Flexbar with larger max read length to get correct overlap stats.\n" << endl;
			}
			
			
			// valid alignment, not neccesarily removal
			
			if(m_randTag && finalRandTag != ""){
				append(seqRead.tag, "_");
				append(seqRead.tag, finalRandTag);
			}
			
			
			// alignment stats
			TString &queryTag = m_queries->at(qIndex).first->tag;
			
			if(m_log == ALL || (m_log == MOD && performRemoval)){
				
				stringstream ss;
				
				if(performRemoval){
					ss << "Sequence removal:";
					
					     if(trimEnd == LEFT  || trimEnd == LEFT_TAIL)  ss << "  left side\n";
					else if(trimEnd == RIGHT || trimEnd == RIGHT_TAIL) ss << "  right side\n";
					else                                               ss << "  any side\n";
				}
				else{ ss << "Sequence detection, no removal:\n"; }
				
				ss << "  query tag        " << queryTag                      << "\n"
				   << "  read tag         " << readTag                       << "\n"
				   << "  read seq         " << readSeq                       << "\n"
				   << "  read pos         " << fstartPosS << "-" << fendPosS << "\n"
				   << "  query pos        " << fstartPosA << "-" << fendPosA << "\n"
				   << "  score            " << scoreMax                      << "\n"
				   << "  overlap          " << foverlapLength                << "\n"
				   << "  errors           " << fgapsR + fgapsA + fmismatches << "\n"
				   << "  allowed errors   " << fallowedErrors                << "\n";
				
				if(performRemoval){
					ss << "  remaining read   "  << seqRead.seq << "\n";
					
					if(m_format == FASTQ)
					ss << "  remaining qual   " << seqRead.qual << "\n";
				}
				
				ss << "\n  Alignment:\n" << endl << finalAlStr;
				
				*m_out << ss.str();
			}
			else if(m_log == TAB){
				
				stringstream ss;
				
				ss << readTag     << "\t" << queryTag        << "\t"
				   << fstartPosA  << "\t" << fendPosA        << "\t" << foverlapLength << "\t"
				   << fmismatches << "\t" << fgapsR + fgapsA << "\t" << fallowedErrors << endl;
				
				*m_out << ss.str();
			}
		}
		else if(m_log == ALL){
			
			stringstream ss;
			
			ss << "No valid alignment:"   << "\n"
			   << "read tag  " << readTag << "\n"
			   << "read seq  " << readSeq << "\n\n" << endl;
			
			*m_out << ss.str();
		}
		
		return ++qIndex;
	}
	
	
	std::string getOverlapStatsString(){
		
		using namespace flexbar;
		
		unsigned long nValues = 0, halfValues = 0, cumValues = 0, lenSum = 0;
		int min = 1000000, max = 0, median = 0, mean = 0;
		
		for (int i = 0; i <= MAX_READLENGTH; ++i){
			unsigned long lenCount = m_rmOverlaps->at(i);
			
			if(lenCount > 0 && i < min) min = i;
			if(lenCount > 0 && i > max) max = i;
			
			nValues += lenCount;
			lenSum  += lenCount * i;
		}
		
		halfValues = nValues / 2;
		
		for (int i = 0; i <= MAX_READLENGTH; ++i){
			cumValues += m_rmOverlaps->at(i);
			
			if(cumValues >= halfValues){
				median = i;
				break;
			}
		}
		
		if(m_modified > 0) mean = lenSum / m_modified;
		
		std::stringstream ss;
		
		ss << "Min, max, mean and median adapter overlap: ";
		ss << min << " / " << max << " / " << mean << " / " << median;
		
		return ss.str();
	}
	
	
	unsigned long getNrPreShortReads() const {
		return m_nPreShortReads;
	}
	
	
	unsigned long getNrModifiedReads() const {
		return m_modified;
	}
	
};

#endif
