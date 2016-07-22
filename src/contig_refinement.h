#ifndef _CONTIG_REFINEMENT_
#define _CONTIG_REFINEMENT_

#include "align_batch.h"
#include "scoring_prot.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
#include <stack>
#include <queue>
#include <deque>
#include <vector>
#include <tuple>
#include <string>
#include <list>



struct AlignmentPrintType {
  std::string seq1, seq2, symbol, header;
  std::unordered_map<int, int> nuc_match;
};

struct ContigType {
  std::string sequence;
  int score; 
  double bit_score;
  double e_value;
  int q_begin, q_end;
  bool valid;
  AlignmentPrintType al_print;
};

struct MerPosType {
  int cid;
  int pos;
  int fw_len;
  int re_len;
};

struct ContigOverlapType  {
  int cid;
  int count;
  int ref_begin, ref_end;
  int target_begin, target_end;
};

class ContigRefinement  {
 public:
  explicit ContigRefinement(void) {}
  ~ContigRefinement(void) {}

  ContigType MakeContigType(
      std::string &seq, const int score, const double bit_score, const double e_value
  );

  void RefineContigs(
      std::string& query, const long int db_size, ScoringProt &score_obj, 
      const int band_size, const int refine_cutoff, const double e_value, const int mer_len, 
      std::list<ContigType>& contigs, std::list<ContigType>& refined_contigs
  );
 private:
  std::vector<ContigType> contig_holder_;
  void IndexContigs(
      const int mer_len, std::unordered_map<std::string, std::list<MerPosType> >& mer_contig 
  );
  void IncorporateContig(
      int mer_len, std::unordered_map<std::string, std::list<MerPosType> >& kmer_map,
      int contig_index, ContigType& current_contig, std::list<int>& incorporated_contigs
  );
  bool TryMergeOverlapedContigs(
      int mer_len, ContigOverlapType& overlap_record, 
      ContigType& ref_contig, int& ref_index
  );
  bool FindConsensusTail(
      std::string& seq1, std::string& seq2, std::string& consensus
  );
  bool FindConsensusHead(
      std::string& seq1, std::string& seq2, std::string& consensus
  );
  void MergeAllContigs(
      std::string query, const long int db_size, 
      ScoringProt &score_obj, int band_size, const int refine_cutoff, double e_value,
      int mer_len, std::unordered_map<std::string, std::list<MerPosType> >& mer_contig,
      std::list<ContigType>& refined_contigs
  );
  std::string FixedWidthString(int len, int num);
};

#endif
