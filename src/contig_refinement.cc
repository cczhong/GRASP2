#include "contig_refinement.h"

using namespace std;

bool _cmp_contig(ContigType& c1, ContigType& c2)  {
  if(c1.score > c2.score)  {
    return true;
  }
  return false;
}

bool _cmp_merpos_fw(MerPosType& m1, MerPosType& m2) {
  if(m1.fw_len > m2.fw_len)  {
    return true;
  }
  return false;
}

bool _cmp_merpos_re(MerPosType& m1, MerPosType& m2) {
  if(m1.re_len > m2.re_len)  {
    return true;
  }
  return false;
}

ContigType ContigRefinement::MakeContigType(
    std::string &seq, const int score, const double bit_score, const double e_value
) {
  ContigType c;
  c.sequence = seq; c.score = score; c.bit_score = bit_score, c.e_value = e_value;
  return c;
}

void ContigRefinement::RefineContigs(
    std::string& query, const long int db_size, ScoringProt &score_obj, 
    const int band_size, const int refine_cutoff, const double e_value, const int mer_len, 
    std::list<ContigType>& contigs, std::list<ContigType>& refined_contigs
) {
  if(contigs.size() <= 0)  {
    return;
  }
  //cout << "RefineContigs called" << endl;
  // sort the contigs based on alignment score
  contigs.sort(_cmp_contig);
  contig_holder_.resize(contigs.size());
  int i = 0;
  for(auto it = contigs.begin(); it != contigs.end(); ++ it) {
    //cout << "contig:  " << it->sequence << endl;
    it->valid = true;
    contig_holder_[i] = *it;
    ++ i;
  }
  // build hash table for the mers
  unordered_map<string, list<MerPosType> > mer_contig;
  IndexContigs(mer_len, mer_contig);
  // merge the contigs based on overlap
  MergeAllContigs(
      query, db_size, score_obj, band_size, refine_cutoff, 
      e_value, mer_len, mer_contig, refined_contigs
  );
  refined_contigs.sort(_cmp_contig);
  //cout << "RefineContigs end" << endl;
  return;
}

void ContigRefinement::IndexContigs(
    const int mer_len, std::unordered_map<std::string, std::list<MerPosType> >& mer_contig 
) {
  //cout << "IndexContigs called" << endl;
  // count all mers in contigs
  for(int i = 0; i < (int) contig_holder_.size(); ++ i) {
    //cout << i << endl;
    for(int j = 0; j < (int) contig_holder_[i].sequence.length() - mer_len; ++ j) {
      //cout << contig_holder_[i].sequence << endl;
      string s = contig_holder_[i].sequence;
      string mer = s.substr(j, mer_len);
      MerPosType mp;
      mp.cid = i; mp.pos = j; 
      mp.fw_len = s.length() - j; mp.re_len = j + mer_len; 
      mer_contig[mer].push_back(mp);
    }
  }
  //cout << "IndexContigs end" << endl;
  return;
}

void ContigRefinement::MergeAllContigs(
    std::string query, const long int db_size, ScoringProt &score_obj, 
    int band_size, const int refine_cutoff, double e_value,
    int mer_len, std::unordered_map<std::string, std::list<MerPosType> >& mer_contig,
    std::list<ContigType>& refined_contigs
) {
  // greedy merge: start with contig with highest score, merge with longest overlap
  int score_cutoff = score_obj.ComputeCutoff(
      query.length(), db_size, e_value
  );
  for(int i = 0; i < (int) contig_holder_.size(); ++ i) {
    if(contig_holder_[i].valid && contig_holder_[i].score >= refine_cutoff)  {
      ContigType c = contig_holder_[i];
      list<int> incorporated_contigs;
      IncorporateContig(mer_len, mer_contig, i, c, incorporated_contigs);
      // redo alignment for the contig
      //cout << "******************" << endl;
      //cout << "q length:  " << query.length() << endl;
      //cout << "index: " << c.q_begin << " " << c.q_end << endl; 
      //cout << "MergeAllContigs: ref:  " << c.sequence << endl;
      //cout << "MergeAllContigs: tgt:  " << q_par << endl;

      //cout << "******************" << endl;     
      //cout << q_par << endl;
      //cout << c.sequence << endl;
      
      AlignBatch aligner;
      pair<string, string> alignment;
      pair<int, int> q_interval, t_interval;
      c.score = aligner.AlignLocalSingle(query, c.sequence, score_obj, alignment, q_interval, t_interval);
      if(c.score >= score_cutoff)  {
        c.sequence = c.sequence.substr(t_interval.first, t_interval.second - t_interval.first + 1);  
        c.q_begin = q_interval.first; c.q_end = q_interval.second;
        c.bit_score = score_obj.ComputeBitScore(c.score);
        c.e_value = score_obj.ComputeEValue(
            query.length(), db_size,
            c.score
        );
        c.al_print.seq1 = alignment.first; c.al_print.seq2 = alignment.second;
        
        refined_contigs.push_back(c); 
        for(auto it_inc = incorporated_contigs.begin(); 
            it_inc != incorporated_contigs.end(); ++ it_inc) {
          contig_holder_[*it_inc].valid = false;
        }
        contig_holder_[i].valid = false;
      }
    }
  }
  return;
}

bool _cmp_contig_overlap(ContigOverlapType& c1, ContigOverlapType& c2)  {
  if(c1.count >= c2.count)  {
    return true;
  }
  return false;
}

void ContigRefinement::IncorporateContig(
    int mer_len, std::unordered_map<std::string, std::list<MerPosType> >& kmer_map,
    int contig_index, ContigType& current_contig, std::list<int>& incorporated_contigs
)  {
  //cout << " IncorporateContig: begin" << endl;
  // count all contigs that share k-mers with the current contig
  unordered_map<int, ContigOverlapType> similar_contigs;
  for(int i = 0; i < (int) current_contig.sequence.length() - mer_len; ++ i) {
    string mer = current_contig.sequence.substr(i, mer_len);
    for(auto it = kmer_map[mer].begin(); it != kmer_map[mer].end(); ++ it) {
      // find if the contig exists
      if(it->cid == contig_index || !contig_holder_[it->cid].valid)  {
        continue;
      }
      auto it_contig = similar_contigs.find(it->cid);
      if(it_contig != similar_contigs.end())  {
        // update the information
        it_contig->second.count ++;
        it_contig->second.ref_end = i;
        it_contig->second.target_end = it->pos;
      } else  {
        // insert the record
        ContigOverlapType cl_record;
        cl_record.cid = it->cid;
        cl_record.count = 1;
        cl_record.ref_begin = cl_record.ref_end = i;
        cl_record.target_begin = cl_record.target_end = it->pos;
        similar_contigs.insert({it->cid, cl_record});
      }
    }
  }
  // filter and sort the overlap contigs
  list<ContigOverlapType> filtered_contigs;
  for(auto it = similar_contigs.begin(); it != similar_contigs.end(); ++ it) {
    if(it->second.ref_end - it->second.ref_begin + 1 == it->second.count &&
       it->second.target_end - it->second.target_begin + 1 == it->second.count)  {
      filtered_contigs.push_back(it->second);
    }
  }
  filtered_contigs.sort(_cmp_contig_overlap);
  // try to merge the contig
  int ref_index = 0;
  for(auto it = filtered_contigs.begin(); it != filtered_contigs.end(); ++ it) {
    //cout << " IncorporateContig:  target contig:" << endl;
    //cout << contig_holder_[it->cid].sequence << endl;
    if(TryMergeOverlapedContigs(mer_len, *it, current_contig, ref_index))  {
      // mark the target contig as invalid
      incorporated_contigs.push_back(it->cid);
    }
    //cout << " IncorporateContig:  merged contig:" << endl;
    //cout << current_contig.sequence << endl;
    //cout << " IncorporateContig:  recruited contig ids: " << endl;
    //for(auto it_cid = incorporated_contigs.begin(); it_cid != incorporated_contigs.end(); ++ it_cid) {
    //  cout << " @@: " << *it_cid << endl; 
    //}
  }
  return;
}

bool ContigRefinement::TryMergeOverlapedContigs(
    int mer_len, ContigOverlapType& overlap_record, 
    ContigType& ref_contig, int& ref_index
) {
  //cout << "**********************************" << endl;
  //cout << "   TryMergeOverlapContigs: begin" << endl;
  //cout << "   TryMergeOverlapContigs: ref_index:  " << ref_index << endl;
  ContigType target_contig = contig_holder_[overlap_record.cid];
  //cout << " TryMergeOverlapContigs: ref sequence: " << ref_contig.sequence << endl;
  //cout << " TryMergeOverlapContigs: tgt sequence: " << target_contig.sequence << endl;
  //cout << " TryMergeOverlapContigs: ref ends: " << ref_contig.q_begin << "  " << ref_contig.q_end << endl;
  //cout << " TryMergeOverlapContigs: tgt ends: " << target_contig.q_begin << "  " << target_contig.q_end << endl;
  //cout << " TryMergeOverlapContigs: overlap positions:  " << overlap_record.ref_begin << "  " << overlap_record.ref_end << "  " << overlap_record.target_begin << " " << overlap_record.target_end << endl;
  //cout << "   TryMergeOverlapContigs: overlap seq1" << endl;
  //cout << ref_contig.sequence.substr(ref_index + overlap_record.ref_begin, overlap_record.ref_end - overlap_record.ref_begin + mer_len + 1) << endl;
  //cout << "   TryMergeOverlapContigs: overlap seq2" << endl;
  //cout << target_contig.sequence.substr(overlap_record.target_begin, overlap_record.target_end - overlap_record.target_begin + mer_len + 1) << endl;
  //cout << "   TryMergeOverlapContigs: indexes:  " << overlap_record.ref_begin << "  " << overlap_record.ref_end << "  " << overlap_record.target_begin << " " << overlap_record.target_end << endl; 
  string ref_head = ref_contig.sequence.substr(0, ref_index + overlap_record.ref_begin);
  string target_head = target_contig.sequence.substr(0, overlap_record.target_begin);
  string ref_tail = ref_contig.sequence.substr(ref_index + overlap_record.ref_end + mer_len,
      ref_contig.sequence.length() - ref_index - overlap_record.ref_end - mer_len
  );
  string target_tail = target_contig.sequence.substr(overlap_record.target_end + mer_len,
      target_contig.sequence.length() - overlap_record.target_end - mer_len
  );
  //cout << "   TryMergeOverlapContigs: ref_head: " << ref_head << endl;
  //cout << "   TryMergeOverlapContigs: target_head: " << target_head << endl;
  //cout << "   TryMergeOverlapContigs: ref_tail: " << ref_tail << endl;
  //cout << "   TryMergeOverlapContigs: target_tail: " << target_tail << endl; 
  string consensus_head, consensus_tail;
  // if consensus can be constructed
  if(!FindConsensusHead(ref_head, target_head, consensus_head) || 
     !FindConsensusTail(ref_tail, target_tail, consensus_tail))  {
    return false;
  }
  //cout << "   TryMergeOverlapContigs: consensus_head:  " << consensus_head << endl;
  //cout << "   TryMergeOverlapContigs: consensus_tail:  " << consensus_tail << endl;
  // if yes, update contig
  ContigType merged_contig;
  merged_contig.sequence = consensus_head + target_contig.sequence.substr(
         overlap_record.target_begin, overlap_record.target_end - overlap_record.target_begin + mer_len
      ) + consensus_tail;
  //cout << " TryMergeOverlapContigs: sequence: " << merged_contig.sequence << endl;    
  //merged_contig.q_begin = ref_contig.q_begin;
  //merged_contig.q_end = ref_contig.q_end;
  if(target_head.length() > ref_head.length())  {
    // update the query begin and reference index
    ref_index += target_head.length() - ref_head.length();
    //if(target_contig.q_begin < merged_contig.q_begin)  {
    //  merged_contig.q_begin = target_contig.q_begin;
    //}
  }
  if(target_tail.length() > ref_tail.length())  {
    if(target_contig.q_end > merged_contig.q_end)  {
      merged_contig.q_end = target_contig.q_end;
    }
  }
  // update the reference contig
  //cout << " TryMergeOverlapContigs: q ends: " << merged_contig.q_begin << " " << merged_contig.q_end << endl;
  ref_contig = merged_contig;
  return true;
}

bool ContigRefinement::FindConsensusHead(
    std::string& seq1, std::string& seq2, std::string& consensus
)  {
  //cout << "   FindConsensusHead:  begin" << endl;
  if(seq2 == "" || (seq1.length() >= seq2.length() && 
     seq1.substr(seq1.length() - seq2.length(), seq2.length()) == seq2))  {
    consensus = seq1;
    //cout << "   FindConsensusHead:  end" << endl;
    return true;
  } else if(seq1 == "" || (seq2.length() >= seq1.length() && 
     seq2.substr(seq2.length() - seq1.length(), seq1.length()) == seq1)) {
    consensus = seq2;
    //cout << "   FindConsensusHead:  end" << endl;
    return true;
  }
  //cout << "   FindConsensusHead:  end" << endl;
  return false;
}

bool ContigRefinement::FindConsensusTail(
    std::string& seq1, std::string& seq2, std::string& consensus
)  {
  //cout << "   FindConsensusTail:  begin" << endl;
  if(seq2 == "" || (seq1.length() >= seq2.length() && 
     seq1.substr(0, seq2.length()) == seq2))  {
    consensus = seq1;
    //cout << "   FindConsensusTail:  end" << endl;
    return true;
  } else if(seq1 == "" || (seq2.length() >= seq1.length() && 
     seq2.substr(0, seq1.length()) == seq1)) {
    consensus = seq2;
    //cout << "   FindConsensusTail:  end" << endl;
    return true;
  }
  //cout << "   FindConsensusTail:  end" << endl;
  return false;
  
}

std::string ContigRefinement::FixedWidthString(int len, int num)  {
  string num_s = to_string(num); 
  if(len <= 0 || len < (int) num_s.length())  {
    cout << "Error: ContigRefinement::FixedWidthString: length is too short to include all digits" << endl;
    cout << "max_num_digits:  " << len << ",  output: " << num_s << endl;
    //exit(0);
  } else  {
    num_s += string(len - num_s.length(), ' ');
  }
  return num_s;
}
