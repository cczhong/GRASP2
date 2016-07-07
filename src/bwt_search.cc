#include "bwt_search.h"

using namespace std;

// sorting the AlignType
bool SortAlignTypeBySize(const AlignType &a, const AlignType &b)  {
  int a_dist = a.q_end - a.q_begin;
  int b_dist = b.q_end - b.q_begin;
  if(a_dist > b_dist || (a_dist == b_dist && a.q_begin < b.q_begin))  {
    return true;
  }  
  return false;
} 

void BWTSearch::SearchExact(BWT &bwt, const char *str, AlignType &pos) {
  
  int len = strlen(str);
  if(len <= 0)  {
    cerr << "Searching empty sequence. An empty list will be returned." << endl;
    return;
  }
  pair<BWTIDX, BWTIDX> range;
  range.first = 0, range.second = bwt.GetSize();
  int i;
  for(i = strlen(str) - 1; i >= 0; -- i) {
    // indicating that there is no solution
    if(range.first >= range.second)  break;
    range = bwt.UpdateRange(str[i], range);    
  }
  // record the positions
  pos.q_begin = i + 1; pos.q_end = strlen(str) - 1; pos.cost = 0;
  pos.bwt_begin = range.first; pos.bwt_end = range.second;
  return;
}

void BWTSearch::SearchAllSubRegions(
    BWT &bwt, const int min_len, const char *str, std::vector<AlignType> &all_pos
) {
  int len = strlen(str);
  if(len <= 0)  {
    cerr << "BWTSearch::Search: Searching empty sequence. An empty list will be returned." << endl; return;
  }
  BWTIDX bwt_size = bwt.GetSize();
  // reverse the string
  string s = str; s = string(s.rbegin(), s.rend());
  //cout << "Query sequence:  " << str << endl;
  // construct a vector of ranges
  vector<AlignType> cand_pos; cand_pos.resize(len);
  vector<AlignType> temp_pos;
  for(int i = 0; i < cand_pos.size(); ++ i) {
    //cout << "=======================" << endl;
    cand_pos[i].bwt_begin = 0; cand_pos[i].bwt_end = bwt_size;
    cand_pos[i].q_begin = cand_pos[i].q_end = len - i - 1; 
    ++ cand_pos[i].q_begin; // this is because we need to subtract 1 for search of the first character
    for(int j = 0; j <= i; ++ j) {
      if(cand_pos[j].bwt_begin >= cand_pos[j].bwt_end) continue;
      BWTIDX x, y; 
      bwt.UpdateRange(s[i], cand_pos[j].bwt_begin, cand_pos[j].bwt_end, x, y);
      if(x >= y) { // this is the end of this search, record the range
        if(cand_pos[j].q_end - cand_pos[j].q_begin + 1 >= min_len)  {
          //cout << "hit recorded" << endl;
          temp_pos.push_back(cand_pos[j]);
        }
      }
      cand_pos[j].q_begin --; cand_pos[j].bwt_begin = x; cand_pos[j].bwt_end = y;
      //cout << "index: " << i << "  " << j  << "  " << s[i] << " query location:  " << cand_pos[j].q_begin << "  " << cand_pos[j].q_end << " " << (long long int) cand_pos[j].bwt_begin << " " << (long long int) cand_pos[j].bwt_end << endl;
    }
  }
  // record the final identified ranges
  for(int i = 0; i < cand_pos.size(); ++ i) {
    if(cand_pos[i].bwt_begin < cand_pos[i].bwt_end && cand_pos[i].q_end - cand_pos[i].q_begin + 1 >= min_len)
      temp_pos.push_back(cand_pos[i]);
  }
  //cout << "size of fragments: " << temp_pos.size() << endl;
  if(temp_pos.size() <= 0)  return;
  // remove duplicates (such that no one range is a subrange of another)
  sort(temp_pos.begin(), temp_pos.end(), SortAlignTypeBySize);
  for(int i = 0; i < temp_pos.size() - 1; ++ i) {
    if(temp_pos[i].q_begin < 0) continue;
    for(int j = i + 1; j < temp_pos.size(); ++ j) {
      // if the ith range contains the jth range, invalidate the jth
      if(temp_pos[i].q_begin <= temp_pos[j].q_begin && temp_pos[i].q_end >= temp_pos[j].q_end)  
        temp_pos[j].q_begin = -1;
    }
  }
  // record all remaining non-redundant set
  for(int i = 0; i < temp_pos.size(); ++ i) {
    if(temp_pos[i].q_begin >= 0)
      all_pos.push_back(temp_pos[i]);
  }
  return;
  /* OBSOLETE
  BWTIDX bwt_size = bwt.GetSize();
  BWTIDX *begin = new BWTIDX [len + 1]; begin[len] = 0;
  BWTIDX *end = new BWTIDX [len + 1]; end[len] = bwt_size;
  int i, j, k;
  // extend the searches one-by-one
  int right_boundary = len;
  for(i = len - 1; i >= 0; -- i) {
    bwt.UpdateRange(str[i], begin[i + 1], end[i + 1], begin[i], end[i]);   
    //cout <<  i << ":  " << begin[i + 1] << "  " << end[i + 1] << "  " << begin[i] << "  " << end[i] << endl; 
    if(begin[i] >= end[i])  {
      // if this is the terminal, record the information 
      if(right_boundary - i >= min_len) {
        AlignType a; 
        a.q_begin = i + 1; a.q_end = right_boundary > len - 1 ? len - 1 : right_boundary, 
        a.bwt_begin = begin[i + 1]; a.bwt_end = end[i + 1];
        all_pos.push_back(a);
      }
      // find the next right bound
      while(i >= 0) {
        bwt.UpdateRange(str[i], 0, bwt_size, begin[i], end[i]);
        if(begin[i] < end[i])  {
          right_boundary = i; break;
        } else  {
          // NOTE!!!: the change of index "i" will affect the outer loop
          -- i;
        }
      }
    }
  }
  // record the last chunk
  if(begin[0] < end[0] && right_boundary >= min_len) {
    AlignType a; 
    a.q_begin = 0; a.q_end = right_boundary > len - 1 ? len - 1 : right_boundary; 
    a.bwt_begin = begin[0]; a.bwt_end = end[0];
    all_pos.push_back(a);
  }
  delete [] begin;
  delete [] end;
  return;
  */  
}

void BWTSearch::Search(
    BWT &bwt, BWT &rev_bwt, 
    const char *str, AlignType &pos
) {
  int len = strlen(str);
  if(len <= 0)  {
    cerr << "BWTSearch::Search: Searching empty sequence. An empty list will be returned." << endl; return;
  }
  int *bound = new int [len];
  CalLowerBound(rev_bwt, str, bound);
  SearchInExact(bwt, rev_bwt, bound, str, pos);  
  AlignType prev_pos = pos;
  if(pos.q_begin >= 0)  {
      string rev_seq = str;
      rev_seq = string(rev_seq.rbegin(), rev_seq.rend());
      CalLowerBound(bwt, rev_seq.c_str(), bound);   
      SearchInExact(rev_bwt, bwt, bound, rev_seq.c_str(), pos);
  }
  int prev_span = prev_pos.q_end - prev_pos.q_begin + 1;
  int span = pos.q_end - pos.q_begin + 1;
  if(prev_span > span || (prev_span == span && prev_pos.cost < pos.cost)) { 
    pos = prev_pos;
  } else  {
    // recall that "pos" was the search result of the reverse BWT
    // need to flip the starting and ending positions for the query
    pos.q_end = len - pos.q_begin - 1; pos.q_begin = 0;
    // TODO: notice that the BWT indexes are for the reverse BWT
    // TODO: and should be converted to exact reference location before returning
  }
  delete [] bound;
  return;
}


void BWTSearch::Enqueue(
    std::priority_queue<ExtInfo, std::vector<ExtInfo>, ExtInfoComp> &candidate, 
    ExtInfo &phase_info
) {
  if(candidate.size() < max_queue_size || phase_info.cost <= candidate.top().cost)
    candidate.push(phase_info);
  return;
}

void BWTSearch::SearchInExact(
    BWT &bwt, BWT &rev_bwt, const int *bound,
    const char *str, AlignType &pos
) {
  
  int len = strlen(str);
  int lower_cost = cost;
  // construct initial stack 
  priority_queue<ExtInfo, std::vector<ExtInfo>, ExtInfoComp> candidate;
  ExtInfo init;
  init.range.first = 0; init.range.second = bwt.GetSize(); init.q_pos = len - 1; init.cost = 0;
  init.cost_bound = bound[len - 1];
  Enqueue(candidate, init);
  pos.bwt_begin = 0; pos.bwt_end = bwt.GetSize(); 
  pos.q_begin = pos.q_end = len - 1; pos.cost = cost;
  // extend each element in the stack with fast mode
  while(!candidate.empty()) {
    // record the best extension so par
    if(candidate.top().q_pos < pos.q_begin || 
        (candidate.top().q_pos == pos.q_begin && candidate.top().cost < pos.cost)
    ) { 
      pos.q_begin = candidate.top().q_pos + 1; pos.cost = candidate.top().cost;
      pos.bwt_begin = candidate.top().range.first; pos.bwt_end = candidate.top().range.second;
    }
       
    if(candidate.top().q_pos < 0)  {
      if(candidate.top().cost < lower_cost) lower_cost = candidate.top().cost;
      candidate.pop();
    } else 
      ExtendInExactFast(bwt, str, candidate, bound, lower_cost);
  }
  // if the entire sequence cannot be aligned, report the longest alignment within 
  // the allowed number of errors
  return;
}

void BWTSearch::ExtendInExactFast(
    BWT &bwt, const char *str,
    std::priority_queue<ExtInfo, std::vector<ExtInfo>, ExtInfoComp> &candidate, 
    const int *bound, int lower_cost
) {

  if(candidate.empty()) return;
  // get the top extension
  ExtInfo current = candidate.top();
  candidate.pop(); 
  if(current.q_pos < 0 || current.cost + bound[current.q_pos] > lower_cost)  return;
  // perform extension
  ExtInfo next;
  // tries to find if the next is a match
  next.range = bwt.UpdateRange(str[current.q_pos], current.range);
  if(next.range.first < next.range.second)  {
    next.q_pos = current.q_pos - 1; next.cost = current.cost;
    next.cost_bound = next.cost;
    if(next.q_pos >= 0) next.cost_bound += bound[next.q_pos];
    if(next.q_pos > 0)  {
      pair<BWTIDX, BWTIDX> peek_ahead = bwt.UpdateRange(str[current.q_pos - 1], next.range);
      if(peek_ahead.first < peek_ahead.second && next.cost_bound <= lower_cost) {
        Enqueue(candidate, next); return;
      }
    } else {
      // can't peek ahead, directly insert
      if(next.cost_bound <= lower_cost) {
        Enqueue(candidate, next); return;
      }
    }
  }
  // otherwise try all possible mismatch/insertion/deletions
  for(int i = 0; i < bwt.alphabet_.alphabet_size_; ++ i) {
    BWTCHAR c = bwt.alphabet_.GetInvCharMap(i);
    // refine the range
    next.range = bwt.UpdateRange(c, current.range);
    if(next.range.first >= next.range.second) continue;    
    // the match/mismatch case; subtract the index
    if(c == str[current.q_pos]) next.cost = current.cost;
    else next.cost = current.cost + m_cost;
    next.q_pos = current.q_pos - 1;
    next.cost_bound = next.cost;
    if(next.q_pos >= 0) next.cost_bound += bound[next.q_pos];
    if(next.cost_bound <= lower_cost) {
      Enqueue(candidate, next);
    }   
    // the insertion case (any character "c" is inserted into the sequence)
    // add gap cost, keep index
    next.cost = current.cost + g_cost;
    next.q_pos = current.q_pos;
    next.cost_bound = next.cost;
    if(next.q_pos >= 0) next.cost_bound += bound[next.q_pos];
    if(next.cost_bound <= lower_cost) {
      Enqueue(candidate, next);
    }
  }
  // the deletion case (skip current character, add gap cost)
  next = current;
  -- next.q_pos; next.cost += g_cost;
  next.cost_bound = next.cost;
  if(next.q_pos >= 0) next.cost_bound += bound[next.q_pos];
  if(next.cost_bound <= lower_cost) {
    Enqueue(candidate, next);
  }
  return;
}

void BWTSearch::CalLowerBound(BWT &rev_bwt, const char *str, int *bound) {
  int cost = 0;
  pair<BWTIDX, BWTIDX> range;
  range.first = 0; range.second = rev_bwt.GetSize();
  for(int i = 0; i < strlen(str); ++ i) {
    range = rev_bwt.UpdateRange(str[i], range);
    if(range.first >= range.second)  {
      range.first = 0; range.second = rev_bwt.GetSize();
      cost += m_cost;
    }
    bound[i] = cost;
  }
  return;
}

// given a read of interest and minimum overlap, find all intervals (in both fw and re FM-indexes) 
// corresponding to the prefix of the reads that perfectly overlap with the given read
void BWTSearch::SearchBeginIntervals(const char* seq, const int min_len, IvInfo &search_info) {
  int n = strlen(seq);
  // if the total length of the read is less than the minimum overlap, return
  if(n < min_len) return;
  // initialize the search for the overlapping region
  pair<BWTIDX, BWTIDX> fw_range, re_range, fw_range_terminal, re_range_terminal;
  fw_range.first = re_range.first = 0; 
  fw_range.second = re_range.second = search_info.bwtF_->GetSize();
  int i;
  for(i = 0; i < min_len - 1; ++ i) {
    //cout << "search sequence: " << &seq[n - i - 1] << endl;
    char c = seq[n - i - 1];
    BWTIDX occbegin = search_info.bwtF_->CountOccurrence(c, fw_range.first);
    BWTIDX occend = search_info.bwtF_->CountOccurrence(c, fw_range.second);
    //cout << "initial range: " << fw_range.first << "  " << fw_range.second << "  " << re_range.first << "  " << re_range.second << endl;
    //cout << "search occurrence:  " << c << " " << occbegin << "  " << occend << endl;
    if(occbegin >= occend)  break;
    //cout << "lexicographical smaller: " << search_info.bwtF_->CountLexicoLess(c, fw_range.first) << "  " << search_info.bwtF_->CountLexicoLess(c, fw_range.second) << endl;
    re_range.first = re_range.first
        + search_info.bwtF_->CountLexicoLess(c, fw_range.second)
        - search_info.bwtF_->CountLexicoLess(c, fw_range.first); 
    re_range.second = re_range.first + occend - occbegin;
    int c_id = search_info.bwtF_->alphabet_.GetCharMap(c);
    fw_range.first = search_info.bwtF_->acc_freq_[c_id + 1] + occbegin;
    fw_range.second = search_info.bwtF_->acc_freq_[c_id + 1] + occend;
    //cout << "updated range: " << fw_range.first << "  " << fw_range.second << "  " << re_range.first << "  " << re_range.second << endl;
  }
  // if no sequence overlap for the given minumum overlap length, return
  if(i < min_len - 1) return;
  // otherwise also try to search the delimitor to detect begin intervals
  // recall that i is less than n - 1 because we do not want the read itself
  // or other wise the read is contained by other reads
  for(i = min_len - 1; i < n - 1; ++ i) {
    // extend the sequence
    //cout << "search sequence: " << &seq[n - i - 1] << " " << n << " " << i << endl;
    //cout << "initial range: " << fw_range.first << "  " << fw_range.second << "  " << re_range.first << "  " << re_range.second << endl;
    char c = seq[n - i - 1];
    BWTIDX occbegin = search_info.bwtF_->CountOccurrence(c, fw_range.first);
    BWTIDX occend = search_info.bwtF_->CountOccurrence(c, fw_range.second);
    //cout << "search occurrence:  " << c << " " << occbegin << "  " << occend << endl;
    //cout << "lexicographical smaller: " << search_info.bwtF_->CountLexicoLess(c, fw_range.first) << "  " << search_info.bwtF_->CountLexicoLess(c, fw_range.second) << endl;
    if(occbegin >= occend)  break;
    re_range.first = re_range.first
        + search_info.bwtF_->CountLexicoLess(c, fw_range.second)
        - search_info.bwtF_->CountLexicoLess(c, fw_range.first);
    re_range.second = re_range.first + occend - occbegin;
    int c_id = search_info.bwtF_->alphabet_.GetCharMap(c);
    fw_range.first = search_info.bwtF_->acc_freq_[c_id + 1] + occbegin;
    fw_range.second = search_info.bwtF_->acc_freq_[c_id + 1] + occend;
    // also search for the delimitor (with the updated fw_range and re_range)
    //cout << "initial terminal range: " << fw_range.first << "  " << fw_range.second << "  " << re_range.first << "  " << re_range.second << endl;
    occbegin = search_info.bwtF_->CountOccurrence(DELIM, fw_range.first);
    occend = search_info.bwtF_->CountOccurrence(DELIM, fw_range.second);
    //cout << "search terminal occurrence:  " << occbegin << "  " << occend << endl;
    if(occbegin >= occend)  continue;
    // note that delimitor is the lexico-smallest char, 
    // no need to add CountLexicoLess nor acc_freq
    re_range_terminal.first = re_range.first; 
    re_range_terminal.second = re_range.first + occend - occbegin;
    fw_range_terminal.first = occbegin;
    fw_range_terminal.second = occend;
    // record such interval
    search_info.intervals_.PushA(fw_range_terminal.first, fw_range_terminal.second);
    search_info.intervals_.PushB(re_range_terminal.first, re_range_terminal.second);
    search_info.intervals_.PushLen(i + 1);
    
    //cout << "overlap recorded !!!" << endl;
  }
  return;
}


void BWTSearch::FindIrreducible(
    IvInfo &search_info, std::vector<BWTIDX> &ir_positions, std::vector<int> &ir_overlap
) {
  
  search_info.intervals_.Reverse();  
  // check boundary conditions
  if(!search_info.intervals_.Check()) {
    cout << "Warning: corrupted intervals, no irreducible edges can be detected." << endl;
    return;
  }
  if(search_info.intervals_.GetSize() <= 0)  return;
  // the stack contains all intervals that ends with the same sequences
  stack<IvSet> candidates;
  candidates.push(search_info.intervals_);
  // recursively check each intervals
  bool used_char[256];  
  while(!candidates.empty()) {
    //cout << "============ handling each interval group ===============  " << candidates.size() << endl;
    IvSet current = candidates.top(); candidates.pop();
    //for(int h = 0; h < current.GetSize(); ++ h) {
    //  cout << "search_info " << h << ": " << current.len_[h] << "  " << current.ivA_[h].first << " " << current.ivA_[h].second << "  " << current.ivB_[h].first << " " << current.ivB_[h].second << endl;
    //}
    // check all presented characters in the alphabet
    // note that for the first time we do not check "$"-extension because it would
    // mean that the extension read is contained
    int i, j, k, n = current.GetSize();
    memset(used_char, 0, 256);
    for(i = 0; i < n; ++ i) { // for each interval
      // for each character in the alphabet (look at the reverse BWT)
      for(k = current.ivB_[i].first; k < current.ivB_[i].second; ++ k) { 
        //cout << "Interval info: " << k << ": " <<  current.len_[i] << "  " << current.ivA_[i].first << " " << current.ivA_[i].second << "  " << current.ivB_[i].first << " " << current.ivB_[i].second << endl;
        char c = (char) search_info.bwtR_->bwt_[k]; // the kth char in the reverse BWT string
        // check if the character has been tested
        if(c == DELIM || used_char[c])  continue;
        used_char[c] = true;
        IvSet next;
        // try to update the intervals by appending such character
        for(j = i; j < n; ++ j) {
          // try to append the character
          pair<BWTIDX, BWTIDX> r1 = search_info.bwtR_->UpdateRange(c, current.ivB_[j]);
          // try to check if the read ends after appending the character
          pair<BWTIDX, BWTIDX> r2 = search_info.bwtR_->UpdateRange(DELIM, r1); 
          //cout << "phase: " << j << ": " << r1.first << " " << r1.second << " " << r2.first << "  " << r2.second << endl; 
          // if the first read leads to a termination, we found an irreducible read
          // in cases where multiple reads end at the same time, take the first position
          // terminate current loop
          if(j == i && r1.second - r1.first == r2.second - r2.first)  {
            //cout << "irreducible read found!!!" << endl;
            ir_positions.push_back(r2.first); 
            ir_overlap.push_back(current.len_[j]);
            break;
          }
          // for other reads that do not terminate, add to interval set next
          if(r1.second - r1.first >= 1 && r2.second - r2.first <= 0)  {
            // we only need to record intervals at the reverse BWT
            //cout << "read extension recorded!!!" << endl;
            next.PushA(-1, -1); next.PushB(r1.first, r1.second); next.PushLen(current.len_[j]);
          }
        }
        if(next.GetSize() > 0)  candidates.push(next);
      }
    }
  }
  return;
}

bool BWTSearch::IsContainedRead(const char* seq, BWT &bwt, AlignType &pos)  {
  // constructing the extended string
  int n = strlen(seq);
  char *extended_seq = new char [n + 3];
  extended_seq[0] = DELIM;
  strcpy(&extended_seq[1], seq);
  extended_seq[n + 1] = DELIM; extended_seq[n + 2] = '\0';
  // searching the entire sequence with and without delimitor against the BWT
  SearchExact(bwt, seq, pos);
  int r1 = pos.bwt_end - pos.bwt_begin;
  SearchExact(bwt, extended_seq, pos);
  int r2 = pos.bwt_end - pos.bwt_begin;
  delete [] extended_seq;
  if(r2 == r1)  return false;
  return true;
}

