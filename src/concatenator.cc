#include "concatenator.h"

using namespace std;

Concatenator::Concatenator(char** const seq, const int n, std::string &concat_seq)  {
  if(n <= 0)  {
    cerr << "Concatenator::FatalError: Empty input sequences; abort." << endl;
    exit(1);
  }  
  // concatenate the sequences
  concat_seq = "";
  concat_seq += DELIM;
  for(int i = 0; i < n; ++ i) {
    concat_seq += seq[i];
    concat_seq += DELIM;
  }
  return;
}

Concatenator::Concatenator(std::list<std::string> &seq, std::string &concat_seq)  {
  if(seq.size() <= 0)  {
    cerr << "Concatenator::FatalError: Empty input sequences; abort." << endl;
    exit(1);
  }  
  concat_seq = "";
  concat_seq += DELIM;
  // concatenate the sequences
  for(auto it = seq.begin(); it != seq.end(); ++ it) {
    concat_seq += *it;
    concat_seq += DELIM;
  }
  return;
}

Concatenator::Concatenator(
    char** const seq, const int n, const long long int max_size, 
    std::vector<std::string> &concat_seqs, std::vector<int> &id_begin
)  {
  if(n <= 0)  {
    cerr << "Concatenator::FatalError: Empty input sequences; abort." << endl;
    exit(1);
  }    
  long long int current_size = 1;
  string current_seq = ""; current_seq += DELIM;
  id_begin.push_back(0);
  for(int i = 0; i < n; ++ i) {
    if(current_size + strlen(seq[i]) + 1 >= max_size)  {
      if(!current_seq.empty())  {
        concat_seqs.push_back(current_seq); id_begin.push_back(i);
        current_seq = ""; current_seq += DELIM; current_size = 1;
      }
    }
    current_seq += seq[i];
    current_seq += DELIM;
    current_size += strlen(seq[i]) + 1;
  }
  concat_seqs.push_back(current_seq);  
  return;
}

Concatenator::Concatenator(
    std::list<std::string> &seq, const long long int max_size, 
    std::vector<std::string> &concat_seqs, std::vector<int> &id_begin
)  {
  if(seq.size() <= 0)  {
    cerr << "Concatenator::FatalError: Empty input sequences; abort." << endl;
    exit(1);
  }    
  long long int current_size = 1;
  string current_seq = ""; current_seq += DELIM;
  id_begin.push_back(0);
  int idx = 0;
  for(auto it = seq.begin(); it != seq.end(); ++ it) {
    if(current_size + it->length() + 1 >= max_size)  {
      if(!current_seq.empty())  {
        concat_seqs.push_back(current_seq); id_begin.push_back(idx);
        current_seq = ""; current_seq += DELIM; current_size = 1;
      }
    }
    current_seq += *it;
    current_seq += DELIM;
    current_size += it->length() + 1;
    ++ idx;
  }
  concat_seqs.push_back(current_seq);
  return;
}
