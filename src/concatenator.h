#ifndef _CONCATENATOR_H_
#define _CONCATENATOR_H_

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cassert>
#include <list>

// the deliminator must be lexicographically the smallest
static char DELIM = (char) 1;

class Concatenator  {
 public:
  explicit Concatenator(char** const seq, const int n, std::string &concat_seq);
  explicit Concatenator(std::list<std::string> &seq, std::string &concat_seq);

  // the following functions break the reads into blocks
  // each block is represented by its concatenated sequence and the id begin
  // that is the kth block contain the id_begin[k]...id_begin[k+1]-1 th reads
  explicit Concatenator(
      char** const seq, const int n, const long long int max_size, 
      std::vector<std::string> &concat_seqs, std::vector<int> &id_begin
  );
  explicit Concatenator(
      std::list<std::string> &seq, const long long int max_size, 
      std::vector<std::string> &concat_seqs, std::vector<int> &id_begin
  );
  ~Concatenator() {}
};

#endif
