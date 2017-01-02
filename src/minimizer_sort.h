#ifndef _MINIMIZER_SORT_H_
#define _MINIMIZER_SORT_H_

#include <iostream>
#include <unordered_map>
#include <string>
#include <cstring>
#include <list>
#include <cassert>

#include "kmer_unitcoder.h"

class MinimizerSort {
 public:
  explicit MinimizerSort() {}
  ~MinimizerSort() {}
  
  void SortSeqs(
      KmerUnitcoder &encoder, const int hash_size, 
      const int num_seqs, char **header, char **seq, int *order
  );
  void SortSeqs(
      KmerUnitcoder &encoder, const int hash_size, 
      const int num_seqs, char **seq, int *order
  );

 protected:
  // recorder the sequences based on the order given
  // in the reordered sequences, the sequence that was originally the order[i]th will be put in the ith
  void ReorderSeqs(const int num_seqs, int *order, char **header, char **seq);
  void ReorderSeqs(const int num_seqs, int *order, char **seq);
};

#endif
