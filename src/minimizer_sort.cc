#include "minimizer_sort.h"

using namespace std;

void MinimizerSort::SortSeqs(
    KmerUnitcoder &encoder, const int hash_size, 
    const int num_seqs, char **header, char **seq, int *order
)  { 
  assert(hash_size > 0);
  int i, j;
  int *count_entries = new int [hash_size];
  for(i = 0; i < hash_size; ++ i) count_entries[i] = 0; 
  int mer_len = encoder.GetMerLen();
  // counting num. entries for each hash bucket
  for(i = 0; i < num_seqs; ++ i) {
    KmerUnitType min_en = encoder.Encode(seq[i]);
    KmerUnitType en = min_en;
    for(j = mer_len; j < strlen(seq[i]); ++ j) {
      en = encoder.RightExt(en, seq[i][j]);
      if(en < min_en) min_en = en;
    }
    int pos = min_en % hash_size;
    ++ count_entries[pos];
  }
  // allocate memory
  int **encode_hash = new int* [hash_size];
  for(i = 0; i < hash_size; ++ i) {
    if(count_entries[i] > 0)  {  
      encode_hash[i] = new int [count_entries[i]];
      count_entries[i] = 0;
    }
  }
  // generate the reordering array
  for(i = 0; i < num_seqs; ++ i) {
    KmerUnitType min_en = encoder.Encode(seq[i]);
    KmerUnitType en = min_en;
    for(j = mer_len; j < strlen(seq[i]); ++ j) {
      en = encoder.RightExt(en, seq[i][j]);
      if(en < min_en) min_en = en;
    }
    int pos = min_en % hash_size;
    encode_hash[pos][count_entries[pos] ++] = i;
  }
  // copy the sequences
  int nn = 0;
  for(i = 0; i < hash_size; ++ i) {
    for(j = 0; j < count_entries[i]; ++ j) {
      order[nn ++] = encode_hash[i][j];
    }
  }

  ReorderSeqs(num_seqs, order, header, seq);
  // release the memory  
  for(i = 0; i < hash_size; ++ i) {
    if(count_entries[i] > 0)  delete [] encode_hash[i];  
  }
  delete [] encode_hash;
  delete [] count_entries;
  return;
}

void MinimizerSort::SortSeqs(
    KmerUnitcoder &encoder, const int hash_size, 
    const int num_seqs, char **seq, int *order
)  { 
  assert(hash_size > 0);
  int i, j;
  int *count_entries = new int [hash_size];
  for(i = 0; i < hash_size; ++ i) count_entries[i] = 0; 
  int mer_len = encoder.GetMerLen();
  // counting num. entries for each hash bucket
  for(i = 0; i < num_seqs; ++ i) {
    KmerUnitType min_en = encoder.Encode(seq[i]);
    KmerUnitType en = min_en;
    for(j = mer_len; j < strlen(seq[i]); ++ j) {
      en = encoder.RightExt(en, seq[i][j]);
      if(en < min_en) min_en = en;
    }
    int pos = min_en % hash_size;
    ++ count_entries[pos];
  }
  // allocate memory
  int **encode_hash = new int* [hash_size];
  for(i = 0; i < hash_size; ++ i) {
    if(count_entries[i] > 0)  {  
      encode_hash[i] = new int [count_entries[i]];
      count_entries[i] = 0;
    }
  }
  // generate the reordering array
  for(i = 0; i < num_seqs; ++ i) {
    KmerUnitType min_en = encoder.Encode(seq[i]);
    KmerUnitType en = min_en;
    for(j = mer_len; j < strlen(seq[i]); ++ j) {
      en = encoder.RightExt(en, seq[i][j]);
      if(en < min_en) min_en = en;
    }
    int pos = min_en % hash_size;
    encode_hash[pos][count_entries[pos] ++] = i;
  }
  // copy the sequences
  int nn = 0;
  for(i = 0; i < hash_size; ++ i) {
    for(j = 0; j < count_entries[i]; ++ j) {
      order[nn ++] = encode_hash[i][j];
    }
  }

  ReorderSeqs(num_seqs, order, seq);
  // release the memory  
  for(i = 0; i < hash_size; ++ i) {
    if(count_entries[i] > 0)  delete [] encode_hash[i];  
  }
  delete [] encode_hash;
  delete [] count_entries;
  return;
}

void MinimizerSort::ReorderSeqs(const int num_seqs, int *order, char **header, char **seq) {
  char **header_holder = new char* [num_seqs];
  char **seq_holder = new char* [num_seqs];
  for(int i = 0; i < num_seqs; ++ i) {
    header_holder[i] = header[i];
    seq_holder[i] = seq[i];
  }
  for(int i = 0; i < num_seqs; ++ i) {
    header[i] = header_holder[order[i]];
    seq[i] = seq_holder[order[i]];
  }
  delete [] header_holder;
  delete [] seq_holder;
  return;
}

void MinimizerSort::ReorderSeqs(const int num_seqs, int *order, char **seq) {
  char **seq_holder = new char* [num_seqs];
  for(int i = 0; i < num_seqs; ++ i) {
    seq_holder[i] = seq[i];
  }
  for(int i = 0; i < num_seqs; ++ i) {
    seq[i] = seq_holder[order[i]];
  }
  delete [] seq_holder;
  return;
}
