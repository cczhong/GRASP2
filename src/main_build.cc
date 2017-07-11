#include "align_batch.h"
#include "bwt.h"
#include "bwt_search.h"
#include "loader.h"
#include "bio_alphabet.h"
#include "kmer_unitcoder.h"
#include "minimizer_sort.h"
#include "string_graph.h"
#include "sequence_search.h"
#include "kmer_unitcoder.h"
#include "scoring_prot.h"
#include "reduced_alphabet.h"
#include "kmer_filtering.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <ctime>
#include <cmath>

using namespace std;

static string workspace_dir;
static string db_file;
static int extd_len;
static int num_threads = 1;
static int neighbor_score = 11;
static string verbose;
static int scoring_matrix = 0;
static int mer_len = 3;

void PrintUsage()  {
  cout << "Usage: grasp-build [peptide_db (FASTA)]" << endl;
  cout << "Use \'--help\' for more options" << endl;
  return;
}

double MyTime (void)
{
    int flag;
    clockid_t cid = CLOCK_REALTIME; // CLOCK_MONOTONE might be better
    timespec tp;
    double timing;
	
    flag = clock_gettime(cid, &tp);
    if (flag == 0) timing = tp.tv_sec + 1.0e-9*tp.tv_nsec;
    else           timing = -17.0;         // If timer failed, return non-valid time
	
    return(timing);
}

void PrintElapsed( double s, double e, const char *task )
{
	double elapsed = e - s ;
	printf ("[%02.0f:%02.0f:%05.2f]\t%s\n", 
			floor(elapsed/3600.0), 
			floor(fmod(elapsed,3600.0)/60.0), 
			fmod(elapsed,60.0),
			task);
	return;
}

string GetFileStem(const string& path)  {
  // going backward untill the '\/' character
  int i;
  for(i = path.length() - 1; i >= 0; -- i) {
    if(path[i] == '/')  break;
  }
  return path.substr(i + 1, path.length() - i - 1);
}

int main(int argc, char** argv)  {
  // reading options
  boost::program_options::options_description desc("List of options");
  desc.add_options()
      ("help", "print the help message")
      ("db_file", boost::program_options::value<string>(&db_file), "short-peptide reads (in FASTA format)")
      ("index", boost::program_options::value<string>(&workspace_dir)->default_value("index"), "working directory for indexing file dump")
      ("extension_len", boost::program_options::value<int>(&extd_len)->default_value(10), "minimum overlap length for path extension")
      ("neighbor_score", boost::program_options::value<int>(&neighbor_score)->default_value(11), "neighbor score for 3-mer seed matches") 
      ("num_threads", boost::program_options::value<int>(&num_threads)->default_value(1), "maximum number of threads to be used")
      ("verbose", boost::program_options::value<string>(&verbose), "print intermediate information (default true)")
  ;
  boost::program_options::positional_options_description pos_opt;
  pos_opt.add("db_file", 1);
  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv).
      options(desc).positional(pos_opt).run(), vm
  );
  boost::program_options::notify(vm);
  if(vm.count("help"))  {
    cout << "Usage: grasp-build [peptide_db (FASTA)]" << endl << endl;
    cout << desc << endl; 
    return 0;
  }
  // check options validity
  boost::filesystem::path abs_workspace = workspace_dir;
  if(!boost::filesystem::exists(db_file))  {
    cout << db_file << endl;
    cout << "Error: grasp-build: db_file does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(!boost::filesystem::is_directory(workspace_dir))  {
    cout << workspace_dir << endl;
    cout << "Error: grasp-build: working space does not exist (please provide full path)." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(extd_len < 6 || extd_len > 20)  {
    cout << "Error: grasp-build: extension length out of range (allowed range: 6-20)." << endl;
    exit(0);
  }
  bool is_verbose = true;
  if(verbose == "False" || verbose == "false" || verbose == "No" || verbose == "no" || verbose == "0")  {
    is_verbose = false;
  }
  
  
  if(is_verbose)  {
    cout << "============================================================" << endl;
    cout << "GRASP2-Build: Begin of program execution." << endl;
  }
  
  BioAlphabet protein_alphabet(PROT);
  ReducedAlphabet reduced_alphabet((enum Alphabet) 10);
  ScoringProt scoring_function(static_cast<enum MatrixName>(scoring_matrix), -10, -1); 
  
  // Load in the peptide sequences to be searched against
  double start_time = MyTime();
  double check_time;
  Loader pepdb_loader;
  int num_seqs = pepdb_loader.CountFastaNumSeqs(db_file.c_str());
  char **seqs = new char* [num_seqs];
  //num_seqs = pepdb_loader.LoadFasta(protein_alphabet, db_file.c_str(), header, seqs);
  pepdb_loader.LoadFasta(protein_alphabet, db_file.c_str(), seqs);
  string concat_seq; 
  //Concatenator concat_obj(seqs, num_seqs, concat_seq);
  
  //cout << "check seq: " << seqs[823927] << endl;
  //cout << "check seq: " << seqs[7792734] << endl;
  //return 0;

  // sort the reads based on minimizers
  KmerUnitcoder min_sort_mer(protein_alphabet, 6);
  MinimizerSort m_sort;
  int *order = new int [num_seqs];
  m_sort.SortSeqs(min_sort_mer, 10000000, num_seqs, seqs, order);
  
  //for(int i = 0; i < num_seqs; ++ i) {
  //  cout << seqs[i] << endl;
  //}
  //return 0;

  vector<string> concat_seqs;
  vector<int> id_begin;
  Concatenator concat_obj(seqs, num_seqs, 100000000, concat_seqs, id_begin);
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Load peptide database done.";
    PrintElapsed(start_time, check_time, "");
  }
  start_time = MyTime();
  //for(auto it = concat_seqs.begin(); it != concat_seqs.end(); ++ it) {
  //  cout << *it << endl << endl;
  //}

  // construct multiple BWT for each block of sequences
  int seq_idx = 0;
  for(auto it = concat_seqs.begin(); it != concat_seqs.end(); ++ it) {
    string bwt_prefix = workspace_dir + "/" + "bwt." + to_string(seq_idx);
    string rev_bwt_prefix = workspace_dir + "/" + "rev_bwt." + to_string(seq_idx);
    // construct forward sequence BWT and write it on hard disk
    BWT bwt;
    bwt.Construct(protein_alphabet, it->c_str());
    bwt.WriteIndex(bwt_prefix);
    bwt.Purge();
    // construct reverse sequence BWT and write it on hard disk
    //string rev_seq = string(it->rbegin(), it->rend());
    //BWT rev_bwt;
    //rev_bwt.Construct(protein_alphabet, rev_seq.c_str());
    //rev_bwt.WriteIndex(rev_bwt_prefix);
    //rev_bwt.Purge();
    ++ seq_idx;
  }
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Construct Burrows-Wheeler transformation done.";
    PrintElapsed(start_time, check_time, "");
  }
  start_time = MyTime();

  // load BWT one-by-one and conduct the search
  seq_idx = 0;
  StringGraph strG;
  BWTSearch bwt_searcher;
  //std::vector<std::vector<TargetOverlapType> > extension;
  vector<vector<TargetOverlapType> *> *extension = new vector<vector<TargetOverlapType> *>;
  extension->resize(num_seqs);
  for(BWTINT i = 0; i < num_seqs; ++ i) {
    (*extension)[i] = new vector<TargetOverlapType>;
  }
  vector<bool> *contained = new vector<bool> (num_seqs, false);
  //TargetOverlapType **extension = new TargetOverlapType* [num_seqs];
  //int *ext_count = new int [num_seqs];
  //cout << "num_threads: " << num_threads << endl;
  for(auto it = concat_seqs.begin(); it != concat_seqs.end(); ++ it) {
    string bwt_prefix = workspace_dir + "/" + "bwt." + to_string(seq_idx);
    //string rev_bwt_prefix = workspace_dir + "/" + "rev_bwt." + to_string(seq_idx);
    // load forward and backward BWTs from index
    BWT reload_bwt;
    reload_bwt.ConstructFromIndex(protein_alphabet, it->c_str(), bwt_prefix);
    //string rev_seq = string(it->rbegin(), it->rend());
    //BWT reload_rev_bwt;
    //reload_rev_bwt.ConstructFromIndex(protein_alphabet,rev_seq.c_str(), rev_bwt_prefix);
    //cout << "Finish loading index" << endl;
    strG.MultiComputeExtension(
        num_threads, extd_len, num_seqs, seqs, 
        id_begin[seq_idx], reload_bwt, extension, contained
    );
    //cout << "Done computing multiple extension" << endl;
    reload_bwt.Purge();
    it->clear();
    ++ seq_idx;
  }
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Compute read overlap done.";
    PrintElapsed(start_time, check_time, "");
  }
  start_time = MyTime();

  // TODO: handle the read ID that are shuffled by minimizer sorting

  /*
  for(int i = 0; i < num_seqs; ++ i) {
    if(i != 6918562) continue;
    cout << "i: " << i << " " << seqs[i] << endl;    
    if((*contained)[i]) continue;
    cout << "source:  " << seqs[i] << endl;
    for(int j = 0; j < (*extension)[i]->size(); ++ j) {
      int tid = (*(*extension)[i])[j].rid;
      if((*contained)[tid])  continue;
      cout << "target:  " << seqs[tid] << endl;
    }
    cout << "=======================" << endl;
  }
  */
  vector<vector<TargetOverlapType> *> *rev_extension = new vector<vector<TargetOverlapType> *>;
  rev_extension->resize(num_seqs);
  for(BWTINT i = 0; i < num_seqs; ++ i) {
    (*rev_extension)[i] = new vector<TargetOverlapType>;
  }
  strG.FillRevExtension(contained, extension, rev_extension);

  string db_stem = GetFileStem(db_file);
  string idx_unitig_file = workspace_dir + "/" + db_stem + ".utg";  

  strG.WriteUnitigsFromExtension(seqs, contained, extension, rev_extension, order, idx_unitig_file);
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Write unitigs done";
    PrintElapsed(start_time, check_time, "");
  }

  // TODO: memory collection for extension and rev_extension
  
  //return 0;
  // Compute the length of the sequence represent by the edge connecting the two sequences
  //strG.ComputeEdgeLen(num_seqs, seqs, extension);
  // remove contained reads
  //vector<bool> contained(num_seqs, false);
  //strG.DetectContained(num_seqs, seqs, extension, contained);
  //strG.RemoveReducibleEdges(num_seqs, contained, extension);
  strG.ImportExtension(num_seqs, contained, extension);
  extension->clear(); delete extension;
  rev_extension->clear(); delete rev_extension;
  contained->clear(); delete contained;

  //return 0;
  //***************************
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Construct string graph done.";
    PrintElapsed(start_time, check_time, "");
  }
  //***************************
  // Post-processing of the string graph
  start_time = MyTime();

  
  strG.CheckGraph(); 

  /*
  while(strG.RemoveTipsBeforeCondense())  {;}
  for(int i = 3; i < 10; i += 3) {
    int num_removed_right = strG.RemoveBubbleRight(i);
    while(strG.RemoveTipsBeforeCondense())  {;}
    int num_removed_left = strG.RemoveBubbleLeft(i);
    while(strG.RemoveTipsBeforeCondense())  {;}
    //cout << "Num edges removed: " << num_removed_right << " " << num_removed_left << endl;
  }
  */

  strG.CondenseGraph(seqs);


  //***************************
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Post-process string graph done.";
    PrintElapsed(start_time, check_time, "");
  }
  //***************************
  // Write the string graph to hard dist
  start_time = MyTime();
  //string db_stem = GetFileStem(db_file);
  //string idx_unitig_file = workspace_dir + "/" + db_stem + ".utg";  
  strG.WriteGraph(protein_alphabet, seqs, order, idx_unitig_file);
  //strG.WriteGraph(protein_alphabet, seqs, idx_unitig_file);
  strG.Purge();
  // Collect memory
  for(int idm = 0; idm < num_seqs; ++ idm) {
    delete [] seqs[idm]; 
  }
  delete [] seqs;
  //***************************
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Write string graph unitigs done.";
    PrintElapsed(start_time, check_time, "");
  }
  
  start_time = MyTime();
  SequenceSearch seq_search; 

  string idx_neighbor_file = workspace_dir + "/" + db_stem + ".knb"; 
  seq_search.IndexKmerNeighbor(
      mer_len, protein_alphabet, scoring_function, 
      neighbor_score, idx_neighbor_file
  );
  
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Constructing and writing k-mer index done. ";
    PrintElapsed(start_time, check_time, "");
    cout << "GRASP2-Build: End of program execution." << endl;
    cout << "============================================================" << endl;
  }
  
  //***************************
  return 0;
}
