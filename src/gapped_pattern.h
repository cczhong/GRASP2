#ifndef _GAPPED_PATTERN_
#define _GAPPED_PATTERN_

#include <string>
#include <cstring>

class GappedPattern  {
 public:
  explicit GappedPattern(void) {
    InitPattern();
  }
  ~GappedPattern(void) {}
  // returns the string without any gap
  std::string GetUngappedStr(const int pid, char *seq)  {
    std::string s = "";    
    if(strlen(seq) >= pattern[pid].length()) {
      for(int i = 0; i < pattern[pid].length(); ++ i) {
        if(pattern[pid][i] == '1') s += seq[i];
      }
    }
    return s;
  }

  std::string GetUngappedStr(const int pid, const std::string &seq, const int begin)  {
    std::string s = "";    
    if(seq.length() - begin >= pattern[pid].length()) {
      for(int i = 0; i < pattern[pid].length(); ++ i) {
        if(pattern[pid][i] == '1') s += seq[begin + i];
      }
    }
    return s;
  }
  // returns the length of the pattern
  int GetPatternLen(const int pid)  {
    if(pid >= pattern.size()) return -1;
    return pattern[pid].length();
  }
  int GetPatternWeight(void)  { return pattern_weight; }
  int GetNumPatterns(void)  { return pattern.size();  }

 private:
  std::vector<std::string> pattern;
  int pattern_weight;
  
  void InitPattern(void)  {

    /*
    pattern.push_back("111010010100110111");
    pattern.push_back("111100110010100001011");
    pattern.push_back("110100001100010101111");
    pattern.push_back("1110111010001111");
    */

    pattern.push_back("1110110111");
    pattern.push_back("111010001000111");
    pattern.push_back("101001100010001011");
    pattern.push_back("101010000001100111");
    pattern.push_back("1100100100010010101");
    pattern.push_back("1011000001001010011");
    pattern.push_back("1110010000010101001");
    pattern.push_back("1101000010010000111");
    pattern.push_back("11000100001101000101");
    pattern.push_back("10110000100100100011");
    pattern.push_back("11000100010000101011");
    pattern.push_back("11001001000000110101");
    pattern.push_back("11010000110000010011");
    pattern.push_back("10010100001010001011");
    pattern.push_back("10100011000010010011");
    pattern.push_back("11000010101000001011");

    pattern_weight = 8;
    return;
  }
  
};

#endif
