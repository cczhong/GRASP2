#ifndef _GAPPED_PATTERN_
#define _GAPPED_PATTERN_

#include <string>

class GappedPattern  {
 public:
  explicit GappedPattern(void) {
    InitPattern();
  }
  ~GappedPattern(void) {}
  // returns the string without any gap
  std::string GetUngappedStr(const int pid, const std::string &s, const int pos);
  // returns the length of the pattern
  int GetPatternLen(const int pid)  {
    if(pid >= pattern.size()) return -1;
    return pattern[pid].length();
  }
  int GetPatternWeight(void)  { return pattern_weight; }
 private:
  std::vector<std::string> pattern;
  int pattern_weight;
  
  void InitPattern(void)  {
    pattern.push_back("111010010100110111");
    pattern.push_back("111100110010100001011");
    pattern.push_back("110100001100010101111");
    pattern.push_back("1110111010001111");
    pattern_weight = 11;
    return;
  }
  
};

#endif
