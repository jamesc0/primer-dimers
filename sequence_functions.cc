#include "sequence_functions.h"

int hash(std::string str) {
  int ret_val = 0;
  int len = str.size();
  for (int i = 0; i < len; ++i) {
    ret_val += pow(hash_base, i) * base_map[str[i]];
  }
  //if (ret_val >= hash_table_size) printf("hashing error, ret_val >= hash_table_size");
  return ret_val;
}

std::string ReverseComplement(std::string src) {
  int len = src.size();
  std::string ret_str = src;
  for (int i = 0; i < len; ++i) {
      ret_str[len - 1 - i] = complement_map[src[i]];
  }
  return ret_str;
}

std::set<std::set<int>> kSubsets(int n, int k) {
  // returns all k-subsets of {0, 1, ..., n-1}
  if (k < 1) printf("ERROR: in kSubsets k must be >= 1");
  if (k > n) printf("ERROR: in kSubsets k must be <= n");
  
  std::set<std::set<int>> ret_set;
  std::set<int> tmp_set;
  std::vector<int> include_indexes;
  for (int i = 0; i < k; ++i) include_indexes.push_back(i);

  while(1) {
    tmp_set.clear();
    for (int i: include_indexes) tmp_set.insert(i);
    ret_set.insert(tmp_set);
    if (include_indexes[0] == (n-1)-(k-1)) break;
    
    // update include_indexes
    for (auto i = include_indexes.size() - 1; i >= 0; --i) {
      if (i == include_indexes.size() - 1 && include_indexes[i] < n-1) {
        ++include_indexes[i];
        break;
      } else if (i < include_indexes.size() - 1 && include_indexes[i] < include_indexes[i+1] - 1) {
        ++include_indexes[i];
        for (auto j = i + 1; j < include_indexes.size(); ++j) {
          include_indexes[j] = include_indexes[j-1]+1;
        }
        break;
      }
    }
  }
  return ret_set;
}

std::set<std::string> kMismatch(std::string input_str, int max_mismatches) {
  // returns a set of strings which have less than or equal to
  // max_mismatches with input_str
  std::set<std::string> ret_set;
  int len = input_str.size();
  if (max_mismatches < 0) {
    printf("ERROR: in kMismatch, max_mismatches must be nonnegative.");
    exit(EXIT_FAILURE);
  } else if (max_mismatches == 0) {
    ret_set.insert(input_str);
  } else {
    if (max_mismatches > len) max_mismatches = len;
    std::set<std::set<int>> s = kSubsets(len, max_mismatches);
    std::vector<int> include_indexes_sorted;
    std::string tmp_str;
    for (std::set<int> include_indexes: s) {
      tmp_str = input_str;
      include_indexes_sorted.clear();
      for (int i: include_indexes) include_indexes_sorted.push_back(i);
      std::sort(include_indexes_sorted.begin(), include_indexes_sorted.end());
      if (include_indexes_sorted[0] == 0) {
        continue; // the 0th base (the 3' end) must match
      }
      for (int i: include_indexes_sorted) tmp_str[i] = 'A';
      for (int j = 0; j < pow(number_of_bases, max_mismatches); j++) {
        ret_set.insert(tmp_str);
        for (auto i = 0u; i < include_indexes_sorted.size(); ++i) {
          if (tmp_str[include_indexes_sorted[i]] != 'G') {
            tmp_str[include_indexes_sorted[i]] = next_base[tmp_str[include_indexes_sorted[i]]];
            for (auto u = 0u; u < i; ++u) tmp_str[include_indexes_sorted[u]] = 'A';
            break;
          }
        }
      }
    }
  }
  return ret_set;
}
