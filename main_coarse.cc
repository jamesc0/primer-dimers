#include <math.h>       // for pow()
#include <stdio.h>      // for printf, scanf, fgets
#include <stdlib.h>     // for malloc()
#include <string.h>     // for strlen(), strcpy()

#include <algorithm>    // for std::accumulate()
#include <cstring>      // for std::strcpy
#include <fstream>      // for std::ifstream
#include <iostream>     // for std::cout
#include <limits>       // for std::numericlimits
#include <map>          // for std::map
#include <set>          // for std::set

const char* input_file_name = "data/Fergus_all.txt";
const unsigned tail_len = 5;
const unsigned max_mismatches = 1;
const unsigned j = 4;
const unsigned minimum_matching_jmers = 2;
const unsigned minimum_lcs_threshold = 0;

const unsigned number_of_bases = 4;

class PrimerClass {
  std::string name_;
  std::string sequence_;
 public:
  void SetName(std::string str) {
    name_ = str;
  }
  void SetSequence(std::string str) {
    sequence_ = str;
  }
  void Print() {
    std::cout << "name: " << name_ << '\n';
    std::cout << "sequence: " << sequence_ << '\n';
  }
  std::string GetName() {
    return name_;
  }
  std::string GetSequence() {
    return sequence_;
  }
};

typedef struct node {
  int primer_index;
  node* next;
} node_t;

std::vector<char> bases = {'A', 'T', 'C', 'G'};
std::map<char, char> complement_map = {{'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}};
std::map<char, char> next_base = {{'A', 'T'}, {'T', 'C'}, {'C', 'G'}, {'G', 'A'}};
std::map<char, int> base_map = {{'A', 0}, {'T', 1}, {'C', 2}, {'G', 3}};

int hash(const std::string &str);
std::string ReverseComplement(const std::string &src);
std::vector<PrimerClass> ReadInputFile(const std::string &input_file_name);
std::set<std::set<int>> kSubsets(int n, int k);
std::set<std::string> kMismatch(std::string input_str, int max_mismatches);
std::vector<node_t*> LoadTailTable(std::vector<PrimerClass> primers, int tail_len,
    int max_mismatches);
std::vector<std::vector<bool>> MatchTails(std::vector<PrimerClass> primers,
    int tail_len, int max_mismatches);
std::vector<std::vector<bool>> LoadJmerTable(std::vector<PrimerClass> primers,
    unsigned j, bool rc = false);
std::vector<std::vector<unsigned>> MatchJmers(std::vector<PrimerClass> primers,
    int j);
unsigned LcsLen(std::string str1, std::string str2);

int main(int argc, char* argv[]) {
	
  // print parameters
  printf("========================================\n");
  printf("Parameters =============================\n");
  printf("========================================\n");
  std::cout << "input_file_name = " << input_file_name << '\n';
  std::cout << "tail_len = " << tail_len << '\n';
  std::cout << "max_mismatches = " << max_mismatches << '\n';
  std::cout << "j (the length of a jmer) = " << j << '\n';
  std::cout << "minimum_matching_jmers = " << minimum_matching_jmers << '\n';
  std::cout << "minimum_lcs_threshold = " << minimum_lcs_threshold << '\n';
  printf("========================================\n");
  printf("\n");

  // load primers
  auto primers = ReadInputFile(input_file_name);
  
  // investigate lcs
  /*
  int lcs_len;
  std::vector<int> all_matches;
  for (auto i = 0u; i < primers.size()/6; ++i) {
    for (auto j = 0u; j < primers.size()/6; ++j) {
      lcs_len = LcsLen(ReverseComplement(primers[i].GetSequence()), primers[j].GetSequence());
      all_matches.push_back(lcs_len);
    }
  }
  printf("========================================\n");
  printf("lcs statistics =========================\n");
  printf("avg matches = %f\n", std::accumulate(all_matches.begin(), all_matches.end(), 0ull)/(double)all_matches.size());
  printf("all_matches.size() = %f\n", (double)all_matches.size());
  std::vector<int> stats(100000, 0);
  for (int match : all_matches) {
    ++stats[match];
  }
  for (auto i = 0u; i < stats.size(); ++i) {
    if (stats[i] != 0) {
      printf("matches = %-10i, occurrences = %-10i, prob = %-10f\n", i, stats[i], stats[i]/(double)all_matches.size());
    }
  }
  printf("========================================\n");
  */
  
  auto tail_hits = MatchTails(primers, tail_len, max_mismatches);
  auto jmer_hits = MatchJmers(primers, j);
  
  // investigate all conditions
  unsigned sample_size = std::min(1000u, static_cast<unsigned>(primers.size()));
  unsigned tail_count = 0;
  unsigned jmer_count = 0;
  unsigned lcs_count = 0;
  unsigned all_count = 0;
  unsigned conditions_met;
  for (unsigned i = 0; i < sample_size; ++i) {
    for (unsigned j = 0; j < sample_size; ++j) {
      conditions_met = 0;
      if (tail_hits[i][j]) {
        ++tail_count;
        ++conditions_met;
      }
      if (jmer_hits[i][j] >= minimum_matching_jmers) {
        ++jmer_count;
        ++conditions_met;
      }
      /*
      if (LcsLen(ReverseComplement(primers[i].GetSequence()), primers[j].GetSequence()) >= minimum_lcs_threshold) {
        ++lcs_count;
        ++conditions_met;
      }
      if (conditions_met == 3) {
        ++all_count;
      }
      */
    }
  }
  printf("========================================\n");
  printf("Sample =================================\n");
  printf("========================================\n");
  std::cout << "proportion of pairs in small sample successful for tail: " << (double)tail_count / (sample_size * sample_size) << '\n';
  std::cout << "proportion of pairs in small sample successful for jmer: " << (double)jmer_count / (sample_size * sample_size) << '\n';
  //std::cout << "proportion of pairs in small sample successful for all: " << (double)all_count / (sample_size * sample_size) << '\n';
  printf("========================================\n");
  printf("\n");

  // final results
  printf("========================================\n");
  printf("Results: primer dimer candidates =======\n");
  printf("========================================\n");
  int count = 0;
  bool found;
  for (auto i = 0u; i < primers.size(); ++i) {
    found = false;
    for (auto j = 0u; j < primers.size(); ++j) {
      if (!tail_hits[i][j]) continue;
      if (jmer_hits[i][j] < minimum_matching_jmers) continue;
      //if (LcsLen(ReverseComplement(primers[i].GetSequence()), primers[j].GetSequence()) < minimum_lcs_threshold) continue;
      if (!found) {
        std::cout << primers[i].GetName() << " : ";
        std::cout << primers[j].GetName();
        found = true;
      } else {
        std::cout << ", " << primers[j].GetName();
      }
      ++count;
      //if (count % 100 == 0) std::cout << "count = " << count << '\n';
    }
    std::cout << '\n';
  }
  printf("========================================\n");
  printf("\n");
  
  std::cout << count << '\n';
  std::cout << "proportion of hits = " << (double)count / (primers.size() * primers.size()) << '\n';

  return 0;
}

int hash(const std::string &str) {
  int ret_val = 0;
  for (auto i = 0u; i < str.size(); ++i) {
    ret_val += pow(number_of_bases, i) * base_map[str[i]];
  }
  return ret_val;
}

std::string ReverseComplement(const std::string &src) {
  int len = src.size();
  std::string ret_str = src;
  for (int i = 0; i < len; ++i) {
    ret_str[len - 1 - i] = complement_map[src[i]];
  }
  return ret_str;
}

std::vector<PrimerClass> ReadInputFile(const std::string &input_file_name) {
  std::ifstream instream(input_file_name);
  if (!instream.is_open()) {
    std::cout << "Could not open input file.\n";
    std::exit(EXIT_FAILURE);
  }
  std::vector<PrimerClass> primers;
  std::string name;
  std::string sequence;
  for (;;) {
    std::getline(instream, name, ',');
    std::getline(instream, sequence, ',');
    std::transform(sequence.begin(), sequence.end(), sequence.begin(),
      ::toupper);
    instream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if (instream.eof()) {
      break;
    } else {
      PrimerClass primer;
      primer.SetName(name);
      primer.SetSequence(sequence);
      primers.push_back(primer);
    }
  }
  instream.close();
  return primers;
}
std::set<std::set<int>> kSubsets(int n, int k) {
  // returns all k-subsets of {0, 1, ..., n-1}
  if (k < 1) printf("ERROR: in kSubsets k must be >= 1");
  if (k > n) printf("ERROR: in kSubsets k must be <= n");

  std::set<std::set<int>> ret_set;
  std::set<int> tmp_set;
  std::vector<int> include_indexes;
  for (int i = 0; i < k; ++i) include_indexes.push_back(i);

  while (1) {
    tmp_set.clear();
    for (int i : include_indexes) tmp_set.insert(i);
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
    for (std::set<int> include_indexes : s) {
      tmp_str = input_str;
      include_indexes_sorted.clear();
      for (int i : include_indexes) include_indexes_sorted.push_back(i);
      std::sort(include_indexes_sorted.begin(), include_indexes_sorted.end());
      if (include_indexes_sorted[0] == 0) {
        continue;  // the 0th base (the 3' end) must match
      }
      for (int i : include_indexes_sorted) tmp_str[i] = 'A';
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

std::vector<node_t*> LoadTailTable(std::vector<PrimerClass> primers, int tail_len,
    int max_mismatches) {
  std::vector<node_t*> table;
  for (int i = 0; i < pow(number_of_bases, tail_len); ++i) table.push_back(nullptr);
  std::string primer;
  std::string tail;
  std::string tail_rc;
  std::set<std::string> similar_sequences;
  int hash_val;
  int primer_len;
  for (int i = 0; i < static_cast<int>(primers.size()); ++i) {
    primer = primers[i].GetSequence();
    primer_len = primer.size();
    tail = primer.substr(primer_len - tail_len, tail_len);
    tail_rc = ReverseComplement(tail);
    similar_sequences = kMismatch(tail_rc, max_mismatches);
    for (std::string sequence : similar_sequences) {
      hash_val = hash(sequence);
      node_t* tmp_node_ptr = (node_t*) malloc(sizeof(node_t));
      tmp_node_ptr->primer_index = i;
      if (table[hash_val] == nullptr) {
        tmp_node_ptr->next = nullptr;
      } else {
        tmp_node_ptr->next = table[hash_val];
      }
      table[hash_val] = tmp_node_ptr;
    }
  }
  return table;
}

std::vector<std::vector<bool>> MatchTails(std::vector<PrimerClass> primers,
    int tail_len, int max_mismatches) {
  std::vector<std::vector<bool>> hit;
  std::vector<node_t*> table = LoadTailTable(primers, tail_len, max_mismatches);
  std::vector<bool> zero_v(primers.size(), false);
  for (auto i = 0u; i < primers.size(); ++i) hit.push_back(zero_v);
  std::string tmp_primer;
  std::string window;
  node_t* tmp_node_ptr;
  for (unsigned i = 0; i < primers.size(); ++i) {
    tmp_primer = primers[i].GetSequence();
    for (unsigned start_index = 0;
        start_index + tail_len != tmp_primer.size();
        ++start_index) {
      window = tmp_primer.substr(start_index, tail_len);
      tmp_node_ptr = table[hash(window)];
      while (tmp_node_ptr != nullptr) {
        hit[i][tmp_node_ptr->primer_index] = true;
        hit[tmp_node_ptr->primer_index][i] = true;
        tmp_node_ptr = tmp_node_ptr->next;
      }
    }
  }
  return hit;
}

std::vector<std::vector<bool>> LoadJmerTable(std::vector<PrimerClass> primers,
    unsigned j, bool rc) {
  std::vector<std::vector<bool>> jmer_table;
  std::vector<bool> false_v(primers.size(), false);
  for (int i = 0; i < pow(number_of_bases, j); ++i) {
    jmer_table.push_back(false_v);
  }
  std::string primer;
  std::string jmer;
  int hash_val;
  for (unsigned i = 0; i < primers.size(); ++i) {
    primer = primers[i].GetSequence();
    if (rc) primer = ReverseComplement(primer);
    std::string jmer;
    for (unsigned start_index = 0; start_index + j <= primer.size(); start_index += j) {
      jmer = primer.substr(start_index, j);
      hash_val = hash(jmer);
      jmer_table[hash_val][i] = true;
    }
  }
  return jmer_table;
}

std::vector<std::vector<unsigned>> MatchJmers(std::vector<PrimerClass> primers,
    int j) {
  std::vector<std::vector<unsigned>> hit;
  auto jmer_table = LoadJmerTable(primers, j, false);
  auto jmer_rc_table = LoadJmerTable(primers, j, true);
  std::vector<unsigned> zero_v(primers.size(), 0);
  for (unsigned i = 0; i < primers.size(); ++i) hit.push_back(zero_v);
  int matches;
  std::vector<int> all_matches;
  for (unsigned i = 0; i < primers.size(); ++i) {
    //if (i%300 == 0) std::cout << i << '\n';
    for (unsigned k = i; k < primers.size(); ++k) {
      matches = 0;
      for (unsigned p = 0; p < jmer_table.size(); ++p) {
        if (jmer_rc_table[p][i] == true &&
            jmer_table[p][k] == true) {
          ++matches;
        }
      }
      hit[i][k] = hit[k][i] = matches;
      all_matches.push_back(matches);
    }
  }

  /*
  printf("avg matches = %f\n", std::accumulate(all_matches.begin(), all_matches.end(), 0ull)/(double)all_matches.size());
  printf("all_matches.size() = %f\n", (double)all_matches.size());
  std::vector<int> stats(jmer_table.size(), 0);
  for (int match : all_matches) {
    ++stats[match];
  }
  for (auto i = 0u; i < stats.size(); ++i) {
    if (stats[i] != 0) {
      printf("matches = %-10i, occurrences = %-10i, prob = %-10f\n", i, stats[i], stats[i]/(double)all_matches.size());
    }
  }
  */

  return hit;
}

unsigned LcsLen(std::string str1, std::string str2) {
  unsigned lcs_len = 0;
  unsigned len1 = str1.size();
  unsigned len2 = str2.size();
  unsigned rows = len2 + 1;
  unsigned cols = len1 + 1;
  std::vector<unsigned> zero_v(cols, 0);
  std::vector<std::vector<unsigned>> table;
  for (unsigned row = 0; row < rows; ++row) table.push_back(zero_v);
  for (unsigned row = 1; row < rows; ++row) {
    for (unsigned col = 1; col < cols; ++col) {
      if (str1.at(col-1) == str2.at(row-1)) {
        table.at(row).at(col) = table.at(row-1).at(col-1)+1;
        lcs_len = std::max(lcs_len, table.at(row).at(col));
      } else {
        table.at(row).at(col) = 0;
      }
    }
  }
  return lcs_len;
}
