#include <math.h>       // for pow()
#include <stdio.h>      // for printf, scanf, fgets
#include <stdlib.h>     // for malloc()
#include <string.h>     // for strlen(), strcpy()

#include <algorithm>    // for std::accumulate()
#include <cstring>      // for std::strcpy
#include <fstream>
#include <iostream>     // for std::cout
#include <limits>       // for numericlimits
#include <map>          // for std::map
#include <set>          // for std::set

const std::string input_file_name = "data/Fergus26genes_05may15_idt.txt";
const int tail_len = 5;
const int max_mismatches = 1;
const int number_of_bases = 4;
const int hash_base = number_of_bases;

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

int hash(std::string str) {
  int ret_val = 0;
  for (auto i = 0u; i < str.size(); ++i) ret_val += pow(hash_base, i) * base_map[str[i]];
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

void PrintHashTableStatistics(node** hash_table, int hash_table_size) {
  int count, max_count = 0, total_count = 0;
  node_t* tmp_node_ptr;
  for (auto i = 0; i < hash_table_size; ++i) {
    count = 0;
    tmp_node_ptr = hash_table[i];
    while (tmp_node_ptr != nullptr) {
      ++count;
      tmp_node_ptr = tmp_node_ptr->next;
    }
    if (count > max_count) max_count = count;
    total_count += count;
  }
  //printf("There are %i total entries in the hash table\n", total_count);
  //printf("The most entries in one index is %i\n", max_count);
  printf(" %-10i|", total_count);
}

void PrintHitStatistics(std::vector<std::vector<int>> hit, int number_of_primers) {
  int hits, max_hits = 0, max_hits_index;
  std::vector<int> all_hits;
  for (int i = 0; i < number_of_primers; ++i) {
    hits = 0;
    for (int j = 0; j < number_of_primers; ++j) {
      if (hit[i][j]) ++hits;
    }
    if (hits > max_hits) {
      max_hits = hits;
      max_hits_index = i;
    }
    all_hits.push_back(hits);
  }
  double avg_hits = std::accumulate(all_hits.begin(), all_hits.end(), 0ull) / (double)all_hits.size();
  //printf("the max_hits for any primer = %i, for primer %i\n", max_hits, max_hits_index);
  //printf("the avg_hits for each primer = %f\n", avg_hits);
  //printf("the avg percentage hits for each primer = %f%%\n", 100*avg_hits/number_of_primers);
  printf(" %-12f|", avg_hits/number_of_primers);
}

std::vector<PrimerClass> ReadInputFile(std::string input_file_name) {
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

std::vector<node_t*> LoadTable(std::vector<PrimerClass> primers, 
    int tail_len, int max_mismatches) {
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
    for (std::string sequence: similar_sequences) {
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

std::vector<std::vector<int>> MatchTails(std::vector<PrimerClass> primers, 
    int tail_len, int max_mismatches) {
  std::vector<std::vector<int>> hit;
  std::vector<node_t*> table = LoadTable(primers, tail_len, max_mismatches);
  std::vector<int> zero_vect(primers.size(), 0);
  for (auto i = 0u; i < primers.size(); ++i) hit.push_back(zero_vect);
  std::string tmp_primer;
  std::string window;
  node_t* tmp_node_ptr;
  for (int i = 0; i < static_cast<int>(primers.size()); ++i) {
    tmp_primer = primers[i].GetSequence();
    for (int start_index = 0; 
        start_index + tail_len != static_cast<int>(tmp_primer.size());
        ++start_index) {
      window = tmp_primer.substr(start_index, tail_len);
      tmp_node_ptr = table[hash(window)];
      while(tmp_node_ptr != nullptr) {
        hit[i][tmp_node_ptr->primer_index] = 1;
        hit[tmp_node_ptr->primer_index][i] = 1;
        tmp_node_ptr = tmp_node_ptr->next;
      }
    }
  }
  return hit;
}

/*
std::vector<std::vector<bool> LoadJmerHashTable(
    std::array<std::vector<bool>, hash_table_size> &jmer_rc_hash_table, 
    char** primers, int number_of_primers, int j, bool rc = false) {
  std::string primer;
  std::string jmer;
  std::string jmer_rc;
  char tmp_primer[primer_len + 1];
  int hash_val;
  for (int i = 0; i < number_of_primers; ++i) {
    strcpy(tmp_primer, primers[i]);
    char* jmer = tmp_primer + strlen(tmp_primer) - j;
    for (; jmer != tmp_primer; --jmer) {
      //std::cout << "jmer = " << jmer << '\n';
      jmer_rc = ReverseComplement(jmer);
      hash_val = hash(jmer_rc);
      jmer_rc_hash_table[hash_val][i] = true;
      tmp_primer[strlen(tmp_primer) - 1] = '\0';
    }
  }
}
*/

int main(int argc, char* argv[]) {
  auto primers = ReadInputFile(input_file_name);

  auto tail_hits = MatchTails(primers, tail_len, max_mismatches);
  
  for (int i = 0; i < primers.size()/60; ++i) {
    for (int j = 0; j < primers.size()/60; ++j) {
      std::cout << tail_hits[i][j] << ' ';
    }
    std::cout << '\n';
  }


  auto jmer_hits = 1;
  return 0;
}
