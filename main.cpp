#include <algorithm>    // for std::accumulate()
#include <cstring>      // for std::strcpy
#include <iostream>     // for std::cout
#include <map>          // for std::map
#include <math.h>       // for pow()
#include <set>          // for std::set
#include <stdio.h>      // for printf, scanf, fgets
#include <stdlib.h>
#include <string.h>     // for strlen(), strcpy()

const char* infile_name = "data/test_data_primers_4000_25.txt";
const char* outfile_name = "data/test_data_primers_4000_25_out.txt";
const int max_number_of_primers = 4000;
const int primer_len = 25;
const int min_tail_len = 5;
const int max_tail_len = 7;

const int number_of_bases = 4;
const int hash_base = number_of_bases;
const int hash_table_size = pow(hash_base, max_tail_len);

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
  if (ret_val >= hash_table_size) printf("hashing error, ret_val >= hash_table_size");
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
    while (tmp_node_ptr != NULL) {
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

void LoadHashTable(node_t** hash_table, char** primers, int number_of_primers, int tail_len, int max_mismatches) {
  std::string primer;
  std::string tail;
  std::string tail_rc;
  std::set<std::string> similar_sequences;
  int hash_val;
  for (int i = 0; i < number_of_primers; ++i) {
    primer = primers[i];
    tail = primer.substr(primer_len - tail_len, tail_len);
    tail_rc = ReverseComplement(tail);
    similar_sequences = kMismatch(tail_rc, max_mismatches);
    for (std::string sequence: similar_sequences) {
      hash_val = hash(sequence);
      node_t* tmp_node_ptr = (node_t*) malloc(sizeof(node_t));
      tmp_node_ptr->primer_index = i;
      if (hash_table[hash_val] == NULL) {
        tmp_node_ptr->next = NULL;
      } else {
        tmp_node_ptr->next = hash_table[hash_val];
      }
      hash_table[hash_val] = tmp_node_ptr;
    }
  }
}

void ClearHashTable(node** hash_table, int hash_table_size) {
  node_t* node_ptr_1;
  node_t* node_ptr_2;
  for (auto i = 0; i < hash_table_size; ++i) {
    node_ptr_1 = hash_table[i];
    while (node_ptr_1 != NULL) {
      node_ptr_2 = node_ptr_1;
      node_ptr_1 = node_ptr_1->next;
      free(node_ptr_2);
      node_ptr_2 = NULL;
    }
    hash_table[i] = NULL;
  }
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

std::vector<std::vector<int>> SlideWindow(char** primers, int number_of_primers, node_t** hash_table, int tail_len) {
  // initialise hit vector
  std::vector<std::vector<int>> hit;
  std::vector<int> zero_vect(number_of_primers, 0);
  for (int i = 0; i < number_of_primers; ++i) hit.push_back(zero_vect);

  char tmp_primer[primer_len + 1];
  node_t* tmp_node_ptr;
  for (int i = 0; i < number_of_primers; ++i) {
    strcpy(tmp_primer, primers[i]);
    char* window = tmp_primer + strlen(tmp_primer) - tail_len;
    for (; window != tmp_primer; --window) {
      if (window == tmp_primer) break;
      tmp_node_ptr = hash_table[hash(window)];
      while(tmp_node_ptr != NULL) {
        hit[i][tmp_node_ptr->primer_index] = hit[tmp_node_ptr->primer_index][i] = 1;
        tmp_node_ptr = tmp_node_ptr->next;
      }
      tmp_primer[strlen(tmp_primer) - 1] = '\0';
    }
  }
  return hit;
}

double P_substring(int substring_len, int target_len) {
  return (target_len - substring_len + 1) * pow(1.0/number_of_bases, substring_len);
}

int main(int argc, char* argv[]) {
  // open files
  FILE* infile = fopen(infile_name, "r");
  FILE* outfile = fopen(outfile_name, "w");
  
  // check if infile exists
  if (infile == NULL) {
    printf("can't open infile\n");
    exit(EXIT_FAILURE);
  }

  // read input file
  int number_of_primers;
  char **primers = (char**) malloc(max_number_of_primers * sizeof(char*));
  fscanf(infile, "%d\n", &number_of_primers);
  
  char *tmp;
  for (int i = 0; i < number_of_primers; ++i) {
    tmp = (char*) malloc(1 + primer_len * sizeof(char));
    fgets(tmp, primer_len+2, infile);
    if (tmp[strlen(tmp) - 1] == '\n') tmp[strlen(tmp) - 1] = '\0';
    primers[i] = tmp;
  }
  printf("number_of_primers = %i\n", number_of_primers);
  printf("primer_len = %i\n", primer_len);

  // create hash table
  node_t** hash_table = (node_t**) malloc(hash_table_size * sizeof(node_t*));
  for (auto i = 0; i < hash_table_size; ++i) hash_table[i] = NULL;
  
  // test probabilies for different tail_len and max_mismatches
  printf("+------------+--------+-----------+-------------+-------------+\n");
  printf("| max        | tail   | total     | actual      | expected    |\n");
  printf("| mismatches | length | entries   | hit         | hit         |\n");
  printf("|            |        | in        | probability | probability |\n");
  printf("|            |        | hashtable |             |             |\n");
  printf("+------------+--------+-----------+-------------+-------------+\n");
  std::vector<std::vector<int>> hit;
  for (auto tail_len = min_tail_len; tail_len <= max_tail_len; ++tail_len) {
    for (auto max_mismatches = 0; max_mismatches <= 3; ++max_mismatches) {
      printf("| %-11i| %-7i|", max_mismatches, tail_len);
      LoadHashTable(hash_table, primers, number_of_primers, tail_len, max_mismatches);
      PrintHashTableStatistics(hash_table, hash_table_size);
      hit.clear();
      hit = SlideWindow(primers, number_of_primers, hash_table, tail_len);
      PrintHitStatistics(hit, number_of_primers);
      ClearHashTable(hash_table, hash_table_size);
      if (max_mismatches == 0) {
        printf(" %-12f|", P_substring(tail_len, primer_len));
      }
      printf("\n");
      printf("+------------+--------+-----------+-------------+-------------+\n");
    }
  }

  // free hash table
  free(hash_table);

  // free memory from input
  for (int i = 0; i < number_of_primers; ++i) free(primers[i]);
  free(primers);

  // close files
  fclose(infile);
  fclose(outfile);

  return 0;
}
