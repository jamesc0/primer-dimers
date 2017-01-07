#include <algorithm>    // for std::accumulate()
#include <bitset>       // for std::bitset()
#include <cstring>      // for std::strcpy
#include <iostream>     // for std::cout
#include <map>          // for std::map
#include <math.h>       // for pow()
#include <set>          // for std::set
#include <stdio.h>      // for printf, scanf, fgets
#include <stdlib.h>
#include <string.h>     // for strlen(), strcpy()

const char* infile_name = "data/test_data_primers_4000_25.txt";
const char* outfile_name = "data/test_data_primers_4000_25_jmer_out.txt";
const int max_number_of_primers = 4000;
const int primer_len = 25;
const int j = 4; // the length of a j-mer

const int number_of_bases = 4;
const int hash_base = number_of_bases;
const int hash_table_size = pow(hash_base, j);

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
  if (ret_val >= hash_table_size) {
    printf("hashing error, ret_val >= hash_table_size");
    printf("hash val = %i\n", ret_val);
  }
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


void LoadJmerHashTable(std::array<std::vector<bool>, hash_table_size> &jmer_hash_table, char** primers, int number_of_primers, int j) {
  std::string primer;
  std::string jmer;
  char tmp_primer[primer_len + 1];
  int hash_val;
  for (int i = 0; i < number_of_primers; ++i) {
    strcpy(tmp_primer, primers[i]);
    char* jmer = tmp_primer + strlen(tmp_primer) - j;
    for (; jmer != tmp_primer; --jmer) {
      //std::cout << "jmer = " << jmer << '\n';
      //std::cout << "jmer_rc = " << jmer_rc << '\n';
      hash_val = hash(jmer);
      jmer_hash_table[hash_val][i] = true;
      tmp_primer[strlen(tmp_primer) - 1] = '\0';
    }
  }
}

void LoadJmerRcHashTable(std::array<std::vector<bool>, hash_table_size> &jmer_rc_hash_table, char** primers, int number_of_primers, int j) {
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
  int hits, max_hits = 0;
  std::vector<int> all_hits;
  for (int i = 0; i < number_of_primers; ++i) {
    hits = 0;
    for (int j = 0; j < number_of_primers; ++j) {
      if (hit[i][j]) ++hits;
    }
    if (hits > max_hits) {
      max_hits = hits;
    }
    all_hits.push_back(hits);
  }
  double avg_hits = std::accumulate(all_hits.begin(), all_hits.end(), 0ull) / (double)all_hits.size();
  //printf("the max_hits for any primer = %i, for primer %i\n", max_hits, max_hits_index);
  //printf("the avg_hits for each primer = %f\n", avg_hits);
  //printf("the avg percentage hits for each primer = %f%%\n", 100*avg_hits/number_of_primers);
  printf(" %-12f|", avg_hits/number_of_primers);
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
    tmp = (char*) malloc(2 + primer_len * sizeof(char));
    fgets(tmp, primer_len+2, infile);
    if (tmp[strlen(tmp) - 1] == '\n') tmp[strlen(tmp) - 1] = '\0';
    primers[i] = tmp;
  }
  printf("number_of_primers = %i\n", number_of_primers);
  printf("primer_len = %i\n", primer_len);

  // create hash tables
  std::array<std::vector<bool>, hash_table_size> jmer_hash_table;
  std::array<std::vector<bool>, hash_table_size> jmer_rc_hash_table;
  std::vector<bool> zero_vect(number_of_primers, false);
  for (auto i = 0u; i < hash_table_size; ++i) {
    jmer_hash_table[i] = zero_vect;
    jmer_rc_hash_table[i] = zero_vect;
  }
  
  // load hash tables
  LoadJmerHashTable(jmer_hash_table, primers, number_of_primers, j);
  LoadJmerRcHashTable(jmer_rc_hash_table, primers, number_of_primers, j);

  // create vector of ints
  std::vector<std::vector<int>> hit;
  std::vector<int> zero_int_vect (number_of_primers, 0);
  for (auto i = 0; i < number_of_primers; ++i) hit.push_back(zero_int_vect);

  // start calculating
  int matches;
  std::vector<int> all_matches;
  for (auto i = 0; i < number_of_primers; ++i) {
    if (i % 100 == 0) printf("%i\n", i);
    for (auto k = i; k < number_of_primers; ++k) {
      matches = 0;
      for (auto p = 0; p < hash_table_size; ++p) {
        if (jmer_rc_hash_table[p][i] == true && 
            jmer_hash_table[p][k] == true) {
          ++matches;
        }
      }
      hit[i][k] = hit[k][i] = matches;
      all_matches.push_back(matches);
      if (i % 500 == 0 && k % 500 == 0) printf("matches = %i\n", matches);
    }
  }
  printf("avg matches = %f\n", std::accumulate(all_matches.begin(), all_matches.end(), 0ull)/(double)all_matches.size());

  // free memory from input
  for (int i = 0; i < number_of_primers; ++i) free(primers[i]);
  free(primers);

  // close files
  fclose(infile);
  fclose(outfile);

  return 0;
}
