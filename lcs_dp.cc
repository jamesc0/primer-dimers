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
const int max_tail_len = 12;

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


std::string ReverseComplement(std::string src) {
  int len = src.size();
  std::string ret_str = src;
  for (int i = 0; i < len; ++i) {
      ret_str[len - 1 - i] = complement_map[src[i]];
  }
  return ret_str;
}

int lcs_dp(std::string str1, std::string str2) {
  int ret_val = 0;
	int len1 = str1.size();
  int len2 = str2.size();
  int rows = len2 + 1;
  int cols = len1 + 1;
  std::vector<int> zero_vect(cols, 0);
  std::vector<std::vector<int>> table;
  for (int row = 0; row < rows; ++row) table.push_back(zero_vect);
  for (int row = 1; row < rows; ++row) {
    for (int col = 1; col < cols; ++col) {
      if (str1.at(col-1) == str2.at(row-1)) {
        table.at(row).at(col) = table.at(row-1).at(col-1)+1;
        ret_val = std::max(ret_val, table.at(row).at(col));
      } else {
        table.at(row).at(col) = 0;
      }
    }
  }

  // print table
  /*
  for (int row = 0; row <= len2; ++row) {
    for (int col = 0; col <= len1; ++col) {
      printf("%i ", table.at(row).at(col));
    }
    printf("\n");
  }
  printf("ret_val = %i\n", ret_val);
  */

  return ret_val;
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

  // create vector of ints
  std::vector<std::vector<int>> hit;
  std::vector<int> zero_int_vect (number_of_primers, 0);
  for (auto i = 0; i < number_of_primers; ++i) hit.push_back(zero_int_vect);

  for (auto i = 0; i < number_of_primers; ++i) {
  //printf("%i\n", i);
    for (auto j = i; j < number_of_primers; ++j) {
      hit[j][i] = hit[i][j] = lcs_dp(ReverseComplement(primers[i]), primers[j]);
    }
  }

  std::vector<int> all_lcs;
  for (auto i = 0; i < number_of_primers; ++i) {
    for (auto j = 0; j < number_of_primers; ++j) {
      all_lcs.push_back(hit[i][j]);
    }
  }

  printf("avg lcs = %f\n", std::accumulate(all_lcs.begin(), all_lcs.end(), 0ull)/(double)all_lcs.size());
  printf("all_lcs.size() = %f\n", (double)all_lcs.size());
  std::vector<int> stats(primer_len + 2, 0);
  for (int lcs : all_lcs) {
    ++stats[lcs];
  }
  for (auto i = 0u; i < stats.size(); ++i) {
    if (stats[i] != 0) {
      printf("lcs = %-10i, occurrences = %-10i, prob = %-10f\n", i, stats[i], stats[i]/(double)all_lcs.size());
    }
  }


  // free memory from input
  for (int i = 0; i < number_of_primers; ++i) free(primers[i]);
  free(primers);

  // close files
  fclose(infile);
  fclose(outfile);

  return 0;
}
