#include <stdio.h>      /* for printf, scanf, fscanf */
#include <string.h>     /* for strlen(), strcpy() */
#include <stdlib.h>
#include <algorithm>    /* for std::accumulate() */
#include <map>          /* for std::map */
#include <math.h>       /* for pow() */
#include <set>          /* for std::set */

const char* infile_name = "data/test_data_primers_4000_25.txt";
const int max_primer_length = 50;
const int max_number_of_primers = 10000;
const int max_tail_len = 5;
const int hash_base = 4;
const int hash_table_size = pow(hash_base, max_tail_len);
const int mismatches = 0;

typedef struct node {
  int primer_index;
  node* next;
} node_t;

std::vector<char> bases = {'A', 'T', 'C', 'G'};
std::map<char, char> complement_map = {{'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}};
std::map<char, int> base_map = {{'A', 0}, {'T', 1}, {'C', 2}, {'G', 3}};

int hash(char* str) {
  int val = 0;
  for (int i = 0; str[i] != '\0'; ++i) val += pow(hash_base, i) * base_map[str[i]];
  val %= hash_table_size;
  return val;
}

void ReverseComplement(char* dst, char* src) {
  int len = strlen(src);
  dst[len] = '\0';
  for (int i = 0; src[i] != '\0'; ++i) dst[len - 1 - i] = complement_map[src[i]];
}

std::set<std::set<int>> kSubsets(int n, int k) {
  /* returns all k-subsets of {0, 1, ..., n-1} */
  if (k < 1) printf("ERROR: in kSubsets k must be >= 1");
  if (k > n) printf("ERROR: in kSubsets k must be < n");
  
  std::set<std::set<int>> ret_set;
  std::set<int> tmp_set;
  std::vector<int> include_indexes;
  for (int i = 0; i < k; ++i) include_indexes.push_back(i);

  while(1) {
    tmp_set.clear();
    for (int i: include_indexes) tmp_set.insert(i);
    ret_set.insert(tmp_set);
    if (include_indexes[0] == (n-1)-(k-1)) break;
    
    /* update include_indexes */
    for (int i = include_indexes.size() - 1; i >= 0; --i) {
      if (i == include_indexes.size() - 1 && include_indexes[i] < n-1) {
        ++include_indexes[i];
        break;
      } else if (i < include_indexes.size() - 1 && include_indexes[i] < include_indexes[i+1] - 1) {
        ++include_indexes[i];
        for (int j = i + 1; j < include_indexes.size(); ++j) {
          include_indexes[j] = include_indexes[j-1]+1;
        }
        break;
      }
    }
  }
  return ret_set;
}

std::set<int> MismatchHash(char* input_string, int mismatches) {
  if (mismatches != 1 && mismatches != 0) printf("ERROR: this can only handle 0 or 1 mismatches right now");
  std::set<int> ret_set;
  char tmp[max_tail_len + 1];
  int len = strlen(input_string);
  if (mismatches == 0) {
    ret_set.insert(hash(input_string));
    return ret_set;
  }
  for (int i = 0; i < len; ++i) {
    /* i is the index of the mismatch */
    strcpy(tmp, input_string);
    for (char c: bases) {
      tmp[i] = c;
      ret_set.insert(hash(tmp));
    }
  }
  return ret_set;
}

void PrintHashTableStatistics(node** hash_table, int hash_table_size) {
  int count, max_count = 0, total_count = 0;
  node_t* tmp_node_ptr;
  for (int i = 0; i < hash_table_size; ++i) {
    count = 0;
    tmp_node_ptr = hash_table[i];
    while (tmp_node_ptr != NULL) {
      ++count;
      tmp_node_ptr = tmp_node_ptr->next;
    }
    if (count > max_count) max_count = count;
    total_count += count;
  }
  printf("There are %i total entries in the hash table\n", total_count);
  printf("The most entries in one index is %i\n", max_count);
}

void LoadHashTable(node_t** hash_table, char** primers, int number_of_primers, int tail_len, int mismatches) {
  char tail_rc[max_tail_len + 1];   /* reverse complement of the tail */
  std::set<int> hash_vals;
  for (int i = 0; i < number_of_primers; ++i) {
    ReverseComplement(tail_rc, primers[i] + strlen(primers[i]) - tail_len);
    hash_vals = MismatchHash(tail_rc, mismatches);
    for (int hash_val: hash_vals) {
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
  for (int i = 0; i < hash_table_size; ++i) {
    node_ptr_1 = hash_table[i];
    while (node_ptr_1 != NULL) {
      node_ptr_2 = node_ptr_1;
      node_ptr_1 = node_ptr_1->next;
      free(node_ptr_2);
    }
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
  double avg_hits = std::accumulate(all_hits.begin(), all_hits.end(), 0ULL) / (double)all_hits.size();
  printf("the max_hits for any primer = %i, for primer %i\n", max_hits, max_hits_index);
  printf("the avg_hits for each primer = %f\n", avg_hits);
  printf("the avg percentage hits for each primer = %f%%\n", 100*avg_hits/number_of_primers);
}

std::vector<std::vector<int>> SlideWindow(char** primers, int number_of_primers, node_t** hash_table, int tail_len) {
  /* initialise hit vector */
  std::vector<std::vector<int>> hit;
  std::vector<int> zero_vect(number_of_primers, 0);
  for (int i = 0; i < number_of_primers; ++i) hit.push_back(zero_vect);
  
  char tmp_primer[max_primer_length + 1];
  node_t* tmp_node_ptr;
  for (int i = 0; i < number_of_primers; ++i) {
    strcpy(tmp_primer, primers[i]);
    char* window;
    window = tmp_primer + strlen(tmp_primer) - tail_len;
    //printf("%s\n", window);
    while (window != tmp_primer) {
      //printf("%s\n", window);
      tmp_node_ptr = hash_table[hash(window)];

      /* travel through the index (which is a linked list) of the hash table */
      while(tmp_node_ptr != NULL) {
        hit[i][tmp_node_ptr->primer_index] = hit[tmp_node_ptr->primer_index][i] = 1;
        tmp_node_ptr = tmp_node_ptr->next;
      }

      --window;
      tmp_primer[strlen(tmp_primer) - 1] = '\0';
    }
  }
  return hit;
}

int main(int argc, char* argv[]) {
  /* open file */
  FILE * infile;
  infile = fopen(infile_name, "r");
  
  /* check if file exists */
  if (infile == NULL) {
    printf("can't open infile\n");
    exit(EXIT_FAILURE);
  }

  /* read input file */
  int number_of_primers;
  char **primers = (char**) malloc(max_number_of_primers * sizeof(char*));
  fscanf(infile, "%d\n", &number_of_primers);
  //printf("number_of_primers = %i\n", *number_of_primers);
  char *tmp;
  for (int i = 0; i < number_of_primers; ++i) {
    tmp = (char*) malloc(1 + max_primer_length * sizeof(char));
    fgets(tmp, max_primer_length, infile);
    if (tmp[strlen(tmp) - 1] == '\n') tmp[strlen(tmp) - 1] = '\0';
    primers[i] = tmp;
  }

  /* close infile */
  fclose(infile);

  /* create hash table */
  node_t** hash_table = (node_t**) malloc(hash_table_size * sizeof(node_t*));
  for (int i = 0; i < hash_table_size; ++i) hash_table[i] = NULL;

  int tail_len = 5;
  LoadHashTable(hash_table, primers, number_of_primers, tail_len, mismatches);
  PrintHashTableStatistics(hash_table, hash_table_size);
  std::vector<std::vector<int>> hit;
  hit = SlideWindow(primers, number_of_primers, hash_table, tail_len);
  PrintHitStatistics(hit, number_of_primers);
  ClearHashTable(hash_table, hash_table_size);

  /* free hash table */
  free(hash_table);

  /* free memory from input*/
  for (int i = 0; i < number_of_primers; ++i) free(primers[i]);
  free(primers);

  return 0;
}
