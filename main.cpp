#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <map>
#include <math.h>

const char* infile_name = "data/test_data_primers_4000_25.txt";
const int base = 4;
const int tail_len = 5;
const int hash_table_size = pow(base, tail_len);
const int max_sequence_length = 50;
const int max_number_of_primers = 50000;

typedef struct node {
  int primer_index;
  node* next;
} node_t;

std::map<char, char> complement_map = {{'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}};
std::map<char, int> base_map = {{'A', 0}, {'T', 1}, {'C', 2}, {'G', 3}};

int hash(char* str) {
  int val = 0;
  int i;
  for (i = 0; str[i] != '\0'; ++i) val += pow(base, i) * base_map[str[i]];
  val %= hash_table_size;
  return val;
}

void ReverseComplement(char* dst, char* src) {
  int len = strlen(src);
  dst[len] = '\0';
  int i;
  for (i = 0; src[i] != '\0'; ++i) dst[len - 1 - i] = src[i];
}

int main(int argc, char* argv[]) {
  int i; /* counter variables */

  /* open file */
  FILE * infile;
  infile = fopen(infile_name, "r");
  
  /* check if file exists */
  if (infile == NULL) {
    printf("can't open infile\n");
    exit(EXIT_FAILURE);
  }
  
  /* read contents of file */
  int number_of_primers;
  fscanf(infile, "%d\n", &number_of_primers);
  printf("number_of_primers = %i\n", number_of_primers);
  char **primers = (char**) malloc(max_number_of_primers * sizeof(char*));
  char *tmp;
  for (i = 0; i < number_of_primers; ++i) {
    tmp = (char*) malloc(max_sequence_length * sizeof(char));
    fgets(tmp, max_sequence_length, infile);
    if (tmp[strlen(tmp) - 1] == '\n') tmp[strlen(tmp) - 1] = '\0';
    primers[i] = tmp;
  }

  /* create hash table */
  node_t** hash_table = (node_t**) malloc(hash_table_size * sizeof(node_t*));
  for (i = 0; i < hash_table_size; ++i) hash_table[i] = NULL;

  /* load hash table */
  int hash_val;               /* output of the hash function */
  char tail_rc[tail_len + 1]; /* reverse complement of the tail */
  for (i = 0; i < number_of_primers; ++i) {
    node_t* tmp_node_ptr = (node_t*) malloc(sizeof(node_t));
    tmp_node_ptr->primer_index = i;
    ReverseComplement(tail_rc, primers[i] + strlen(primers[i]) - tail_len);
    hash_val = hash(tail_rc);
    if (hash_table[hash_val] == NULL) {
      tmp_node_ptr->next = NULL;
    } else {
      tmp_node_ptr->next = hash_table[hash_val];
    }
    hash_table[hash_val] = tmp_node_ptr;
  }

  /* free hash table */
  node_t* node_ptr_1;
  node_t* node_ptr_2;
  for (i = 0; i < hash_table_size; ++i) {
    node_ptr_1 = hash_table[i];
    while (node_ptr_1 != NULL) {
      node_ptr_2 = node_ptr_1;
      node_ptr_1 = node_ptr_1->next;
      free(node_ptr_2);
    }
  }
  free(hash_table);

  /* free memory from input*/
  for (i = 0; i < number_of_primers; ++i) free(primers[i]);
  free(primers);

  /* close infile */
  fclose(infile);

  return 0;
}
