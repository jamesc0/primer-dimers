#ifndef SEQUENCE_FUNCTIONS_H
#define SEQUENCE_FUNCTIONS_H

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
const int max_tail_len = 5;

const int number_of_bases = 4;
const int hash_base = number_of_bases;
const int hash_table_size = pow(hash_base, max_tail_len);

std::vector<char> bases = {'A', 'T', 'C', 'G'};
std::map<char, char> complement_map = {{'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}};
std::map<char, char> next_base = {{'A', 'T'}, {'T', 'C'}, {'C', 'G'}, {'G', 'A'}};
std::map<char, int> base_map = {{'A', 0}, {'T', 1}, {'C', 2}, {'G', 3}};

#endif
