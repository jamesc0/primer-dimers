#include <stdio.h>      // for printf()
#include <stdlib.h>     // for malloc()

#include <algorithm>
#include <fstream>      // for std::ifstream
#include <iostream>     // for std::cout
#include <map>          // for std::map

std::vector<char> bases = {'A', 'T', 'C', 'G'};
std::map<char, char> complement_map = {{'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}};

std::string ReverseComplement(const std::string &src) {
  int len = src.size();
  std::string ret_str = src;
  for (int i = 0; i < len; ++i) {
    ret_str[len - 1 - i] = complement_map[src[i]];
  }
  return ret_str;
}

bool ValidSequence(std::string str) {
  for (char c : str) {
    if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
      return false;
    }
  }
  return true;
}

int DpAlignment(std::string seq1, std::string seq2) {
  std::vector<std::vector<std::pair<unsigned, char>>> table;
  // char is the direction you follow back to the (0, 0).
  // L is left, U is down, D is diagonal.
  
  std::cout << "Input sequences:\n";
  std::cout << "seq1 = " << seq1 << '\n';
  std::cout << "seq2 = " << seq2 << '\n';
  std::cout << '\n';

  seq1 = ReverseComplement(seq1);
  std::cout << "Aligning the following sequences:\n";
  std::cout << "seq1_rc = " << seq1 << '\n';
  std::cout << "seq2 = " << seq2 << '\n';
  std::cout << '\n';
  
  unsigned cols = seq1.size() + 1;
  unsigned rows = seq2.size() + 1;
  std::vector<std::pair<unsigned, char>> zero_row(cols, {0, '\0'});
  for (unsigned i = 0; i < rows; ++i) table.push_back(zero_row);
  for (unsigned i = 0; i < cols; ++i) table.at(0).at(i) = {i, 'L'};
  for (unsigned i = 0; i < rows; ++i) table.at(i).at(0) = {i, 'U'};

  // write table
  unsigned min;
  for (unsigned row = 1; row < rows; ++row) {
    for (unsigned col = 1; col < cols; ++col) {
      if (seq1.at(col - 1) == seq2.at(row - 1)) {
        table.at(row).at(col) = {table.at(row - 1).at(col - 1).first, 'D'};
        min = table.at(row).at(col).first;
      } else {
        table.at(row).at(col) = {table.at(row - 1).at(col - 1).first + 1, 'D'};;
        min = table.at(row).at(col).first;
      }
      if (table.at(row - 1).at(col).first + 1 < min) {
        table.at(row).at(col) = {table.at(row - 1).at(col).first + 1, 'U'};
        min = table.at(row).at(col).first;
      }
      if (table.at(row).at(col - 1).first + 1 < min) {
        table.at(row).at(col) = {table.at(row).at(col - 1).first + 1, 'L'};
        min = table.at(row).at(col).first;
      }
    }
  }
  
  // print table
  for (unsigned j = 0; j < rows; ++j) {
    for (unsigned i = 0; i < cols; ++i) {
      printf("%2d ", table.at(j).at(i).first);
      printf("%c ", table.at(j).at(i).second);
    }
    std::cout << '\n';
  }
  
  // trace path back to origin and visualise
  std::string out_top, out_mid, out_bot;
  for (unsigned col = cols - 1, row = rows - 1; col + row != 0; ) {
    //std::cout << row << ',' << col << ',' << table.at(row).at(col).second << '\n';
    if (row == 0 || table.at(row).at(col).second == 'L') {
      std::cout << "start\n";
      std::cout << "1\n";
      std::cout << row << '\n';
      std::cout << col << '\n';
      out_top.append(1, seq1.at(col - 1));
      out_bot.append(1, ' ');
      out_mid.append(1, ' ');
      --col;
    } else if (col == 0 || table.at(row).at(col).second == 'U') {
      out_top.append(1, ' ');
      out_bot.append(1, seq2.at(row - 1));
      out_mid.append(1, ' ');
      --row;
    } else if (table.at(row).at(col).second == 'D') {
      out_top.append(1, seq1.at(col - 1));
      out_bot.append(1, seq2.at(row - 1));
      if (seq1.at(col - 1) == seq2.at(row - 1)) {
        out_mid.append(1, '|');
      } else {
        out_mid.append(1, ' ');
      }
      --row;
      --col;
    }
    else {
      std::cout << "error\n";
      exit(EXIT_FAILURE);
    }
  }
  std::reverse(std::begin(out_top), std::end(out_top));
  std::reverse(std::begin(out_mid), std::end(out_mid));
  std::reverse(std::begin(out_bot), std::end(out_bot));
  std::cout << out_top << "\n";
  std::cout << out_mid << "\n";
  std::cout << out_bot << "\n";

	return table.at(rows - 1).at(cols - 1).first;
}

int main(int argc, char* argv[]) {
  if (argc < 2) std::cout << "please supply one input argument, the input file path\n";
  std::string input_file_name;
  input_file_name = argv[1];
  std::cout << "input_file_name = " << input_file_name << '\n';

  // check file exists
  std::ifstream instream(input_file_name);
  if (!instream.is_open()) {
    std::cout << "Could not open input file.\n";
    std::exit(EXIT_FAILURE);
  }

  // read file
  std::string seq1, seq2;
  std::getline(instream, seq1, '\n');
  std::getline(instream, seq2, '\n');
  std::transform(seq1.begin(), seq1.end(), seq1.begin(), ::toupper);
  std::transform(seq2.begin(), seq2.end(), seq2.begin(), ::toupper);
  instream.close();  
  
  // check valid sequence
  if (!ValidSequence(seq1)) {
    std::cout << "seq1 = " << seq1 << " is invalid\n";
  }
  if (!ValidSequence(seq2)) {
    std::cout << "seq2 = " << seq2 << " is invalid\n";
  }

  // print input strings
  std::cout << '\n';
  std::cout << "Input sequences:\n";
  std::cout << "seq1 = " << seq1 << '\n';
  std::cout << "seq2 = " << seq2 << '\n';
  std::cout << '\n';

  // align using dynamic programming
  DpAlignment(seq1, seq2);

  return 0;
}
