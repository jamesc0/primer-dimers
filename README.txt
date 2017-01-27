An algorithm which, given a list of primers, returns primer dimer candidates.

Compile using

make main

Run using

./main data/data.txt > out.txt

replacing data/data.txt with the name of the input file.
The input file should be formatted like the file at data/data.txt
The parameters can be changed in lines 12-17 of the main.cc code

The algorithm take 5.6 seconds to run on 200 candidate primers.
The running time should grow quadratically with the number of candidate primers.
