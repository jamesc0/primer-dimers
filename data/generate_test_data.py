import random
import time

number_of_primers =200
primer_length = 25

outfile_name = "test_data_primers_{0}_{1}.txt".format(number_of_primers, primer_length)
base_dict = {0: "A", 1: "T", 2: "C", 3: "G"}

random.seed(time.time())
with open(outfile_name, "w") as outfile:
    outfile.write(str(number_of_primers))
    outfile.write("\n")
    for primer in range(number_of_primers):
        if primer != 0:
            outfile.write("\n")
        for base in range(primer_length):
            outfile.write(base_dict[random.randint(0,3)])
