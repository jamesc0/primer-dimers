CC = g++

CPPFLAGS=-std=c++11 -Wall -lm

.PHONY : clean

%.o : %.cc %.h
	$(CC) $(CPPFLAGS) -c $<

clean :
	rm -f *.o a.out main jmer_counting
