CC=g++
CPPFLAGS=-O2 -msse4 -std=c++0x 

tagcluster: 
	$(CC) $(CPPFLAGS) -o tagcluster tagcluster.cpp iodata.c

.PHONY : clean

clean: 
	rm -f *.o tagcluster
