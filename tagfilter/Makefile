CC=g++
CPPFLAGS=-O2 -msse4 -msse3 -I /home/xinyin/include/
LDFLAGS=-L /home/xinyin/lib/ -lpcrecpp
EXE_NAME=tagfilter

srcs = $(wildcard *.cpp)
objs = $(srcs:.cpp=.o)

$(EXE_NAME) : $(objs)
	$(CC) $(CPPFLAGS) -o $(EXE_NAME) $(objs) $(LDFLAGS)

.PHONY : clean

clean: 
	rm -f *.o $(EXE_NAME)
