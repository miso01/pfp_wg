# compilation flags
CC=gcc
CFLAGS=-std=c99 -Wall -Wextra

CXX=g++
CXX_FLAGS=-std=c++11 -Wall -Wextra -DNDEBUG

EXECS=tfm_index_construct.x tfm_index_invert.x

.PHONY: build test clean release

build: $(EXECS)

test: build
	./tfm_index_construct.x -w 4 -p 50 data/yeast.raw
	./tfm_index_invert.x 	data/yeast.raw
	cmp 					data/yeast.raw data/yeast.raw.untunneled \
	&& echo "Output is correct."

clean:
	rm -f *.o *.x
	rm -f data/yeast.raw.*

# libs
%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

tfm_index_construct.x: tfm_index_construct.cpp gsacak.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lsdsl -ldivsufsort64 -lz -ldl -pthread -ldivsufsort

tfm_index_invert.x: tfm_index_invert.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lsdsl -ldivsufsort -ldivsufsort64
