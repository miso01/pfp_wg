# compilation flags
CC=gcc
CFLAGS=-std=c99 -Wall -Wextra

CXX=g++
CXX_FLAGS=-std=c++11 -Wall -Wextra -DNDEBUG

EXECS=pfwg.x newscan.x tfm_index_construct.x tfm_index_invert.x

.PHONY: build test clean release

build: $(EXECS)

test: build
	./newscan.x -w 4 -p 50	data/yeast.raw
	./tfm_index_construct.x data/yeast.raw.parse data/yeast.raw.tunnel
	./pfwg.x -w 4 			data/yeast.raw
	./tfm_index_invert.x 	data/yeast.raw
	cmp 					data/yeast.raw data/yeast.raw.untunneled \
	&& echo "Output is correct."

clean:
	rm -f *.o *.x
	rm -f data/yeast.raw.*

# libs
%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

# executables
newscan.x: newscan.cpp utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ \
		-lz -ldl

pfwg.x: pfwg.cpp gsacak.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ \
		-pthread -ldl -lsdsl -ldivsufsort -ldivsufsort64

tfm_index_construct.x: tfm_index_construct.cpp gsacak.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lsdsl -ldivsufsort64

tfm_index_invert.x: tfm_index_invert.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lsdsl -ldivsufsort -ldivsufsort64
