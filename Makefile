# compilation flags
CC=gcc
CFLAGS=-std=c99 -Wall -Wextra -g

CXX=g++
CXX_FLAGS=-std=c++11 -Wall -Wextra -g

EXECS=tfm_index_construct.x tfm_index_invert.x

.PHONY: build test clean release small_test

build: $(EXECS)

test: build
	./tfm_index_construct.x -w 4 -p 50 -i data/yeast.raw -o data/yeast.wg
	./tfm_index_invert.x data/yeast.wg data/yeast.raw.untunneled
	cmp data/yeast.raw.untunneled data/yeast.raw && echo "Output is correct."

small_test: build
	./tfm_index_construct.x -w 4 -p 50 -i data/yeast.small -o data/yeast.wg
	./tfm_index_invert.x data/yeast.wg data/yeast.small.untunneled
	cmp data/yeast.small.untunneled data/yeast.small && echo "Output is correct."

clean:
	rm -f data/yeast.raw.* data/yeast.wg* *.x data/yeast.small.*

tfm_index_construct.x: tfm_index_construct.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lsdsl

tfm_index_invert.x: tfm_index_invert.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lsdsl
