# compilation flags
CC=gcc
CFLAGS=-std=c99 -Wall -Wextra -g

CXX=g++
CXX_FLAGS=-std=c++11 -Wall -Wextra -g

EXECS=tfm_index_construct.x tfm_index_invert.x

.PHONY: build test clean release

build: $(EXECS)

test: build
	./tfm_index_construct.x -w 4 -p 50 -i data/yeast.raw -o data/yeast.wg
	./tfm_index_invert.x data/yeast.raw
	cmp data/yeast.raw.untunneled data/yeast.raw && echo "Output is correct."

clean:
	rm -f data/yeast.raw.* data/yeast.wg* *.x

# g++ -std=c++11 -Wall -Wextra -o tfm_index_construct.x tfm_index_construct.cpp gsacak.o utils.o -lsdsl
tfm_index_construct.x: tfm_index_construct.cpp
	$(CXX) $(CXX_FLAGS) -DNDEBUG -o $@ $^ -lsdsl

# g++ -std=c++11 -Wall -Wextra -o tfm_index_invert.x tfm_index_invert.cpp -lsdsl
tfm_index_invert.x: tfm_index_invert.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lsdsl
