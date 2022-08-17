# compilation flags
CC=/usr/bin/cc
CFLAGS=-O3 -Wall -std=c99 -g

CXX=/usr/bin/c++
CXX_FLAGS=-std=c++11 -Wall -Wextra -DNDEBUG
CXX_OPT_FLAGS=-O3 -ffast-math -funroll-loops -msse4.2 -march=native -DHAVE_CXA_DEMANGLE

EXECS=pfwg.x newscanNT.x tfm_index_construct.x tfm_index_invert.x

.PHONY: build test clean release

build: $(EXECS)

test: build
	./newscanNT.x -w 4 -p 50	data/yeast.raw
	./tfm_index_construct.x 	data/yeast.raw.parse data/yeast.raw.tunnel
	./pfwg.x -w 4 				data/yeast.raw
	./tfm_index_invert.x 		data/yeast.raw
	cmp 						data/yeast.raw data/yeast.raw.untunneled \
	&& echo "Output is correct."

clean:
	rm -f *.o *.x
	rm -f data/yeast.raw.*

# libs
%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

# executables
newscanNT.x: newscan.cpp utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ \
		-lz -ldl -DNOTHREADS

pfwg.x: pfwg.cpp gsacak.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ \
		-pthread -ldl -lsdsl -ldivsufsort -ldivsufsort64

tfm_index_construct.x: tfm_index_construct.cpp gsacak.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lsdsl -ldivsufsort -ldivsufsort64

tfm_index_invert.x: tfm_index_invert.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lsdsl -ldivsufsort -ldivsufsort64
