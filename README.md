Original work:                              \
https://gitlab.com/manzai/Big-BWT           \
https://github.com/simongog/sdsl-lite       \
https://github.com/waYne1337/BWT-Tunneling  \
https://github.com/felipelouza/gsa-is       \

# How to install and run
1. Install system-wide SDSL (https://github.com/simongog/sdsl-lite)
2. Run `make build test clean`

# Used for code analysis
```
cppcheck --enable=unusedFunction tfm_index_construct.cpp
```

# Profiling
https://youtu.be/fDlE93hs_-U
```
g++ -std=c++11 -O3 -march=native -g myprog.cpp -o myprog

# memory leaks
valgrind ./myprog

# cache and branch prediction
valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes --cachegrind-out-file=chg.out ./tfm_index_construct.x -w 2 -p 11 -i data/yeast.small -o data/yeast.wg

cg_annotate chg.out $PWD/tfm_index_construct.cpp

# system calls
strace -c ./tfm_index_construct.x -w 2 -p 11 -i data/yeast.small -o data/yeast.wg

# on which line the time is spent
# sudo apt install linux-tools-common linux-tools-generic
sudo perf stat ./tfm_index_construct.x -w 2 -p 11 -i data/yeast.small -o data/yeast.wg
sudo perf record ./tfm_index_construct.x -w 2 -p 11 -i data/yeast.small -o data/yeast.wg
sudo perf report
sudo perf annotate

# what is the difference between callgrind and cachegrind?
valgrind --tool=callgrind --branch-sim=yes --cacheuse=yes --callgrind-out-file=clg.out ./tfm_index_construct.x -w 2 -p 11 -i data/yeast.small -o data/yeast.wg

callgrind_annotate clg.out tfm_index_construct.cpp | less

kcachegrind clg.out
```