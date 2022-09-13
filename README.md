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