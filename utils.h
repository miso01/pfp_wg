#include <stdint.h>
#include <stdio.h>

// special symbols used by the construction algorithm:
//   they cannot appear in the input file
//   the 0 symbol is used in the final BWT file as the EOF char

#define Dollar 2
#define EndOfWord 1 // word delimiter for the plain dictionary file
#define EndOfDict 0 // end of dictionary delimiter

// number of bytes to represent integers larger than 32 bit
#define IBYTES 5 // bytes used to represent a large integer (at most 8)
#define SABYTES IBYTES // bytes used to write a suffix array value in the output .sa file

// multi segment file functions
typedef struct multiFile {
    FILE *f;    // file currently opened
    char *base; // basename
    char *ext;  // file extension
    int cur;    // current index
    int nsegs;  // number of segments
} mFile;
