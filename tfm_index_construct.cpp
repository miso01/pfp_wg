#include <algorithm>
#include <assert.h>
#include <ctime>
#include <deque>
#include <errno.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sdsl/util.hpp>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>
#include <vector>
#include <zlib.h>

#include "tfm_index.hpp"

extern "C" {
#include "gsacak.h"
#include "utils.h"
}

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace sdsl;
using namespace __gnu_cxx;

// =============== algorithm limits ===================
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX - 1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;

// clang-format off
uint8_t asc2dnacat[] = {
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0,
           /*                                        -     */
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
           /*    A  B  C  D        G  H        K     M  N  */
    /*  80 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
           /*       R  S  T     V  W  X  Y */
    /*  96 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
           /*    a  b  c  d        g  h        k     m  n  */
    /* 112 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
           /*       r  s  t     v  w  x  y                 */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
// clang-format on

// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
    string str;
    occ_int_t occ;
    word_int_t rank = 0;
};

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
    string inputFileName = "";
    size_t w = 10;         // sliding window size and its default
    size_t p = 100;        // modulus for establishing stopping w-tuples
    bool SAinfo = false;   // compute SA information
    bool is_fasta = false; // read a fasta file
    bool compress = false; // parsing called in compress mode
    int th = 0;            // number of helper threads
    int verbose = 0;       // verbosity level
};

// -----------------------------------------------------------------
// class to maintain a window in a string and its KR fingerprint
struct KR_window {
    int wsize;
    int *window;
    int asize;
    const uint64_t prime = 1999999973;
    uint64_t hash;
    uint64_t tot_char;
    uint64_t asize_pot; // asize^(wsize-1) mod prime

    KR_window(int w) : wsize(w) {
        asize = 256;
        asize_pot = 1;
        for (int i = 1; i < wsize; i++)
            asize_pot =
                (asize_pot * asize) % prime; // ugly linear-time power algorithm
        // alloc and clear window
        window = new int[wsize];
        reset();
    }

    // init window, hash, and tot_char
    void reset() {
        for (int i = 0; i < wsize; i++)
            window[i] = 0;
        // init hash value and related values
        hash = tot_char = 0;
    }

    uint64_t addchar(int c) {
        int k = tot_char++ % wsize;
        // complex expression to avoid negative numbers
        hash +=
            (prime - (window[k] * asize_pot) % prime
            );                             // remove window[k] contribution
        hash = (asize * hash + c) % prime; //  add char i
        window[k] = c;
        // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
        return hash;
    }
    // debug only
    string get_window() {
        string w = "";
        int k = (tot_char - 1) % wsize;
        for (int i = k + 1; i < k + 1 + wsize; i++)
            w.append(1, window[i % wsize]);
        return w;
    }

    ~KR_window() { delete[] window; }
};
// -----------------------------------------------------------

static void save_update_word(
    string &w, unsigned int minsize, map<uint64_t, word_stats> &freq,
    FILE *tmp_parse_file, FILE *last, FILE *sa, uint64_t &pos
);

// compute 64-bit KR hash of a string
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
uint64_t kr_hash(string s) {
    uint64_t hash = 0;
    // const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
    const uint64_t prime =
        27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for (size_t k = 0; k < s.size(); k++) {
        int c = (unsigned char)s[k];
        assert(c >= 0 && c < 256);
        hash = (256 * hash + c) % prime; //  add char k
    }
    return hash;
}

// save current word in the freq map and update it leaving only the
// last minsize chars which is the overlap with next word
static void save_update_word(
    string &w, unsigned int minsize, map<uint64_t, word_stats> &freq,
    FILE *tmp_parse_file, FILE *last, FILE *sa, uint64_t &pos
) {
    assert(pos == 0 || w.size() > minsize);
    if (w.size() <= minsize)
        return;
    // get the hash value and write it to the temporary parse file
    uint64_t hash = kr_hash(w);
    if (fwrite(&hash, sizeof(hash), 1, tmp_parse_file) != 1)
        die("parse write error");

    // update frequency table for current hash
    if (freq.find(hash) == freq.end()) {
        freq[hash].occ = 1; // new hash
        freq[hash].str = w;
    } else {
        freq[hash].occ += 1; // known hash
        if (freq[hash].occ <= 0) {
            cerr << "Emergency exit! Maximum # of occurence of dictionary word "
                    "(";
            cerr << MAX_WORD_OCC << ") exceeded\n";
            exit(1);
        }
        if (freq[hash].str != w) {
            cerr << "Emergency exit! Hash collision for strings:\n";
            cerr << freq[hash].str << "\n  vs\n" << w << endl;
            exit(1);
        }
    }
    // output char w+1 from the end
    if (fputc(w[w.size() - minsize - 1], last) == EOF)
        die("Error writing to .last file");
    // compute ending position +1 of current word and write it to sa file
    // pos is the ending position+1 of the previous word and is updated here
    if (pos == 0)
        pos = w.size() - 1; // -1 is for the initial $ of the first word
    else
        pos += w.size() - minsize;
    if (sa)
        if (fwrite(&pos, IBYTES, 1, sa) != 1)
            die("Error writing to sa info file");
    // keep only the overlapping part of the window
    w.erase(0, w.size() - minsize);
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is immediately written to the parse file
uint64_t process_file(Args &arg, map<uint64_t, word_stats> &wordFreq) {
    // open a, possibly compressed, input file
    string fnam = arg.inputFileName;
    // open the 1st pass parsing file
    FILE *g = open_aux_file(arg.inputFileName.c_str(), EXTPARS0, "wb");
    // open output file containing the char at position -(w+1) of each word
    FILE *last_file = open_aux_file(arg.inputFileName.c_str(), EXTLST, "wb");
    // if requested open file containing the ending position+1 of each word
    FILE *sa_file = NULL;
    if (arg.SAinfo)
        sa_file = open_aux_file(arg.inputFileName.c_str(), EXTSAI, "wb");

    // main loop on the chars of the input file
    int c;
    uint64_t pos = 0; // ending position +1 of previous word in the original
                      // text, used for computing sa_info
    assert(
        IBYTES <= sizeof(pos)
    ); // IBYTES bytes of pos are written to the sa info file
    // init first word in the parsing with a Dollar char
    string word("");
    word.append(1, Dollar);
    KR_window krw(arg.w);
    std::string line;
    if (arg.is_fasta) {
        gzFile fp;
        kseq_t *seq;
        int l;
        fp = gzopen(fnam.c_str(), "r");
        seq = kseq_init(fp);
        while ((l = kseq_read(seq)) >= 0) {
            for (size_t i = 0; i < seq->seq.l; i++) {
                c = std::toupper(seq->seq.s[i]);
                if (c <= Dollar) {
                    cerr << "Invalid char found in input file: no additional "
                            "chars will be read\n";
                    break;
                }
                word.append(1, c);
                uint64_t hash = krw.addchar(c);
                if (hash % arg.p == 0) {
                    save_update_word(
                        word, arg.w, wordFreq, g, last_file, sa_file, pos
                    );
                }
            }
            if (c <= Dollar)
                break;
        }
        kseq_destroy(seq);
        gzclose(fp);
    } else {
#ifdef GZSTREAM
        igzstream f(fnam.c_str());
#else
        ifstream f(fnam);
#endif
        if (!f.rdbuf()->is_open()) { // is_open does not work on igzstreams
            perror(__func__);
            throw new std::runtime_error("Cannot open input file " + fnam);
        }
        while ((c = f.get()) != EOF) {
            if (c <= Dollar) {
                cerr << "Invalid char found in input file: no additional chars "
                        "will be read\n";
                break;
            }
            word.append(1, c);
            uint64_t hash = krw.addchar(c);
            if (hash % arg.p == 0) {
                // end of word, save it and write its full hash to the output
                // file cerr << "~"<< c << "~ " << hash << " ~~ <" << word << ">
                // ~~ <" << krw.get_window() << ">" <<  endl;
                save_update_word(
                    word, arg.w, wordFreq, g, last_file, sa_file, pos
                );
            }
        }
        f.close();
    }
    // virtually add w null chars at the end of the file and add the last word
    // in the dict
    word.append(arg.w, Dollar);
    save_update_word(word, arg.w, wordFreq, g, last_file, sa_file, pos);
    // close input and output files
    if (sa_file)
        if (fclose(sa_file) != 0)
            die("Error closing SA file");
    if (fclose(last_file) != 0)
        die("Error closing last file");
    if (fclose(g) != 0)
        die("Error closing parse file");
    if (pos != krw.tot_char + arg.w)
        cerr << "Pos: " << pos << " tot " << krw.tot_char << endl;
    return krw.tot_char;
}

// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b) { return *a <= *b; }

// given the sorted dictionary and the frequency map write the dictionary and
// occ files also compute the 1-based rank for each hash
void writeDictOcc(
    Args &arg, map<uint64_t, word_stats> &wfreq,
    vector<const string *> &sortedDict
) {
    assert(sortedDict.size() == wfreq.size());
    FILE *fdict;
    // open dictionary and occ files
    if (arg.compress)
        fdict = open_aux_file(arg.inputFileName.c_str(), EXTDICZ, "wb");
    else
        fdict = open_aux_file(arg.inputFileName.c_str(), EXTDICT, "wb");
    FILE *focc = open_aux_file(arg.inputFileName.c_str(), EXTOCC, "wb");

    word_int_t wrank = 1; // current word rank (1 based)
    for (auto x : sortedDict) {
        const char *word = (*x).data(); // current dictionary word
        int offset = 0;
        size_t len = (*x).size(); // offset and length of word
        assert(len > (size_t)arg.w);
        if (arg.compress) { // if we are compressing remove overlapping and
                            // extraneous chars
            len -= arg.w;   // remove the last w chars
            if (word[0] == Dollar) {
                offset = 1;
                len -= 1;
            } // remove the very first Dollar
        }
        size_t s = fwrite(word + offset, 1, len, fdict);
        if (s != len)
            die("Error writing to DICT file");
        if (fputc(EndOfWord, fdict) == EOF)
            die("Error writing EndOfWord to DICT file");
        uint64_t hash = kr_hash(*x);
        auto &wf = wfreq.at(hash);
        assert(wf.occ > 0);
        s = fwrite(&wf.occ, sizeof(wf.occ), 1, focc);
        if (s != 1)
            die("Error writing to OCC file");
        assert(wf.rank == 0);
        wf.rank = wrank++;
    }
    if (fputc(EndOfDict, fdict) == EOF)
        die("Error writing EndOfDict to DICT file");
    if (fclose(focc) != 0)
        die("Error closing OCC file");
    if (fclose(fdict) != 0)
        die("Error closing DICT file");
}

void remapParse(Args &arg, map<uint64_t, word_stats> &wfreq) {
    // open parse files. the old parse can be stored in a single file or in
    // multiple files
    mFile *moldp = mopen_aux_file(arg.inputFileName.c_str(), EXTPARS0, arg.th);
    FILE *newp = open_aux_file(arg.inputFileName.c_str(), EXTPARSE, "wb");

    // recompute occ as an extra check
    vector<occ_int_t> occ(wfreq.size() + 1, 0); // ranks are zero based
    uint64_t hash;
    while (true) {
        size_t s = mfread(&hash, sizeof(hash), 1, moldp);
        if (s == 0)
            break;
        if (s != 1)
            die("Unexpected parse EOF");
        word_int_t rank = wfreq.at(hash).rank;
        occ[rank]++;
        s = fwrite(&rank, sizeof(rank), 1, newp);
        if (s != 1)
            die("Error writing to new parse file");
    }
    if (fclose(newp) != 0)
        die("Error closing new parse file");
    if (mfclose(moldp) != 0)
        die("Error closing old parse segment");
    // check old and recomputed occ coincide
    // for(auto& x : wfreq)
    //  assert(x.second.occ == occ[x.second.rank]);
}

void print_help(char **argv, Args &args) {
    cout << "Usage: " << argv[0] << " <input filename> [options]" << endl;
    cout << "  Options: " << endl
         << "\t-w W\tsliding window size, def. " << args.w << endl
         << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
         << "\t-h  \tshow help and exit" << endl
         << "\t-s  \tcompute suffix array info" << endl;
#ifdef GZSTREAM
    cout << "If the input file is gzipped it is automatically extracted\n";
#endif
    exit(1);
}

void parseArgs(int argc, char **argv, Args &arg) {
    int c;
    extern char *optarg;
    extern int optind;

    puts("==== Command line:");
    for (int i = 0; i < argc; i++)
        printf(" %s", argv[i]);
    puts("");

    string sarg;
    while ((c = getopt(argc, argv, "p:w:fsht:v")) != -1) {
        switch (c) {
        case 's':
            arg.SAinfo = true;
            break;
        case 'c':
            arg.compress = true;
            break;
        case 'w':
            sarg.assign(optarg);
            arg.w = stoi(sarg);
            break;
        case 'p':
            sarg.assign(optarg);
            arg.p = stoi(sarg);
            break;
        case 'f':
            arg.is_fasta = true;
            break;
        case 't':
            sarg.assign(optarg);
            arg.th = stoi(sarg);
            break;
        case 'v':
            arg.verbose++;
            break;
        case 'h':
            print_help(argv, arg);
            exit(1);
        case '?':
            cout << "Unknown option. Use -h for help." << endl;
            exit(1);
        }
    }
    // the only input parameter is the file name
    if (argc == optind + 1) {
        arg.inputFileName.assign(argv[optind]);
    } else {
        cout << "Invalid number of arguments" << endl;
        print_help(argv, arg);
    }
    // check algorithm parameters
    if (arg.w < 4) {
        cout << "Windows size must be at least 4\n";
        exit(1);
    }
    if (arg.p < 10) {
        cout << "Modulus must be at leas 10\n";
        exit(1);
    }
    if (arg.th != 0) {
        cout << "The NT version cannot use threads\n";
        exit(1);
    }
}

bool is_gzipped(std::string fname) {
    FILE *fp = fopen(fname.c_str(), "rb");
    int byte1 = 0, byte2 = 0;
    fread(&byte1, sizeof(char), 1, fp);
    fread(&byte2, sizeof(char), 1, fp);
    fclose(fp);
    return (byte1 == 0x1f && byte2 == 0x8b);
}

// -------------------------------------------------------
// type used to represent an entry in the SA
// this is currently 32 bit for gsacak and 64 bit for gsacak-64
// note that here we use sacak (SA computation for a single string of 32 bit
// symbols)
typedef uint_t sa_index_t;

void printUsage(char **argv) {
    cerr << "USAGE: " << argv[0] << " INFILE TFMOUTFILE" << endl;
    cerr << "INFILE:" << endl;
    cerr << "  Parse file to construct tunneled FM index from, nullbytes are "
            "permitted"
         << endl;
    cerr << "TFMOUTFILE:" << endl;
    cerr
        << "  File where to store the serialized tunneled FM index of the parse"
        << endl;
};

// little hack to get extra information from memory managemant
const format_type leet_format = (format_type)1337;
namespace sdsl {
template <>
void write_mem_log<leet_format>(ostream &out, const memory_monitor &m) {

    // get all memory events
    auto events = m.completed_events;
    std::sort(events.begin(), events.end());
    auto e = events.begin();

    // scan events and detect peak and time of both suffix array and xbwt
    // construction

    // scan for fm index construction
    int64_t fm_mem_peak = 0;
    auto fm_start = m.start_log;
    auto fm_end = fm_start;
    while (e != events.end() && e->name != "FINDMINDBG") {
        for (auto alloc : e->allocations) {
            fm_mem_peak = (std::max)(fm_mem_peak, alloc.usage);
            fm_end = alloc.timestamp;
        }
        ++e;
    }

    // scan for min dbg algorithm
    int64_t dbg_mem_peak = 0;
    auto dbg_start = fm_end;
    auto dbg_end = dbg_start;
    while (e != events.end() && e->name != "TFMINDEXCONSTRUCT") {
        for (auto alloc : e->allocations) {
            dbg_mem_peak = (std::max)(dbg_mem_peak, alloc.usage);
            dbg_end = alloc.timestamp;
        }
        ++e;
    }

    // scan for tunneled fm index construction
    int64_t tfm_mem_peak = 0;
    auto tfm_start = dbg_end;
    auto tfm_end = tfm_start;
    while (e != events.end()) {
        for (auto alloc : e->allocations) {
            tfm_mem_peak = (std::max)(tfm_mem_peak, alloc.usage);
            tfm_end = alloc.timestamp;
        }
        ++e;
    }

    // print results
    out << "fm_mem_peak\t" << fm_mem_peak << endl;
    out << "fm_time\t"
        << chrono::duration_cast<chrono::milliseconds>(fm_end - fm_start)
               .count()
        << endl;

    out << "min_dbg_mem_peak\t" << dbg_mem_peak << endl;
    out << "min_dbg_time\t"
        << chrono::duration_cast<chrono::milliseconds>(dbg_end - dbg_start)
               .count()
        << endl;

    out << "tfm_mem_peak\t" << tfm_mem_peak << endl;
    out << "tfm_time\t"
        << chrono::duration_cast<chrono::milliseconds>(tfm_end - tfm_start)
               .count()
        << endl;
};
}; // namespace sdsl

void compute_BWT(uint32_t *Text, long n, long k, string filename) {
    sa_index_t *SA = (sa_index_t *)malloc(n * sizeof(*SA));
    printf("Computing SA of size %ld over an alphabet of size %ld\n", n, k);
    int depth = sacak_int(Text, SA, n, k);
    printf("SA computed with depth: %d\n", depth);

    // transform SA->BWT inplace and write remapped last array, and possibly
    // sainfo
    sa_index_t *BWTsa = SA; // BWT overlapping SA
    assert(n > 1);
    // first BWT symbol
    assert(SA[0] == n);
    // 2nd, 3rd etc BWT symbols
    for (long i = 0; i < n; i++) {
        if (SA[i] == 0) {
            assert(i == 1); // Text[0]=$abc... is the second lex word
            BWTsa[i] = 0;   // eos in BWT, there is no phrase in D corresponding
                            // to this symbol so we write dummy values
        } else {
            BWTsa[i] = Text[SA[i] - 1];
        }
        // if(BWTsa[i]==0) cout << i << endl;
    }
    printf("BWT constructed\n");

    FILE *fout = fopen(filename.c_str(), "wb");
    fwrite(BWTsa, sizeof(BWTsa[0]), n, fout);
    fclose(fout);
    printf("BWT written to file\n");
}

uint32_t *load_parse(const string &infile, size_t &psize) {
    FILE *fp = fopen(infile.c_str(), "r");

    fseek(fp, 0, SEEK_END);
    size_t n;
    n = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    n = n / 4;
    psize = n;

    uint32_t *parse = new uint32_t[n + 1];
    size_t s = fread(parse, sizeof(*parse), n, fp);
    if (s != n) {
        printf("Parse loading error!");
        exit(1);
    }
    parse[n] = 0; // sacak needs this
    return parse;
}

size_t compute_sigma(const uint32_t *parse, const size_t psize) {
    uint32_t max = 0;
    for (size_t i = 0; i < psize; i++) {
        if (max < parse[i])
            max = parse[i];
    }
    return (max + 1);
}

void calculate_word_frequencies(Args &arg, map<uint64_t, word_stats> &wordFreq) {
    try {
        process_file(arg, wordFreq);
    } catch (const std::bad_alloc &) {
        cout << "Out of memory (parsing phase)... emergency exit\n";
        die("bad alloc exception");
    }

    
    if (wordFreq.size() > MAX_DISTINCT_WORDS) {
        cerr << "Emergency exit! The number of distinc words (" << wordFreq.size()
             << ")\n";
        cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS
             << ")\n";
        exit(1);
    }
}

int main(int argc, char **argv) {
    Args arg;
    parseArgs(argc, argv, arg);

    map<uint64_t, word_stats> wordFreq;
    calculate_word_frequencies(arg, wordFreq); // + <fn>.last <fn>.parse_old

    // create array of dictionary words
    vector<const string *> dictArray;
    uint64_t totDWord = wordFreq.size();
    dictArray.reserve(totDWord);
    for (auto &x : wordFreq) {
        dictArray.push_back(&x.second.str);
    }
    assert(dictArray.size() == totDWord);
    sort(dictArray.begin(), dictArray.end(), pstringCompare);
    // write plain dictionary and occ file, also compute rank for each hash
    writeDictOcc(arg, wordFreq, dictArray); // + <fn>.dict <fn>.occ
    dictArray.clear(); // reclaim memory

    remapParse(arg, wordFreq); // + <fn>.parse
    
    // construct tunneled fm index
    tfm_index<> tfm;

    cache_config config(true, "./", util::basename(arg.inputFileName));
    size_t psize;
    uint32_t *parse = load_parse(arg.inputFileName + ".parse", psize);
    size_t sigma = compute_sigma(parse, psize);
    compute_BWT(parse, psize + 1, sigma, arg.inputFileName + ".bwt"); // <fn>.bwt
    delete parse;

    construct_tfm_index(tfm, arg.inputFileName + ".bwt", psize + 1, config);
    store_to_file(tfm, arg.inputFileName + ".tunnel"); // <fn>.tunnel

    return 0;
}
