#include <iostream>
#include <map>
#include <deque>
#include <utility>
#include <vector>

#include <iostream>
#include <fstream>

#include <sdsl/util.hpp>

#include "tfm_index.hpp"

extern "C" {
#include "gsacak.h"
#include "utils.h"
}

using namespace std;
using namespace sdsl;

// -------------------------------------------------------
// type used to represent an entry in the SA
// this is currently 32 bit for gsacak and 64 bit for gsacak-64
// note that here we use sacak (SA computation for a single string of 32 bit symbols)
typedef uint_t sa_index_t;

void printUsage(char **argv) {
    cerr << "USAGE: " << argv[0] << " INFILE TFMOUTFILE" << endl;
    cerr << "INFILE:" << endl;
    cerr << "  Parse file to construct tunneled FM index from, nullbytes are permitted" << endl;
    cerr << "TFMOUTFILE:" << endl;
    cerr << "  File where to store the serialized tunneled FM index of the parse" << endl;
};

//little hack to get extra information from memory managemant
const format_type leet_format = (format_type) 1337;
namespace sdsl {
    template<>
    void write_mem_log<leet_format>(ostream &out, const memory_monitor &m) {

        //get all memory events
        auto events = m.completed_events;
        std::sort(events.begin(), events.end());
        auto e = events.begin();

        //scan events and detect peak and time of both suffix array and xbwt construction

        //scan for fm index construction
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

        //scan for min dbg algorithm
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

        //scan for tunneled fm index construction
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

        //print results
        out << "fm_mem_peak\t" << fm_mem_peak << endl;
        out << "fm_time\t" << chrono::duration_cast<chrono::milliseconds>(fm_end - fm_start).count() << endl;

        out << "min_dbg_mem_peak\t" << dbg_mem_peak << endl;
        out << "min_dbg_time\t" << chrono::duration_cast<chrono::milliseconds>(dbg_end - dbg_start).count() << endl;

        out << "tfm_mem_peak\t" << tfm_mem_peak << endl;
        out << "tfm_time\t" << chrono::duration_cast<chrono::milliseconds>(tfm_end - tfm_start).count() << endl;
    };
};

void compute_BWT(uint32_t *Text, long n, long k, string filename) {
    sa_index_t *SA = (sa_index_t*)malloc(n*sizeof(*SA));
    printf("Computing SA of size %ld over an alphabet of size %ld\n", n, k);
    int depth = sacak_int(Text, SA, n, k);
    printf("SA computed with depth: %d\n", depth);

    // transform SA->BWT inplace and write remapped last array, and possibly sainfo
    sa_index_t *BWTsa = SA; // BWT overlapping SA
    assert(n>1);
    // first BWT symbol
    assert(SA[0] == n);
    // 2nd, 3rd etc BWT symbols
    for (long i = 0; i < n; i++) {
        if(SA[i] == 0) {
            assert(i == 1);  // Text[0]=$abc... is the second lex word
            BWTsa[i] = 0;   // eos in BWT, there is no phrase in D corresponding to this symbol so we write dummy values
        }
        else {
            BWTsa[i] = Text[SA[i]-1];
        }
        //if(BWTsa[i]==0) cout << i << endl;
    }
    printf("BWT constructed\n");

    FILE *fout = fopen(filename.c_str(), "wb");
    fwrite(BWTsa, sizeof(BWTsa[0]), n, fout);
    fclose(fout);
    printf("BWT written to file\n");
}

uint32_t* load_parse(const string &infile, size_t &psize) {
    FILE *fp = fopen(infile.c_str(), "r");

    fseek(fp, 0, SEEK_END);
    size_t n;
    n = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    n = n/4;
    psize = n;

    uint32_t *parse = new uint32_t[n+1];
    size_t s = fread(parse, sizeof(*parse), n, fp);
    if (s != n) { printf("Parse loading error!"); exit(1); }
    parse[n] = 0; //sacak needs this
    return parse;
}


size_t compute_sigma(const uint32_t* parse, const size_t psize) {
    uint32_t max = 0;
    for (size_t i = 0; i < psize; i++) {
        if (max < parse[i]) max = parse[i];
    }
    return (max + 1);
}


int main(int argc, char **argv) {
    if (argc < 3) { printUsage(argv); return 1; }

    string infile = argv[argc - 2];
    string outfile = argv[argc - 1];

    //construct tunneled fm index
    tfm_index<> tfm;

    cache_config config(true, "./", util::basename(infile));
    size_t psize;
    uint32_t *parse = load_parse(infile, psize);
    size_t sigma = compute_sigma(parse, psize);
    compute_BWT(parse, psize+1, sigma, infile + ".bwt");
    delete parse;
    construct_tfm_index(tfm, infile + ".bwt", psize+1, config);

    store_to_file(tfm, outfile);
    return 0;
}
