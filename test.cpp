#include <cstdint>
#include <cstdlib>
#include <iostream>

extern "C" {
#include "gsacak.c"
}

using namespace std;

int main() {
    // 214414413310
    // mississippi$

    // 134420314411 -> 356642536633
    // ipssm$pissii
    uint32_t t[] = {2,1,4,4,1,4,4,1,3,3,1,0};

    uint32_t *in = (uint32_t *)malloc((12+2) * sizeof(*in));
    for (size_t i = 0; i < 12; i++) in[i] = t[i]+2;
    in[12] = 1;
    in[13] = 0;

    uint32_t *SA = (uint32_t *)malloc((12+2) * sizeof(*SA));

    gsacak_int(in, SA, NULL, NULL, 12+2, 5+2);

    for (size_t i = 0; i < 12 + 2; i++) {
        if (SA[i] == 0) { cout << 0 << " "; continue; }
        if (in[SA[i]-1] == 1) continue;
        if (in[SA[i]-1] == 2) continue;
        else cout << in[SA[i] - 1]-2 << " ";
    }
    cout << endl;
}