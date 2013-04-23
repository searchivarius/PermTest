/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Leonid Boytsov, http://boytsov.info
 *
 * For details see:
 *   Leonid Boytsov, Anna Belova, Peter Westfall, 2013, 
 *   Deciding on an Adjustment for Multiplicity in IR Experiments.
 *   In Proceedings of SIGIR 2013.
 *
 */
#include <stdlib.h>

//#define STANDARD_RAND

void init_genrand(unsigned long s);
unsigned long genrand_int32(void);

inline void InitRand(unsigned long seed) {
#ifdef STANDARD_RAND
    srand(seed);
#else
    init_genrand(seed);
#endif
}

inline unsigned long GenRand() {
#ifdef STANDARD_RAND
    return rand();
#else
    return genrand_int32();
#endif
}
