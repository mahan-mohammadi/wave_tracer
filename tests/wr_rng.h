/* Deterministic xorshift64* generator for reproducible randomized tests. */
#ifndef WR_RNG_H
#define WR_RNG_H
#include <stdint.h>

typedef struct { uint64_t s; } wr_rng;

static inline uint64_t wr_rng_next(wr_rng *r)
{
    uint64_t x = r->s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    r->s = x;
    return x * 0x2545F4914F6CDD1DULL;
}
/* uniform in [0,1) with 53 random bits */
static inline double wr_rng_uniform(wr_rng *r) { return (double)(wr_rng_next(r) >> 11) * (1.0 / 9007199254740992.0); }
static inline double wr_rng_range(wr_rng *r, double a, double b) { return a + (b - a) * wr_rng_uniform(r); }
#endif
