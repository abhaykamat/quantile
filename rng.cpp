
#include <rng.h>

#define M           397
#define MATRIX_A    0x9908b0df
#define UPPER_MASK  0x80000000
#define LOWER_MASK  0x7fffffff
static uint32_t mag01[2]= { 0x0, MATRIX_A };

void UniformRNG::reset(
    uint32_t seed
)
{
    if(seed==0) seed = 5489;
    states[0] = (seed & 0xFFFFFFFF);

    for(state=1; state<N; ++state) {
        uint32_t s = states[state-1];
        s = (1812433253*(s^(s>>30))+state);
        states[state] = (s & 0xFFFFFFFF);
    }
}

void UniformRNG::refill()
{
    int32_t k;
    uint32_t y;
    for(k=0; k<N-M; ++k) {
        y = (states[k]&UPPER_MASK) | (states[k+1]&LOWER_MASK);
        states[k] = states[k+M] ^ (y>>1) ^ mag01[y & 0x1];
    }

    while(k<N-1) {
        y = (states[k]&UPPER_MASK) | (states[k+1]&LOWER_MASK);
        states[k] = states[k+(M-N)] ^ (y>>1) ^ mag01[y & 0x1];
        ++k;
    }

    y = (states[N-1]&UPPER_MASK) | (states[0]&LOWER_MASK);
    states[N-1] = states[M-1] ^ (y>>1) ^ mag01[y & 0x1];
    state = 0;
}

