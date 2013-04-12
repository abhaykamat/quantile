#ifndef __RNG_H__
    #define __RNG_H__

    #include <stddef.h>
    #include <stdint.h>

    struct RNG
    {
    private:

        enum { N = 624 };

        int32_t  state;
        uint32_t states[N];

        void refill();

        uint32_t mersenneSample32()
        {
            if(N<=state) {
                refill();
            }

            uint32_t y =  states[state++];
            y ^= (y >> 11);
            y ^= (y <<  7) & 0x9D2C5680;
            y ^= (y << 15) & 0xEFC60000;
            y ^= (y >> 18);
            return y;
        }

    public:

        void reset(uint32_t seed=0);

        RNG(uint32_t seed=0) { reset(seed); }
        ~RNG()               {              } 

        uint32_t sample32()
        {
            return mersenneSample32();
        }

        uint32_t stronger_sample32()
        {
            uint32_t s = mersenneSample32();
            return (s>>13) ^ mersenneSample32();
        }

        uint64_t sample64()
        {
            uint64_t s = sample32();
            return (s<<32) + sample32();
        }

        double sample()
        {
            return sample64() * 0x1p-64;
        }

        double sample(
            double  min,
            double  max
        )
        {
            return min + sample()*(max-min);
        }
    };

#endif // __RNG_H__

