#ifndef __RNG_H__
    #define __RNG_H__

    #include <math.h>
    #include <stdint.h>
    #include <stdint.h>

    struct RNG
    {
        virtual void reset(uint32_t seed) = 0;
        virtual double sample() = 0;
    };

    struct UniformRNG: public RNG
    {
    private:

        int32_t  state;
        enum { N = 624 };
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

        virtual void reset(
            uint32_t seed = 0
        );

        UniformRNG(
            uint32_t seed = 0
        )
        {
            reset(seed);
        }

        virtual ~UniformRNG()
        {
        } 

        uint32_t sample32()
        {
            return mersenneSample32();
        }

        uint32_t strongSample32()
        {
            uint32_t s0 = mersenneSample32();
            uint32_t s1 = mersenneSample32();
            uint32_t s2 = mersenneSample32();
            return (s0<<13) ^ (s1) ^ (s2>>13);
        }

        uint64_t sample64()
        {
            uint64_t s = sample32();
            return (s<<32) | sample32();
        }

        virtual double sample()
        {
            return sample64() * 0x1p-64;
        }

        double sample(
            double  min,
            double  max
        )
        {
            return min + sample()*(max - min);
        }
    };

    // Normal distribution
    struct GaussRNG: public RNG
    {
    private:
        UniformRNG urng0;
        UniformRNG urng1;
        double stashedValue;
        bool haveStashedValue;

    public:

        virtual void reset(
            uint32_t seed = 0
        )
        {
            urng0.reset(3475*seed + 0x123);
            urng1.reset(7543*seed + 0x321);
        }

        GaussRNG(
            uint32_t seed = 0
        )
            :   stashedValue(0.0),
                haveStashedValue(false)
        {
            reset(seed);
        }

        virtual ~GaussRNG()
        {
        }

        // Marsaglia polar method for the normal distribution
        double sample()
        {
            if(haveStashedValue) {
                haveStashedValue = false;
                return stashedValue;
            }

            while(1) {
                double x0  = 2.0*urng0.sample() - 1.0; 
                double y0  = 2.0*urng1.sample() - 1.0;
                double r2 = x0*x0 + y0*y0;
                if(r2<=1.0) {
                    double scaleFactor = ::sqrt(-2.0*log(r2)/r2);
                    stashedValue = y0*scaleFactor;
                    haveStashedValue = true;
                    return x0*scaleFactor;
                }
            }
        }
    };

    // Positive half of the normal distribution
    struct HalfGaussRNG: public RNG
    {
    private:
        GaussRNG gaussRNG;

    public:

        virtual void reset(
            uint32_t seed = 0
        )
        {
            gaussRNG.reset(seed);
        }

        HalfGaussRNG(
            uint32_t seed = 0
        )
        {
            reset(seed);
        }

        virtual ~HalfGaussRNG()
        {
        }

        virtual double sample()
        {
            double x = gaussRNG.sample();
            return fabs(x);
        }
    };

#endif // __RNG_H__

