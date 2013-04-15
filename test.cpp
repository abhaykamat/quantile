
#include <rng.h>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdint.h>
#include <mpreal.h>
#include <boost/math/distributions/normal.hpp>
using namespace mpfr;

// Get us a 256-bit quantile implementation
namespace MyQuantile {
    static const int kBits(256);
    #include <quantile.h>
}

// Run a simulation of a half-normal distribution, count how many samples "fit"
static double sampleProportion(
    double variance,
    double threshold,
    double nbSamples = 10 * 1000 * 1000
)
{
    int count = 0;
    HalfGaussRNG hgrng;
    for(int i=0; i<nbSamples; ++i) {
        double x = variance * hgrng.sample();
        count += (x<threshold ? 1 : 0);
    }
    return count / (double)nbSamples;
}

// Compute variance of a half-normal distribution so that "wantedProportion" of samples below "threshold"
static double findVariance(
    double threshold,
    double proportion
)
{
    double y = 0.5*(proportion + 1.0);
    return threshold / MyQuantile::quantile(y).toDouble();
}

// Test if findVariance works
static void testFindVariance(
    double wantedProportion,
    double threshold
)
{
    // We want "wantedProportion" of samples below "threshold"
    double computedVariance = findVariance(threshold, wantedProportion);

    // See if that works as advertized
    double sampledProportion = sampleProportion(
        computedVariance,
        threshold
    );

    printf(
        "threshold=%14.8f wantedProportion = %14.8f, computedVariance=%25.20f. Simulation shows: sampledProportion = %14.8f\n",
        threshold,
        wantedProportion,
        computedVariance,
        sampledProportion
    );
}

static void t0()
{
    testFindVariance(0.07123,  0.04321); // want  7.123 % of samples below  0.04321
    testFindVariance(0.50000, 10.00000); // want 50.000 % of samples below 10.00000
    testFindVariance(0.01000,  1.00000); // want  1.000 % of samples below  1.00000
    testFindVariance(0.90000,  1.00000); // want 90.000 % of samples below  1.00000
    testFindVariance(0.99000,  1.00000); // want 99.000 % of samples below  1.00000
    testFindVariance(0.99900,  1.00000); // want 99.900 % of samples below  1.00000
    testFindVariance(0.99990,  1.00000); // want 99.990 % of samples below  1.00000
    testFindVariance(0.99999,  1.00000); // want 99.999 % of samples below  1.00000
    printf("\n");
}

static void t1()
{
    static const int nbSamples = 20;
    const boost::math::normal dist(0.0, 1.0);

    for(int i=0; i<=nbSamples; ++i) {

        double mu = i/(double)nbSamples;

        double y0 = 1e-12;
        double y1 = 1.0 - y0;
        double  y = y0 + mu*(y1 - y0);

        mpreal  x =  MyQuantile::quantile(y);
        mpreal yp =  MyQuantile::cdf(x);
        mpreal err =  fabs(yp - y);

        printf(
            "cdf(quantile(%.40f)) = %s    err= %s\n",
            y,
            yp.toString().c_str(),
            err.toString().c_str()
        );

        double q = quantile(dist, y);
        printf(
            "       myq = %s\n"
            "     boost = %.40f\n"
            "  cdf(myq) = %s\n"
            "cdf(boost) = %s\n"
            "\n",
            x.toString().c_str(),
            q,
            MyQuantile::cdf(x).toString().c_str(),
            MyQuantile::cdf(q).toString().c_str()
        );
    }
    printf("\n");
}

int main(
    int  ,
    char *[]
)
{
    printf("\n");
    t0();
    t1();
    return 0;
}

