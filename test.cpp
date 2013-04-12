
#include <rng.h>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdint.h>
#include <mpreal.h>
using namespace mpfr;

// Get us a 1024-bit quantile implementation
namespace Quantile1024 {
    static const int kBits(1024);
    #include <quantile.h>
}

// Marsaglia polar method for the normal distribution
static double normalSample(
    RNG *rng
)
{
    static double stashedValue;
    static bool haveStashedValue = false;
    if(haveStashedValue) {
        haveStashedValue = false;
        return stashedValue;
    }

    while(1) {
        double x0  = 2.0*rng->sample() - 1.0; 
        double y0  = 2.0*rng->sample() - 1.0;
        double r2 = x0*x0 + y0*y0;
        if(r2<=1.0) {
            double scaleFactor = ::sqrt(-2.0*log(r2)/r2);
            stashedValue = y0*scaleFactor;
            haveStashedValue = true;
            return x0*scaleFactor;
        }
    }
}

// Positive half of the normal distribution
static double halfNormalSample(
    RNG *rng
)
{
    double x = normalSample(rng);
    return fabs(x);
}

// Run a simulation of a half-normal distribution, count how many samples "fit"
static double sampleProportion(
    double variance,
    double threshold,
    double nbSamples = 10 * 1000 * 1000
)
{
    RNG rng;
    int count = 0;
    for(int i=0; i<nbSamples; ++i) {
        double x = variance * halfNormalSample(&rng);
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
    return threshold / Quantile1024::quantile(y).toDouble();
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
    for(int i=0; i<=nbSamples; ++i) {
        double y0 = 1e-12;
        double y1 = 1.0 - y0;
        double mu = i/(double)nbSamples;
        double  y = y0 + mu*(y1 - y0);
        mpreal  x =  Quantile1024::quantile(y);
        mpreal yp =  Quantile1024::cdf(x);
        mpreal err =  fabs(yp - y);
        printf(
            "cdf(quantile(%.40f)) = %s    err= %s\n",
            y,
            yp.toString().c_str(),
            err.toString().c_str()
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

