/*

    To compute a quantile with a given precision of 2048 bits, do this:

        namespace Quantile128 {
            static const int kBits(128);
            #include <quantile.h>
        }

        ...

        double myQuantile = Quantile128::quantile(0.1234);

    Technique to compute the quantile was lifted from P.J Acklam's work:

            http://home.online.no/~pjacklam/notes/invnorm

    Basic idea : quick lo-fidelity first pass of the quantile using
    multi-branch rational approximators, followed by a refinement
    step using Halley's method.

*/

static const mpreal kOne (1, kBits);
static const mpreal kTwo (2, kBits);
static const mpreal kFour(4, kBits);
static const mpreal kZero(0, kBits);
static const mpreal kMinusOne(-1, kBits);
static const mpreal kMinusTwo(-2, kBits);
static const mpreal kPi = const_pi(kBits);
static const mpreal kNAN = sqrt(kMinusOne);
static const mpreal kPointFive = kOne / kTwo;
static const mpreal kSqrtOfTwoPi = sqrt(kTwo * kPi);
static const mpreal kPlusInf = const_infinity( 1, kBits);
static const mpreal kMinusInf = const_infinity(-1, kBits);
static const mpreal kOneOverSqrtOfTwo = kOne / sqrt(kTwo);
static const mpreal kPointTwentyFive = kPointFive * kPointFive;
static const mpreal kMinusOneOverSqrtOfTwo = kMinusOne / sqrt(kTwo);

static inline mpreal cdf(
    const mpreal &x
)
{
    return kPointFive*(kOne + erf(x*kOneOverSqrtOfTwo));
}

static inline mpreal refine(
    const mpreal &x0,
    const mpreal &y0
)
{
    mpreal y = y0;
    for(int i=0; i<10; ++i) {
        mpreal e = kPointFive * erfc(y*kMinusOneOverSqrtOfTwo) - x0;
        mpreal u = e * kSqrtOfTwoPi * exp(y*y*kPointFive);
        mpreal nextY = y - u/(kOne + y*u*kPointFive);
        if(nextY==y) break;
        y = nextY;
    }

    return y;
}

static inline mpreal quantileL(
    const mpreal &x
)
{
    static const mpreal c1(-7.784894002430293e-03, kBits);
    static const mpreal c2(-3.223964580411365e-01, kBits);
    static const mpreal c3(-2.400758277161838e+00, kBits);
    static const mpreal c4(-2.549732539343734e+00, kBits);
    static const mpreal c5( 4.374664141464968e+00, kBits);
    static const mpreal c6( 2.938163982698783e+00, kBits);

    static const mpreal d1( 7.784695709041462e-03, kBits);
    static const mpreal d2( 3.224671290700398e-01, kBits);
    static const mpreal d3( 2.445134137142996e+00, kBits);
    static const mpreal d4( 3.754408661907416e+00, kBits);

    mpreal q = sqrt(kMinusTwo*log(x));
    mpreal num = ((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6;
    mpreal den = (((d1*q + d2)*q + d3)*q + d4)*q + kOne;
    return (num/den);
}

static inline mpreal quantileM(
    const mpreal &x
)
{
    static const mpreal a1(-3.969683028665376e+01, kBits);
    static const mpreal a2( 2.209460984245205e+02, kBits);
    static const mpreal a3(-2.759285104469687e+02, kBits);
    static const mpreal a4( 1.383577518672690e+02, kBits);
    static const mpreal a5(-3.066479806614716e+01, kBits);
    static const mpreal a6( 2.506628277459239e+00, kBits);

    static const mpreal b1(-5.447609879822406e+01, kBits);
    static const mpreal b2( 1.615858368580409e+02, kBits);
    static const mpreal b3(-1.556989798598866e+02, kBits);
    static const mpreal b4( 6.680131188771972e+01, kBits);
    static const mpreal b5(-1.328068155288572e+01, kBits);

    mpreal q = x - kPointFive;
    mpreal r = q * q;

    mpreal num = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q;
    mpreal den = (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + kOne);
    return (num/den);
}

static inline mpreal quantileR(
    const mpreal &x
)
{
    return -quantileL(kOne - x);
}

static inline mpreal quantile(
    const mpreal &x
)
{
    // Get the weirdos out of the way
    if(sgn(x)<0)      return kNAN;
    if(iszero(x))     return kMinusInf;
    if(kPointFive==x) return kZero;
    if(kOne==x)       return kPlusInf;
    if(kOne<x)        return kNAN;

    // Pick a branch
    static const mpreal kLo(0.02425, kBits);
    static const mpreal kHi = kOne - kLo;
    mpreal r = (
        (x < kLo) ? quantileL(x) :
        (x < kHi) ? quantileM(x) :
                    quantileR(x)
    );
    return refine(x, r);
}

