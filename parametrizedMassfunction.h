#define PI 3.14159265359

namespace parametrizedMassfunction
{
static double Massfunc(double p)
{
    return 1. / (1. + pow(p, 4));
}

static double momentum(double p_t, double p_n)
{
    return sqrt(p_t * p_t + p_n * p_n);
}

static double pn(int n, double beta)
{
    return PI * (2. * n + 1) / beta;
}
}
