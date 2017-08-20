
#include <iostream>
#include <iomanip>
#include <vector>
#include <exception>
#include <stdexcept>
#include "numfunction.h"


// constructor of numerical function
numfunction::numfunction(std::size_t N, double m, double Lam) : data(N), data2(N), coords(N)
{
    cout << "#using logarithmic numfunction" << endl;
    if(N < 2)
        throw std::invalid_argument("insufficient number of components!");
    mu = m;
    Lambda = Lam;
    num = N;
    delta = Lambda / mu;
    isInterpolated = false;

    for(size_t i = 0; i < num; ++i)
    {
        coords[i] = computeCoord(i);
    }
}
std::size_t numfunction::index(double p) const
{
    //	double x = trunc(((log(p / mu)) / (log(delta))) * (num - 1));
    //	return x;

    // for many cases this is much faster than the explicit calculation of the index above
    int k;
    int n = coords.size();
    int klo = 0;
    int khi = n - 1;
    while(khi - klo > 1)
    {
        k = (khi + klo) >> 1;
        if(coords[k] > p)
            khi = k;
        else
            klo = k;
    }
    return klo;
}

double numfunction::coord(size_t i) const
{
    return coords[i];
}

double numfunction::computeCoord(size_t i) const
{
    return mu * pow(delta, (double) i / (double) (num - 1));
}

double inline numfunction::optimizedSplint(double p)
{
    size_t klo = index(p);
    size_t khi = klo + 1;
    double h = coords[khi] - coords[klo];
    double a = (coords[khi] - p) / h;
    double b = (p - coords[klo]) / h;
    return a * data[klo] + b * data[khi]
           + ((a * a * a - a) * data2[klo] + (b * b * b - b) * data2[khi]) * (h * h) / 6.0;
}

double numfunction::Fderivative(size_t i) const
{
    if(i == num - 1)
    {
        i = num - 2;
    }

    return (data[i + 1] - data[i]) / (coord(i + 1) - coord(i));
}

double numfunction::Sderivative(size_t i) const
{
    if(i == 0)
    {
        i = 1;
    }
    if(i == num - 1)
    {
        i = num - 2;
    }

    return (data[i + 1] - 2 * data[i] + data[i - 1]) / (pow((coord(i + 1) - coord(i - 1)), 2));
}

double numfunction::splint(double p)
{
    int k;
    DP h, b, a;

    int n = coords.size();
    int klo = 0;
    int khi = n - 1;
    while(khi - klo > 1)
    {
        k = (khi + klo) >> 1;
        if(coords[k] > p)
            khi = k;
        else
            klo = k;
    }
    h = coords[khi] - coords[klo];
    a = (coords[khi] - p) / h;
    b = (p - coords[klo]) / h;
    return a * data[klo] + b * data[khi]
           + ((a * a * a - a) * data2[klo] + (b * b * b - b) * data2[khi]) * (h * h) / 6.0;
}

void numfunction::updateInterpolation()
{
    NR::spline(coords, data, 1e30, 1e30, data2);

    extrapolateK = (log(data[num - 2] / data[num - 1])) / (log(coord(num - 1) / coord(num - 2)));
    extrapolateZ = (data[num - 2] * pow(coord(num - 2), extrapolateK));

    isInterpolated = true;
}

// function call for numerical function
double numfunction::operator()(double p)
{
    p = fabs(p);

    if(!isInterpolated)
    {
        updateInterpolation();
    }

    if(p > Lambda)
    {
        return (extrapolateZ / (pow(p, extrapolateK)));
    }
    else if(p < mu)
    {
        return data[0];
    }
    else
    {
        //		return optimizedSplint( p );
        //		return splint(p);
        double result;
        NR::splint(coords, data, data2, p, result);
        return result;
    }
}
