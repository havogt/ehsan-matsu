#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "nr.h"
using namespace std;

// numerical function
class numfunction
{
    double Lambda;
    double mu;
    std::size_t num;
    double delta;
    NRVec<double> data;
    NRVec<double> data2;  // used for the second derivates in spline
    NRVec<double> coords; // need to save coords for NR::spline
    bool isInterpolated;  // keeps track if the interpolation data is up to date;
    double extrapolateZ;
    double extrapolateK;

public:
    numfunction(std::size_t, double = 0, double = 100);
    double operator()(double p);
    std::size_t size() const { return num; }
    double coord(size_t i) const; // that is p
    double computeCoord(size_t i) const;
    std::size_t index(double p) const; // that is i
    double& operator[](const std::size_t i)
    {
        isInterpolated = false;
        if(i >= 0 && i < num)
            return data[i];
        else
        {
            cout << "Error: reading out of bounds!";
            exit(1);
        }
    }
    double operator[](std::size_t i) const
    {
        if(i >= 0 && i < num)
            return data[i];
        else
        {
            cout << "Error: reading out of bounds!";
            exit(1);
        }
    }
    double Fderivative(size_t i) const;
    double Sderivative(size_t i) const;
    void updateInterpolation();
    double optimizedSplint(double p);
    double splint(double p);
    void save(string filename)
    {
        ofstream fout;
        fout.open(filename.c_str());
        for(size_t i = 0; i < num; ++i)
        {
            fout << coord(i) << "\t" << data[i] << endl;
            //                                   fout<< data[i] << endl;
        }
        fout.close();
    }
};
