#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#define PI 3.14159265359
using namespace std;
#include "numfunction.h"

void gauleg_d(double x1, double x2, double x[], double w[], int n)
{
#define EPS 3.0e-15
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;

    m = (n + 1) / 2;
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);
    for(i = 0; i < m; i++)
    {
        z = cos(PI * (i + 3. / 4.) / (n + 0.5));
        do
        {
            p1 = 1.0;
            p2 = 0.0;
            for(j = 0; j < n; j++)
            {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while(fabs(z - z1) > EPS);
        x[i] = xm - xl * z;
        x[n - 1 - i] = xm + xl * z;
        w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
}

double MassFunction(double x,
                    double y,
                    double alfa,
                    std::vector<numfunction>& M,
                    int n,
                    int m,
                    double phi,
                    double beta,
                    double delta) // here you define the mass function
{
    double theta = sqrt(x * x + y * y + 2 * x * y * cos(alfa)); // shifted momentum.
    double Mtheta = M[n](theta);                                // for the not extrapolated case.
    double NegativeMtheta = Mtheta;                             // according to the talks.
    double Mx = M[m](x);           // the internal mass function for the not extrapolated case.
    double M2 = (Mtheta * Mtheta); // the shifted mass to the power of two.
    double OmegaN = (2 * n + 1) * PI / beta;                   // the internal matsubara frequency.
    double OmegaM = (2 * m + 1) * PI / beta;                   // the external matsubara frequency.
    double NegativeOmegaN = -1 * OmegaN;                       // the internal matsubara frequency.
    double C2Positive = pow(y, 2) + pow((OmegaN - OmegaM), 2); // denominator in the Coulomb kernel
    double C2Negative
        = pow(y, 2) + pow((NegativeOmegaN - OmegaM),
                          2); // denominator in the Coulomb kernel for the negative frequencies.
    double CK = y / (pow(C2Positive, 2));         // Coulomb kernel.
    double NegativeCK = y / (pow(C2Negative, 2)); // Coulomb kernel for negative frequencies.
    double Coef = 1. / PI;
    double Scalar = (x * x + x * y * cos(alfa)) + (OmegaM * OmegaN); // the Scalar product.
    double NegativeScalar
        = (x * x + x * y * cos(alfa))
          + (OmegaM * NegativeOmegaN); // the Scalar product for the negative frequencies.
    double ScalarD = (x * x) + (OmegaM * OmegaM);           // the denominator of scalar part.
    double Den = sqrt(pow(theta, 2) + pow(OmegaN, 2) + M2); // denominator of the equation (Energy).
    double NegativeDen = sqrt(pow(theta, 2) + pow(NegativeOmegaN, 2)
                              + M2); // denominator of the equation for negative frequencies.
    double SP = Scalar / ScalarD;    // scalar part complete.
    double NegativeSP = NegativeScalar / ScalarD; // scalar part with negative numarator.
    double Main = (Mtheta - (Mx * SP)) / Den;     // main part of the integral
    double NegativeMain
        = (Mtheta - (Mx * NegativeSP)); // main part for negative matsubara frequencies.

    double Solution = Coef * ((CK * Main) + (NegativeCK * NegativeMain));

    return Solution;
}


double integral(double y, double beta, int n, std::vector<numfunction>& M, double phi)
{
    double OmegaN = (2 * n + 1) * PI / beta;
    double Mass = M[n](y);
    double M2 = Mass * Mass;
    double Den = sqrt(pow(y, 2.) + pow(OmegaN, 2.) + M2);
    double Coef = -3 / (PI * beta);
    double Conds = 2. * y * Coef * (Mass / Den);

    return Conds;
}
double condestatecalculation(std::vector<numfunction>& M,
                             int Mats,
                             double* wx,
                             double* wy,
                             double* x,
                             double* y,
                             size_t ng,
                             size_t mg,
                             double beta,
                             double phi)
{
    double ro = 0;
    for(size_t k = 0; k < ng; ++k)
    {
        for(int n = 0; n < Mats; ++n)
        {
            ro += wx[k] * integral(x[k], beta, n, M, phi);
        }
    }

    return ro;
}


void save(int Mats, std::vector<numfunction>& f, char* filename, int beta)
{
    ofstream ausgabe;
    ausgabe.open(filename);
    for(int m = 0; m < Mats; m++)
    {

        for(size_t i = 0; i < f[m].size(); i++)
        {
            double phi = 0;
            int mtilde = m - Mats;
            double Mtilde = ((PI * ((2 * mtilde) + 1)) - phi) / beta;


            ausgabe << f[m].coord(i) << "\t" << m << "\t" << Mtilde << "\t" << f[m][i] << endl;
        }
    }
    ausgabe.close();
}


void test(char* filename, size_t maxiter, int Mats, double beta)
{
    size_t ng = 30.;
    size_t mg = 25.;
    double epsilon = 1e-5;
    double wx[ng];
    double x[ng];
    double wy[mg];
    double y[mg];
    int gridPoints = 30.;
    double irCutoff = 1e-2;
    double uvCutoff = 20.;
    double q0 = log10(irCutoff), q1 = log10(uvCutoff);
    double z0 = 1e-2, z1 = 2 * PI;
    std::vector<numfunction> fnew(Mats, numfunction(gridPoints, irCutoff, uvCutoff));
    std::vector<numfunction> fold(Mats, numfunction(gridPoints, irCutoff, uvCutoff));
    gauleg_d(q0, q1, x, wx, ng);
    gauleg_d(z0, z1, y, wy, mg);
    double update = 0.0001;
    double phi = 0.;
    double delta = 0.;


    for(int m = 0; m < Mats; ++m)
    {
        for(size_t i = 0; i < fold[m].size(); ++i)
        {
            fold[m][i] = 1. / (1. + pow(fold[m].coord(i), 4));
        }
    }

    // update weights for log integration
    for(size_t k = 0; k < ng; ++k)
    {
        x[k] = pow(10, x[k]);
        wx[k] *= x[k] * log(10);
    }

    for(size_t l = 0; l < maxiter; l++)
    {
        //   for(int m=0; m<Mats; m++)
        //
        //   {fold[m].updateInterpolation();
        //   }
        // #pragma omp parallel for

        for(int m = 0; m < Mats; ++m)
        {
            for(size_t i = 0; i < fold[m].size(); i++)
            {
                double h = fnew[0].coord(i);
                double r = 0.;

                for(size_t k = 0; k < ng; ++k)
                {
                    double rf = 0.;

                    for(size_t j = 0; j < mg; ++j)
                    {

                        for(int n = 0; n < Mats; ++n)
                        {
                            double F = MassFunction(h, x[k], y[j], fold, n, m, phi, beta, delta);
                            rf += wy[j] * F;
                        }
                    }
                    r += rf * wx[k];
                }

                fnew[m][i] = r;
            }
        }


        double cp = 0;
        double cd = 0;
        // do convergence test
        for(int m = 0; m < Mats; ++m)
        {
            for(size_t i = 0; i < fold[m].size(); i++)
            {
                double cprime = fabs(fold[m][i] - fnew[m][i]);

                if(cprime > cd)
                {
                    cd = cprime;
                }
            }
            for(size_t i = 0; i < fold[m].size(); i++)
            {
                double c = fabs((fold[m][i] - fnew[m][i]) / fold[m][i]);
                if(c > cp)
                {
                    cp = c;
                }
            }
        }
        cout << "cp"
             << "..." << cp << "..."
             << "number"
             << "...." << cd << "...."
             << "Abs" << l << endl;

        if(cd < 1e-6)
        {
            cout << "iteration"
                 << "..." << l << endl;
            break;
        }

        if(cp < epsilon)
        {
            cout << "iteration"
                 << "..." << l << endl;
            break;
        }


        if(l % 100 == 0)
        {
            save(Mats, fnew, filename, beta);
        }

        for(int m = 0; m < Mats; m++)
        {
            for(size_t i = 0; i < fold[m].size(); i++)
            {
                fold[m][i] = (1 - update) * fold[m][i] + update * fnew[m][i];
            }
        }
    }

    save(Mats, fnew, filename, beta);
    cout << "condenstate"
         << "..." << condestatecalculation(fnew, Mats, wx, wy, x, y, ng, mg, beta, phi) << endl;
}
int main(int argc, char** argv)
{
    if(argc < 5)
    {
        cout << "Error!\nUsage: " << argv[0] << " <filename> <maxiter> <Mats> <beta>" << endl;
        exit(1);
    }
    test(argv[1], atoi(argv[2]), atoi(argv[3]), atof(argv[4]));
    return 0;
}
