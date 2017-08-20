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


struct MassfunctionStructure
{
    double MassFunction;
};

MassfunctionStructure MassFunction(double x,
                                   double y,
                                   double alfa,
                                   std::vector<numfunction>& M,
                                   int n,
                                   int m,
                                   double beta,
                                   int Mats,
                                   double delta)
{
    MassfunctionStructure result;

    double theta = sqrt(x * x + y * y + 2 * x * y * cos(alfa)); // momentum part (y=q,x=k)
    double Mtheta = M[n](theta); // M defined by the momentum part and number part
    double Mthetaint = M[m](x);  // M defined by the momentum part and number part
    double M2 = (Mtheta * Mtheta);
    int ntildeNeg = -n - 1;
    double Ntilde = (PI * ((2 * n) + 1)) / beta;
    double Mtilde = (PI * ((2 * m) + 1)) / beta;
    double NtildeNeg = PI * ((2 * ntildeNeg) + 1) / beta;
    double xi = (x * x) + (y * y) + (2 * x * y * cos(alfa))
                + pow(Ntilde, 2); // this is theta to the power of two ;
    double xiNeg = (x * x) + (y * y) + (2 * x * y * cos(alfa))
                   + pow(NtildeNeg, 2); // this is theta negative to the power of two ;
    double omega = (((x * x) + (x * y) * cos(alfa)) + (Mtilde * Ntilde))
                   / ((x * x) + (pow(Mtilde, 2))); // this is the scalar product part
    double omegaNeg = (((x * x) + (x * y) * cos(alfa)) + (Mtilde * NtildeNeg))
                      / ((x * x) + pow(Mtilde, 2)); // this is the scalar product negative part
    double zetta
        = pow((2 * PI * ((n - m) / beta)),
              2); // this is the matsubara part of the kernel we just shifted the momentum part
    double zettaNeg = pow(
        (2 * PI * ((ntildeNeg - m) / beta)),
        2); // this is the matsubara part of the kernel negative we just shifted the momentum part
    double xiprime = (pow(((y * y) + zetta), 2) + delta); // this is the coloumb kernel
    double xiprimeNeg = (pow(((y * y) + zettaNeg), 2)
                         + delta); // this is the coloumb kernel with the negative part
    double B = (y / (xiprime)) * (1. / (beta * sqrt(xi + M2))) * (Mtheta - (Mthetaint * omega));
    double C
        = (y / (xiprimeNeg)) * (1. / (beta * sqrt(xiNeg + M2))) * (Mtheta - (Mthetaint * omegaNeg));
    double mass = B + C;

    result.MassFunction = mass;

    return result;
};

double integral(double y, double beta, int n, std::vector<numfunction>& M)
{
    int N = 3.;
    double Ntilde = (PI * ((2 * n) + 1)) / beta;
    double NtildeNeg = (PI * (-(2 * n) - 1)) / beta;
    double Mass = M[n](y);
    double M2 = Mass * Mass;
    double xi = (y * y) + (pow(Ntilde, 2));       // this is momentum to the power of two ;
    double xiNeg = (y * y) + (pow(NtildeNeg, 2)); // this is momentum to the power of two ;
    double A = -(N / (PI * beta)) * (y * Mass) * ((1. / sqrt(M2 + xi)) + (1. / sqrt(M2 + xiNeg)));

    return A;
}

double condestatecalculation(std::vector<numfunction>& M,
                             int Mats,
                             double* wx,
                             double* wy,
                             double* x,
                             double* y,
                             size_t ng,
                             size_t mg,
                             double beta)
{
    double ro = 0;
    for(size_t k = 0; k < ng; ++k)
    {
        for(int n = 0; n < Mats; ++n)
        {
            ro += wx[k] * integral(x[k], beta, n, M);
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
    size_t ng = 30;
    size_t mg = 25;
    double epsilon = 1e-5;
    double wx[ng];
    double x[ng];
    double wy[mg];
    double y[mg];
    int gridPoints = 30;
    double irCutoff = 1e-2;
    double uvCutoff = 1e2;
    double q0 = log10(irCutoff), q1 = log10(uvCutoff);
    double z0 = 1e-2, z1 = 2 * PI;
    std::vector<numfunction> fnew(2 * Mats, numfunction(gridPoints, irCutoff, uvCutoff));
    std::vector<numfunction> fold(2 * Mats, numfunction(gridPoints, irCutoff, uvCutoff));
    gauleg_d(q0, q1, x, wx, ng);
    gauleg_d(z0, z1, y, wy, mg);
    double update = 0.0003;
    double phi = 0;
    double delta = 0;
    double gtilde = 0;

    for(int m = 0; m < Mats; ++m)
    {
        for(size_t i = 0; i < fold[m].size(); ++i)
        {
            fold[m][i] = 1. / (1. + pow(fold[m].coord(i), 4));

            //  cout << fold[m].coord(i) << "\t ";
        }
    }
    cout << endl;

    // update weights for log integration
    for(size_t k = 0; k < ng; ++k)
    {
        x[k] = pow(10, x[k]);
        wx[k] *= x[k] * log(10);
    }

    for(size_t l = 0; l < maxiter; l++)
    {
#pragma omp parallel for
        for(int m = 0; m < Mats; ++m)
        {
            fold[m].updateInterpolation();

            for(size_t i = 0; i < fold[m].size(); i++)
            {
                double h = fnew[m].coord(i);
                double r = 0;

                for(size_t k = 0; k < ng; ++k)
                {
                    double rk = 0;

                    for(size_t j = 0; j < mg; ++j)
                    {
                        for(int n = 0; n < Mats; ++n)
                        {
                            MassfunctionStructure result
                                = MassFunction(h, x[k], y[j], fold, n, m, beta, Mats, delta);

                            rk += wy[j] * result.MassFunction;
                        }
                    }
                    r += rk * wx[k];
                }

                fnew[m][i] = r / (PI);
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
         << "..." << condestatecalculation(fnew, Mats, wx, wy, x, y, ng, mg, beta) << endl;
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
