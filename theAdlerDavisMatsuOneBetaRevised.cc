#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define PI 3.14159265359
using namespace std;
#include "numfunction.h"
#include "matsextrapolation.h"
#include "parametrizedMassfunction.h"
using namespace parametrizedMassfunction;
ofstream ausgabe_condensate;

const char path[255] = "MatsuBaraRes/";
const char type[255] = ".txt";

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


string makeFilname(char* name, const char* suffix)
{
    char buffer1[256]; // <- danger, only storage for 256 characters.
    strncpy(buffer1, path, sizeof(buffer1));
    strncat(buffer1, name, sizeof(buffer1));
    strncat(buffer1, suffix, sizeof(buffer1));
    strncat(buffer1, type, sizeof(buffer1));
    return string(buffer1);
}


double MassFunction(double x,
                    double y,
                    double alfa,
                    std::vector<numfunction>& M,
                    int n,
                    int m,
                    double phi,
                    double beta,
                    double delta)
{
    double theta = sqrt(x * x + y * y - 2 * x * y * cos(alfa));
    double Mtheta = extrapolateMats(n,
                                    [&](int n_mats)
                                    {
                                        return M[n_mats](theta);
                                    },
                                    M.size());
    double Mass = extrapolateMats(n,
                                  [&](int n_mats)
                                  {
                                      return M[n_mats](x);
                                  },
                                  M.size());
    double M2 = (Mtheta * Mtheta);
    int mtilde = m;
    int ntilde = n;
    int ntildeNeg = -n - 1;
    double Ntilde = ((PI * ((2 * ntilde) + 1)) - phi) / beta;
    double Mtilde = ((PI * ((2 * mtilde) + 1)) - phi) / beta;
    double NtildeNeg = ((PI * ((2 * ntildeNeg) + 1)) - phi) / beta;
    double xi = ((x * x) + (y * y) - (2 * fabs(x * y) * cos(alfa)))
                + (pow(Ntilde, 2)); // this is theta to the power of two ;
    double xiNeg = ((x * x) + (y * y) - (2 * fabs(x * y) * cos(alfa)))
                   + (pow(NtildeNeg, 2)); // this is theta negative to the power of two ;
    double omega
        = (((x * x) - fabs(x * y) * cos(alfa)) + (Mtilde * Ntilde))
          / ((x * x) + (pow(Mtilde, 2))); // this is the scalar product part in the nominator
    double omegaNeg
        = (((x * x) - fabs(x * y) * cos(alfa)) + (Mtilde * NtildeNeg))
          / ((x * x)
             + (pow(Mtilde, 2))); // this is the scalar product negative part in the nominator
    double zetta
        = pow((2 * PI * ((ntilde - mtilde) / beta)),
              2); // this is the matsubara part of the kernel we just shifted the momentum part
    double zettaNeg = pow(
        (2 * PI * ((ntildeNeg - mtilde) / beta)),
        2); // this is the matsubara part of the kernel negative we just shifted the momentum part
    double xiprime = (pow(((y * y) + zetta), 2) + delta); // this is the coloumb kernel
    double xiprimeNeg = (pow(((y * y) + zettaNeg), 2)
                         + delta); // this is the coloumb kernel with the negative part
    double B = (y / (xiprime)) * (Mtheta / (beta * sqrt(xi + M2)));
    double C = (y / (xiprime)) * (omega / (beta * sqrt(xi + M2)));
    double BNeg = (y / (xiprimeNeg)) * (Mtheta / (beta * sqrt(xiNeg + M2)));
    double CNeg = (y / (xiprimeNeg)) * (omegaNeg / (beta * sqrt(xiNeg + M2)));
    double First = C + CNeg;
    double Second = B + BNeg;
    double result = (1. / PI) * (Second - (Mass * First));

    if(::isnan(C))
    {
        cout << "For C"
             << "theta"
             << "....." << theta << "..." << Mtheta << "..."
             << "Mtheta"
             << "..." << xiprime << "coloumb kernel" << endl;
        exit(1);
    }

    return result;
}

double integral(double y, double beta, int n, std::vector<numfunction>& M, double phi)
{
    int N = 3.;
    int ntilde = n;
    int ntildeNeg = -n - 1;
    double Ntilde = ((PI * ((2 * ntilde) + 1)) - phi) / beta;
    double NtildeNeg = ((PI * ((2 * ntildeNeg) + 1)) - phi) / beta;
    double Mass = extrapolateMats(n,
                                  [&](int n_mats)
                                  {
                                      return M[n_mats](y);
                                  },
                                  M.size());
    double M2 = Mass * Mass;
    double xi = (y * y) + (pow(Ntilde, 2));       // this is momentum to the power of two ;
    double xiNeg = (y * y) + (pow(NtildeNeg, 2)); // this is momentum to the power of two ;
    double A = -(N / (PI * beta)) * (y * Mass) * ((1 / sqrt(M2 + xi)) + (1 / sqrt(M2 + xiNeg)));

    return A;
}


double condensateCalculation(std::vector<numfunction>& M,
                             int number,
                             double* wx,
                             double* wy,
                             double* x,
                             double* y,
                             size_t ng,
                             size_t mg,
                             double beta,
                             double phi) // this is for the calculation of condestate
{
    double ro = 0;
    for(size_t k = 0; k < ng; ++k)
    {
        for(int n = 0; n < number; ++n)
        {
            ro += wx[k] * integral(x[k], beta, n, M, phi);
        }
    }

    return ro;
}


void save(int Mats,
          std::vector<numfunction>& f,
          char* filename1,
          int beta,
          double phi,
          int steps,
          double A,
          double B,
          int number,
          double delta)

{


    ofstream ausgabe;
    ausgabe.open(makeFilname(filename1, "_rawData"));
    ausgabe.precision(10);
    ausgabe.setf(std::ios::scientific);

    ausgabe << "#" << beta << "\t" << Mats << "\t" << A << "\t" << B << "\t" << number << "\t"
            << delta << endl;

    for(int m = 0; m < Mats; m++)
    {


        for(size_t i = 0; i < f[m].size(); i++)
        {
            int mtilde = m;
            double Mtilde = ((PI * ((2 * mtilde) + 1)) - phi) / beta;


            ausgabe << f[m].coord(i) << "\t" << m << "\t" << Mtilde << "\t" << f[m][i] << "\t"
                    << sqrt(pow(f[m].coord(i), 2) + pow(Mtilde, 2)) << "\t" << beta << endl;
        }
    }
    ausgabe.close();

    double Mtilde = ((PI * ((2 * number) + 1)) - phi) / beta;
    ausgabe.open(makeFilname(filename1, "_extrapolatedData"));

    for(int m = 0; m < number; m++)
    {
        for(int j = 0; j < steps; j++)
        {
            double x = pow(10, log10(A) + log10(B / A) / (steps - 1) * j);

            auto reducedMassFunction = [&](int n_mats)
            {
                return f[n_mats](x);
            };

            double F = extrapolateMats(m, reducedMassFunction, Mats);

            ausgabe << x << "\t" << m << "\t" << Mtilde << "\t" << F << "\t" << pn(m, beta) << "\n"
                    << endl;
        }
        ausgabe << "\n";
    }
    ausgabe.close();
}


void Matsu(char* filename1, char* filename3, size_t maxiter, int Mats)
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
    double uvCutoff = 20.;
    double q0 = log10(irCutoff), q1 = log10(uvCutoff);
    double z0 = 1e-2, z1 = 2 * PI;
    std::vector<numfunction> fnew(Mats, numfunction(gridPoints, irCutoff, uvCutoff));
    std::vector<numfunction> fold(Mats, numfunction(gridPoints, irCutoff, uvCutoff));
    int number = 2 * Mats;
    gauleg_d(q0, q1, x, wx, ng);
    gauleg_d(z0, z1, y, wy, mg);
    double update = 0.0005;
    double phi = 0.;
    double delta = 0.;
    double beta = 15;

    cout << "Starting solving for beta = " << beta << "\t"
         << "T"
         << "\t" << 803 / beta << "MeV"
         << "\n";

    for(int m = 0; m < Mats; ++m)
    {
        for(size_t i = 0; i < fold[m].size(); ++i)
        {
            fold[m][i] = 0.1;
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

#pragma omp parallel for collapse(2)
        for(int m = 0; m < Mats; ++m)
        {
            for(size_t i = 0; i < gridPoints; i++)
            {
                double h = fnew[0].coord(i);
                double r = 0.;

                for(size_t k = 0; k < ng; ++k)
                {
                    double rf = 0;

                    for(size_t j = 0; j < mg; ++j)
                    {
                        for(int n = 0; n < number; ++n)
                        {
                            // auto reducedMassFunction =[&](int n_mats){return MassFunction
                            // (h,x[k],y[j],fold,n_mats,m,phi,beta,delta) ;};

                            // double F = extrapolateMats(n,reducedMassFunction,Mats);

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
                double c = fabs((fold[m][i] - fnew[m][i]) / (fold[m][i]));
                if(c > cp)
                {
                    cp = c;
                }
            }
        }
        cout << "cp"
             << "..." << cp << "..."
             << "ABS"
             << "..." << cd << "...."
             << "Iter"
             << "...." << l << endl;

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


        if(l % 10 == 0)
        {

            save(Mats, fnew, filename1, beta, phi, gridPoints, irCutoff, uvCutoff, number, delta);
        }

        for(int m = 0; m < Mats; m++)
        {
            for(size_t i = 0; i < fold[m].size(); i++)
            {
                fold[m][i] = (1 - update) * fold[m][i] + update * fnew[m][i];
            }

            fold[m].updateInterpolation();
        }
    }

    save(Mats, fnew, filename1, beta, phi, gridPoints, irCutoff, uvCutoff, number, delta);
    cout << "condensate"
         << "..." << condensateCalculation(fnew, number, wx, wy, x, y, ng, mg, beta, phi) << "\t"
         << beta << "\t" << (803 / beta) << "\t"
         << "MeV"
         << "\t" << maxiter << "\t" << Mats << "\t" << number << endl;
    ausgabe_condensate.open(filename3, std::fstream::app);
    ausgabe_condensate << beta << "\t" << fold[0](irCutoff) << "\t"
                       << condensateCalculation(fnew, number, wx, wy, x, y, ng, mg, beta, phi)
                       << "\t" << endl;
    ausgabe_condensate.flush();
    ausgabe_condensate.close();
}

int main(int argc, char** argv)
{
    if(argc < 5)
    {
        cout << "Error!\nUsage: " << argv[0] << " <filename1><filename3><maxiter><Mats>" << endl;
        exit(1);
    }

    Matsu(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));

    return 0;
}
