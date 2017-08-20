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

double m_tilde(int m, double phi, double beta)
{
    return (PI * ((2. * (double) m) + 1.) - phi) / beta;
}

double xiprime_f(int n, int l, double beta, double q, double delta)
{
    double omegaN_minus_omegaL_squared
        = pow((2. * PI * ((n - l) / beta)),
              2); // this is the matsubara part of the kernel we just shifted the momentum part
    double xiprime
        = pow(((q * q) + omegaN_minus_omegaL_squared), 2) + delta; // this is the coloumb kernel
    return xiprime;
}

double xi_f(double p, double q, double alfa, double OmegaN)
{
    double xi = (p * p) + (q * q) - (2 * fabs(p * q) * cos(alfa))
                + OmegaN * OmegaN; // |q-p|^2 + OmegaN^2 ;
    return xi;
}

double omega_f(double p, double q, double alfa, double OmegaL, double OmegaN)
{
    double omega
        = (((p * p) - fabs(p * q) * cos(alfa)) + (OmegaL * OmegaN))
          / ((p * p) + (pow(OmegaL, 2))); // this is the scalar product part in the nominator
    return omega;
}

double MassFunction(double p,    // outer momentum
                    double q,    // inner momentum
                    double alfa, // angle between p and q
                    std::vector<numfunction>& M,
                    int n, // inner Matsubara momentum
                    int l, // outer Matsubara momentum
                    double phi,
                    double beta,
                    double delta)
{
    double theta = sqrt(p * p + q * q - 2 * p * q * cos(alfa));
    // Mtheta = M(q-p,n)
    double Mtheta = extrapolateMats(n,
                                    [&](int n_mats)
                                    {
                                        return M[n_mats](theta);
                                    },
                                    M.size());
    double Mp = extrapolateMats(l,
                                [&](int n_mats)
                                {
                                    return M[n_mats](p);
                                },
                                M.size());
    double M2 = (Mtheta * Mtheta);
    int n_neg = -n - 1;

    double OmegaL = m_tilde(l, phi, beta);
    double OmegaN_neg = m_tilde(n_neg, phi, beta);
    double OmegaN = m_tilde(n, phi, beta);

    double xi = xi_f(p, q, alfa, OmegaN);
    double xi_neg = xi_f(p, q, alfa, OmegaN_neg);

    double omega = omega_f(p, q, alfa, OmegaL, OmegaN);
    double omega_neg = omega_f(p, q, alfa, OmegaL, OmegaN_neg);

    double xi_prime = xiprime_f(n, l, beta, q, delta);
    double xi_prime_neg = xiprime_f(n_neg, l, beta, q, delta);

    double u_plus = (q / xi_prime) * (Mtheta / (beta * sqrt(xi + M2)));
    double u_neg = (q / xi_prime_neg) * (Mtheta / (beta * sqrt(xi_neg + M2)));

    double v_plus = (q / xi_prime) * (omega / (beta * sqrt(xi + M2)));
    double v_neg = (q / xi_prime_neg) * (omega_neg / (beta * sqrt(xi_neg + M2)));

    double v_tot = v_plus + v_neg;
    double u_tot = u_plus + u_neg;

    double result = (1. / PI) * (u_tot - (Mp * v_tot));

    if(::isnan(v_plus))
    {
        cout << "For C"
             << "theta"
             << "....." << theta << "..." << Mtheta << "..."
             << "Mtheta"
             << "..." << xiprime_f << "coloumb kernel" << endl;
        exit(1);
    }

    return result;
}


double integral(double y, double beta, int n, std::vector<numfunction>& M, double phi)
{
    double N = 3.;
    int ntildeNeg = -n - 1;
    double Ntilde = m_tilde(n, phi, beta);
    double NtildeNeg = m_tilde(ntildeNeg, phi, beta);
    double Mass = extrapolateMats(n,
                                  [&](int n_mats)
                                  {
                                      return M[n_mats](y);
                                  },
                                  M.size());
    double M2 = Mass * Mass;
    double xi = (y * y) + (pow(Ntilde, 2));       // this is momentum to the power of two ;
    double xiNeg = (y * y) + (pow(NtildeNeg, 2)); // this is momentum to the power of two ;
    double A = -(N / (PI * beta)) * (y * Mass) * ((1. / sqrt(M2 + xi)) + (1. / sqrt(M2 + xiNeg)));

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
                             double phi) // this is for the calculation of condensate
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
            double Mtilde = m_tilde(m, phi, beta);

            ausgabe << f[m].coord(i) << "\t" << m << "\t" << Mtilde << "\t" << f[m][i] << "\t"
                    << sqrt(pow(f[m].coord(i), 2) + pow(Mtilde, 2)) << "\t" << beta << endl;
        }
    }
    ausgabe.close();

    double Mtilde = m_tilde(number, phi, beta);
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

            ausgabe << x << "\t" << m << "\t" << Mtilde << "\t" << F << "\t"
                    << m_tilde(m, phi, beta) << "\n" << endl;
        }
        ausgabe << "\n";
    }
    ausgabe.close();
}


void Matsu(char* filename1, char* filename3, size_t maxiter, int Mats, double beta)
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
    int number = 4 * Mats;
    gauleg_d(q0, q1, x, wx, ng);
    gauleg_d(z0, z1, y, wy, mg);
    double update = 0.0005;
    double phi = 0.;
    double delta = 0.;

    cout << "Starting solving for beta = " << beta << ", T = " << 803. / beta << "MeV\n";

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

    for(size_t iter = 0; iter < maxiter; iter++)
    {
#pragma omp parallel for collapse(2)
        for(int l = 0; l < Mats; ++l)
        {
            for(size_t i = 0; i < gridPoints; i++)
            {
                double p = fnew[0].coord(i);
                double r = 0.;

                for(size_t k = 0; k < ng; ++k)
                {
                    double q = x[k];
                    double rf = 0;
                    for(size_t j = 0; j < mg; ++j)
                    {
                        double alpha = y[j];
                        for(int n = 0; n < number; ++n)
                        {
                            double F = MassFunction(p, q, alpha, fold, n, l, phi, beta, delta);

                            rf += wy[j] * F;
                        }
                    }
                    r += rf * wx[k];
                }
                fnew[l][i] = r;
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
             << "...." << iter << endl;

        if(iter == 14)
        {
            if(fabs(cp - 9.31688) > 1e-4)
            {
                std::cout << "error" << std::endl;
                exit(1);
            }
        }

        if(cd < 1e-6)
        {
            cout << "iteration"
                 << "..." << iter << endl;

            break;
        }

        if(cp < epsilon)
        {
            cout << "iteration"
                 << "..." << iter << endl;

            break;
        }

        if(iter % 10 == 0)
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

    ofstream ausgabe_condensate;
    ausgabe_condensate.open(filename3, std::fstream::app);
    ausgabe_condensate << beta << "\t" << fold[0](irCutoff) << "\t"
                       << condensateCalculation(fnew, number, wx, wy, x, y, ng, mg, beta, phi)
                       << "\t" << endl;
    ausgabe_condensate.flush();
    ausgabe_condensate.close();
}

int main(int argc, char** argv)
{
    if(argc < 6)
    {
        cout << "Error!\nUsage: " << argv[0] << " <filename1> <filename3> <maxiter> <Mats> <beta>"
             << endl;
        exit(1);
    }

    Matsu(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atof(argv[5]));

    return 0;
}
