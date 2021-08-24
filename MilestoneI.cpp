// File to solve milestone I of Ast5220

#include <cmath>
#include <iostream>
#include <fstream>

// Global parameters
#define pi 3.14159265 // Pi
#define h 0.67        // 
#define T_CMB0 2.7255 // Present CMB temperature [K]
#define OM_b0 0.05    // Present baryon density parameter
#define OM_CDM0 0.267 // Present CDM density parameter
#define OM_k0 0       // Present k density parameter
#define N_eff 0       // Effective number of massles neutrinoes
#define OM_gam0 2 * pow(pi, 2) / 30 * pow(T_CMB0, 4) * 8 * pi / (3 * pow(H_0, 2))  // Present radiation density parameter
#define OM_mu0 N_eff * 7 / 8 * pow(4 / 11, 4 / 3) * OM_gam0                        // Present nutrino density parameter
#define OM_mu0 0      // Present nutrino density parameter
#define OM_lam0 1 - (OM_k0 + OM_b0 + OM_CMB0 + OM_gam0 + OM_mu0)                   // Present cosmological constant density parameter

using namespace std;

double gammln(double);


void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i, its, j;
    double ai;
    double p1, p2, p3, pp, z, z1;

    for (i = 1; i <= n; i++)
    {
        if (i == 1)
        {
            z = (1.0 + alf) * (3.0 + 0.92 * alf) / (1.0 + 2.4 * n + 1.8 * alf);
        }
        else if (i == 2)
        {
            z += (15.0 + 6.25 * alf) / (1.0 + 0.9 * alf + 2.5 * n);
        }
        else
        {
            ai = i - 2;
            z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alf /
                                                       (1.0 + 3.5 * ai)) *
                 (z - x[i - 2]) / (1.0 + 0.3 * alf);
        }
        for (its = 1; its <= MAXIT; its++)
        {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 1; j <= n; j++)
            {
                p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1 + alf - z) * p2 - (j - 1 + alf) * p3) / j;
            }
            pp = (n * p1 - (n + alf) * p2) / z;
            z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) <= EPS)
                break;
        }
        if (its > MAXIT)
            cout << "too many iterations in gaulag" << endl;
        x[i] = z;
        w[i] = -exp(gammln(alf + n) - gammln((double)n)) / (pp * n * p2);
    }
}
// end function gaulag

#undef EPS
#undef MAXIT