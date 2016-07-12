#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>


/* GNU Scientifc Library */
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_rng.h>


#define p2(a) ((a)*(a))
#define p3(a) ((a)*(a)*(a))

/* mathematical constants */
#define pi 			M_PI
#define sqrt2		M_SQRT2
#define sqrt3       1.73205080756887719
#define fourpithird 4.18879032135009765

/* physical constants cgs */
#define c			GSL_CONST_CGSM_SPEED_OF_LIGHT
#define k_B 		GSL_CONST_CGSM_BOLTZMANN
#define m_p 		GSL_CONST_CGSM_MASS_PROTON
#define m_e			GSL_CONST_CGSM_MASS_ELECTRON
#define Grav        GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT

/* unit conversions */
#define Msol2cgs    (1.98892e33)
#define kpc2cgs 	(3.08568025e21)
#define K2eV        (1.5*8.617343e-5)
#define DEG2RAD		(pi / 180)


double extern gsl_sf_hyperg_2F1(double, double, double, double);
double extern gsl_sf_gamma(double);


#if defined(FREEBETA)
/* return M(<= R) of a beta profile with beta=free parameter */
double Mass_profile(const double r, const double rho0, const double rc,
        const double rcut, const double beta, const bool Is_Cuspy)
{
    /*
     * M(<r) = int rho(r) 4 pi r**2 dr. Beautifully analtical solution
     * onlyif beta == 2/3 or beta == -1. If beta is a free parameter
     * we get a lovely hypergeometrical function 2F1(a,b;c;x)
     *
     * These are built-in in gsl, and in scipy.special using cephes library.
     * In scipy this works, but both in gsl and when using cephes it breaks.
     *
     * gsl_sf_hyperg_2F1 only works for x in [-1, 1]
     * and, sadly, hyp2f1 in cephes math library breaks with overflow error
     *
     * Therefore we use some dark magic gleaned from
     * https://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf
     * Specifially table 13 case 1, and the corresponding equation 4.16 works
     * Lets pray to the gods of mathemagics that this voodoo works :-)...
     */

    // Probably a mighty bad idea to shut the gsl error handler up ...
    // gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    // yay, this is no longer needed because it works :-)!

    // gsl_set_error_handler(old_handler);

    // Because Wolfram Alpha says this is the solution to integrating
    // int r**2(1+r**2/a**2)**(-3/2 * beta) dr
    //      = 1/3 r**3 2F1(1.5, 1.5*beta; 2.5; -r**2/rc**2) + constant
    // TODO: int r**2(1+r**2/a**2)**(-3/2 * beta)*(1+r**4/b**4)**-1 dr
    // WA: "no result found in terms of standard mathematical functions"
    double a = 1.5;
    double b = 1.5*beta;
    // cannot be called c because c is defined as speed of light!
    double c_hg = 2.5;
    double z = -p2(r)/p2(rc);
    //  Probably a mighty bad idea ...
    // gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double magic = NAN;


    // Because black magic voodoo ArXiV 1502.05624 equation 5, 6
    // TODO: add original reference referenced by ref above..
    if (z < -2.0)
    {
        magic = pow(-z, -a)
            * (gsl_sf_gamma(c_hg)*gsl_sf_gamma(b-a))
            / (gsl_sf_gamma(b)*gsl_sf_gamma(c_hg-a))
            * gsl_sf_hyperg_2F1(a, 1.0-c_hg+a, 1.0-b+a, 1.0/z)
            + pow(-z, -b)
            * (gsl_sf_gamma(c_hg)*gsl_sf_gamma(a-b))
            / (gsl_sf_gamma(a)*gsl_sf_gamma(c_hg-b))
            * gsl_sf_hyperg_2F1(b, 1.0-c_hg+b, 1.0-a+b, 1.0/z);
    }
    // Because black magic voodoo ArXiV 1502.05624 equation 3, 4
    // TODO: add original reference referenced by ref above..
    else if (z < -1.0)
    {
        // printf("1-z=%g\n", 1.0-z);
        // only if |arg(1-z)| < pi. Solved by the z < -2 above
        magic = pow(1.0-z, -a)
            * (gsl_sf_gamma(c_hg)*gsl_sf_gamma(b-a))
            / (gsl_sf_gamma(b)*gsl_sf_gamma(c_hg-a))
            * gsl_sf_hyperg_2F1(a, c_hg-b, a-b+1.0, 1.0/(1.0-z))
            + pow(1.0-z, -b)
            * (gsl_sf_gamma(c_hg)*gsl_sf_gamma(a-b))
            / (gsl_sf_gamma(a)*gsl_sf_gamma(c_hg-b))
            * gsl_sf_hyperg_2F1(b, c_hg-a, b-a+1.0, 1.0/(1.0-z));
    }
    else
    {
        magic = gsl_sf_hyperg_2F1(a, b, c_hg, z);
    }

    // printf("magic=%g\n", magic);
    double Mr = 4 * pi * rho0 * p3(r)/3 * magic;
    // printf("Mr=%g\n\n", Mr*g2msun);

    // gsl_set_error_handler(old_handler);

    return Mr;
}
#else
/* return M(<= R) of a beta profile with beta=2/3 */
double Mass_profile(const double r, const double rho0, const double rc,
        const double rcut, const bool Is_Cuspy)
{
    // beta == 2/3
    const double r2 = p2(r);
    const double rc2 = p2(rc);
    const double rcut2 = p2(rcut);

    double Mr = rho0 * rc2*rcut2*rcut/(8*(p2(rcut2)+p2(rc2))) * // fourth order
        ( sqrt2 *( (rc2-rcut2) *( log(rcut2 - sqrt2*rcut*r+r2)
                                - log(rcut2 + sqrt2*rcut*r+r2))
                   - 2 * (rc2 + rcut2) * atan(1 - sqrt2*r/rcut)
                   + 2 * (rc2 + rcut2) * atan(sqrt2 * r/rcut + 1))
          - 8 * rc * rcut * atan(r/rc));

    return 4 * pi * Mr;
}
#endif // FREEBETA


/* Double beta profile at rc and rcut */
double Gas_density_profile(const double r, const double rho0,
        const double rc, const double rcut, const double beta, const bool Is_Cuspy)
{
    // NB Donnert 2014 w/o additional cut-off at r200 for numerical stability!
    // Not Donnert+ 2016 in prep, where rho/=(1+p4(r/rcut))
    double rho = rho0 * pow((1 + p2(r/rc)), -3.0/2*beta) / (1 + p3(r/rcut) * (r/rcut));

    return rho;
}

#define TABLESIZE 1024

static gsl_spline *U_Spline = NULL;
static gsl_interp_accel *U_Acc = NULL;
#pragma omp threadprivate(U_Spline, U_Acc)

int main(int argc, char *argv[])
{

    // #pragma omp parallel for
    for (int r = 0; r < 2000; r++) {
        printf("r=%d, rho=%g\n", r,
                Gas_density_profile(r, 1, 25, 750, 2.0/3, 0));
    }



    return 0;
}
