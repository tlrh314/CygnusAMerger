/* Numerical Recipes chapter 3 section 1: searching an ordered table */

#include "globals.h"
#define SPLINETABLE 2048

// move to proto.h
int locate(const double);
int hunt(const double);
double rawinterp_linear(int, double);
double rawinterp_spline(int, double, double, double);
double interp(double);
void set_y2(const double *xv, const double *yv, double yp1, double ypn);

// leave here
extern struct Spline_interp
{
    int n, mm, jsav, cor, dj;
    const double *xx, *yy;
    // NB only for Spline interpolation, not for linear
    double *y2;
    bool debug;
    bool y2_set;
    // int n = sizeof( xx ) / sizeof( xx[0] );
} Spline;

// move to aux.c
struct Spline_interp Spline;

// leave here
int hunt(const double x)
{
    int jl=Spline.jsav, jm, ju, inc = 1;
    if (Spline.debug)
        printf("%d %d %d %d\n", jl, jm, ju, inc);

    if (Spline.n < 2 || Spline.mm < 2 || Spline.mm > Spline.n) {
        printf("Hunt size error\n");
        return Spline.n+1;  // to ensure runtime error will be thrown
    }

    // bool: True if ascending order of table, else False
    bool ascnd = (Spline.xx[Spline.n-1] >= Spline.xx[0]);
    if (Spline.debug)
        printf("ascnd=%d\n", ascnd);

    // input guess not useful, use bisection
    if (jl < 0 || jl > Spline.n-1) {
        jl = 0;
        ju = Spline.n-1;
        if (Spline.debug)
            printf("hunt input not useful\n");
    } else {
        if (x >= Spline.xx[jl] == ascnd) {  // hunt up
            if (Spline.debug)
                printf("hunting up\n");
            for(;;) {
                ju = jl + inc;
                if (Spline.debug)
                    printf("ju=%d\n", ju);
                if (ju >= Spline.n-1)  // outside end of table
                    { ju = Spline.n-1; break; }
                else if(x < Spline.xx[ju] == ascnd) break;  // found bracket
                else {
                    jl = ju;
                    inc += inc;
                }
            }
        } else {  // hunt down
            if (Spline.debug)
                printf("hunting down\n");
            ju = jl;
            for(;;) {
                jl = jl - inc;
                if (jl <= 0) { jl = 0; break; } // outside begin of table
                else if (x >= Spline.xx[jl] == ascnd) { break; } // found bracket
                else {  // Not done, so double increment and try again.
                    ju = jl;
                    inc += inc;
                }
            }
        }
    }  // Done hunting. Now finish using bisection
    if(Spline.debug)
        printf("Now bisection\nju=%d\njl=%d\n", ju, jl);
    while (ju - jl > 1) {
        jm = (ju+jl) >> 1;
        if (Spline.debug) printf("jm =%d\n", jm);
        if (x >= Spline.xx[jm] == ascnd)
            jl = jm;
        else
            ju = jm;
    }

    // Decide whether to hunt or locate next time
    Spline.cor = abs(jl-Spline.jsav) > Spline.dj ? 0 : 1;
    Spline.jsav = jl;
    if (Spline.debug)
        printf("cor=%d, jsav=%d\n", Spline.cor, Spline.jsav);

    return max(0, min(Spline.n-Spline.mm, jl-((Spline.mm-2)>>1)));
}

// leave here
int locate(const double x)
{
    int ju, jm, jl;
    if(Spline.n < 2 || Spline.mm < 2 || Spline.mm > Spline.n) {
        printf("Locate size error\n");
        return Spline.n+1;  // to ensure runtime error will be thrown
    }

    // bool: True if ascending order of table, else False
    bool ascnd = (Spline.xx[Spline.n-1] >= Spline.xx[0]);
    jl = 0;  // lower bound
    ju = Spline.n-1;  // upper bound
    if(Spline.debug) {
        printf("n = %d\n", Spline.n);
        printf("x[n-1] = %f\n", Spline.xx[Spline.n-1]);
        printf("x[0] = %f\n", Spline.xx[0]);
        printf("ascnd = %d\n", ascnd);
    }
    while (ju - jl > 1) {
        jm = (ju+jl) >> 1;  // compute midpoint
        if (Spline.debug) printf("jm =%d\n", jm);
        // replace lower/upper bound
        if (x >= Spline.xx[jm] == ascnd)
            jl = jm;
        else
            ju = jm;
    }
    // decide whether to use hunt or locate next time
    Spline.cor = abs(jl-Spline.jsav) > Spline.dj ? 0 : 1;
    Spline.jsav = jl;
    if (Spline.debug)
        printf("cor=%d, jsav=%d\n", Spline.cor, Spline.jsav);

    return max(0, min(Spline.n-Spline.mm, jl-((Spline.mm-2)>>1)));
}

double rawinterp_linear(int j, double x)
{
    if(Spline.xx[j] == Spline.xx[j+1]) {
        printf("Table is defective, but we can recover\n");
        return Spline.yy[j];
    } else {
        return Spline.yy[j] + ((x-Spline.xx[j])/(Spline.xx[j+1]-Spline.xx[j]))*(Spline.yy[j+1]-Spline.yy[j]);
    }
}

void set_y2(const double *xv, const double *yv, double yp1, double ypn)
{
    int i, k;
    double p, qn, sig, un;
    double u[SPLINETABLE-1] = { 0 };

    // The lower boundary condition is set either to be natural,
    //     or to have specified first derivative.
    if (yp1 > 0.99e99)
        Spline.y2[0] = u[0] = 0.0;
    else {
        Spline.y2[0] = -0.5;
        u[0] = (3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1);
    }

    // This is the decomposition loop of the tridiagonal algorithm.
    // y2 and u are used for temporary storage of the decomposed factors.
    for (i=1; i<Spline.n-1; i++) {
        sig = (xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]);
        p = sig*Spline.y2[i-1] + 2.0;
        Spline.y2[i] = (sig-1.0)/p;
        u[i] = (yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]);
        u[i] = (6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p;
    }

    // The upper boundary condition is set either to be natural,
    //     or to have a specified first derivative.
    if (ypn > 0.99e99)
        qn = un = 0.0;
    else {
        qn = 0.5;
        un = (3.0/(xv[Spline.n-1]-xv[Spline.n-2]))*
            (ypn-(yv[Spline.n-1]-yv[Spline.n-2])/(xv[Spline.n-1]-xv[Spline.n-2]));
    }
    Spline.y2[Spline.n-1] = (un-qn*u[Spline.n-2])/(qn*Spline.y2[Spline.n-2]+1.0);
    // This is the backsubstitution loop for the tridiagonal algorithm.
    for (k=Spline.n-2; k>=0; k--) {
        Spline.y2[k] = Spline.y2[k]*Spline.y2[k+1]+u[k];
    }

    Spline.y2_set = true;
}

double rawinterp_spline(int jl, double x, double yp1, double ypn)
{
    set_y2(&Spline.xx[0], &Spline.yy[0], yp1, ypn);

    int klo = jl, khi = jl+1;
    double y, h, b, a;

    h = Spline.xx[khi] - Spline.xx[klo];
    if (h == 0.0) { printf("Bad input to rawinterp_spline\n"); exit(0); }
    a = (Spline.xx[khi]-x)/h;
    b = (x - Spline.xx[klo])/h;

    if (Spline.debug)
        printf("klo=%d\nkhi=%d\nh=%d\na=%g\nb=%g\n", klo, khi, h, a, b);

    y = a*Spline.yy[klo]+b*Spline.yy[khi] +
        ((a*a*a-a)*Spline.y2[klo]+(b*b*b-b)*Spline.y2[khi])*(h*h/6.0);
    return y;
}

double interp(double x)
{
    int jlo = Spline.cor ? hunt(x) : locate(x);
    if (Spline.debug) printf("jlo=%d\n", jlo);
    // return rawinterp_linear(jlo, x);
    return rawinterp_spline(jlo, x, 1.e99, 1.e99);
}

double rawintegrate_spline(int jl, double x)
{
    if (! Spline.y2_set)
        set_y2(&Spline.xx[0], &Spline.yy[0], 1.e99, 1.e99);

    int klo = jl, khi = jl+1;
    double integrated_y, h;

    h = Spline.xx[khi] - Spline.xx[klo];
    if (h == 0.0) { printf("Bad input to routine splint\n"); exit(0); }

    if (Spline.debug)
        printf("klo=%d\nkhi=%d\nh=%d\n", klo, khi, h);

    integrated_y = h/2*(Spline.yy[klo]+Spline.yy[khi])
        - p3(h)/24*(Spline.y2[klo]+Spline.y2[khi]);
    return integrated_y;
}

double integrate(double x)
{
    int jlo = Spline.cor ? hunt(x) : locate(x);
    if (Spline.debug) printf("jlo=%d\n", jlo);
    return rawintegrate_spline(jlo, x);
}

// lives in setup.c
double Gas_density_profile_freebeta(const double r, const double rho0,
        const double rc, const double rcut, const double beta, const bool Is_Cuspy)
{
    // NB Donnert 2014 w/o additional cut-off at r200 for numerical stability!
    // Not Donnert+ 2016 in prep, where rho/=(1+p4(r/rcut))
    double rho = rho0 * pow((1 + p2(r/rc)), -3.0/2*beta) /
        (1 + p3(r/rcut) * (r/rcut));

    return rho;
}


int main(int argc, char *argv[])
{
    Spline.n = SPLINETABLE;
    Spline.mm = 2;  // For Linear and Spline interpolation
    Spline.jsav = 0;
    Spline.cor = 0;
    Spline.dj = min(1, (int)pow((double)Spline.n, 0.25));
    Spline.debug = false;
    Spline.y2_set = false;

    if (Spline.debug)
        printf("n=%d, mm=%d, jsav=%d, cor=%d, dj=%d\n\n", Spline.n, Spline.mm,
                Spline.jsav, Spline.cor, Spline.dj);


    // move to setup
    double xx[SPLINETABLE] = { 0 };
    double yy[SPLINETABLE] = { 0 };
    double y2[SPLINETABLE] = { 0 };

    double r = 0.0, r_lower = 1.0, r_upper = 5000.0;
    for (int i=1; i<SPLINETABLE; i++)
    {
        // r = r_lower + i*(r_upper-r_lower)/(SPLINETABLE-1);
        r = pow(10, log10(r_lower)+i*(log10(r_upper)-log10(r_lower))/(SPLINETABLE-1));
        xx[i] = r;
        yy[i] = 4*pi*1.0*p2(r)*Gas_density_profile_freebeta(r, 1, 25, 1200, 0.5, 1);
        printf("r=%g, rho=%g\n", xx[i], yy[i]);
    }

    Spline.xx = xx;
    Spline.yy = yy;
    Spline.y2 = y2;

    printf("END OF TABLE\n");

    r = 0.0, r_lower = 1.0, r_upper = 5000.0;
    double integrated = 0.0;
    for (int i=1; i<16*SPLINETABLE; i++)
    {
        r = pow(10, log10(r_lower)+i*(log10(r_upper)-log10(r_lower))/(16*SPLINETABLE-1));
        integrated += integrate(r);
        printf("r=%g, spline=%g\n", r, integrated/16);
    }


    // printf("\ninterp(41) = %g\n", interp(41.5));

#undef SPLINETABLE
    return 0;
}
