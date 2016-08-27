import numpy
import matplotlib
from matplotlib import pyplot
matplotlib.use("Qt4Agg", warn=False)

def p2(a):
    return (a)*(a)

def p3(a):
    return (a)*(a)*(a)

def rho(r, rho0, rc, rcut, beta):
    return rho0 * (1+p2(r/rc))**(-1.5*beta) * (1+p2(r/rcut)*(r/rcut))**(-1)

def mass(r, rho0, rc, rcut):
    r2 = p2(r)
    rc2 = p2(rc)
    rcut2 = p2(rcut)
    sqrt2 = numpy.sqrt(2)

    A = (rc2 - rcut2)*(numpy.log(rcut2 - sqrt2*rcut*r + r2) \
        - numpy.log(rcut2 + sqrt2*rcut*r + r2))
    Bplus = 2*(rc2 + rcut2)*numpy.arctan(1 + sqrt2*r/rcut)
    Bmin = 2*(rc2 + rcut2)*numpy.arctan(1 - sqrt2*r/rcut)

    M_gas_below_r = (4 * numpy.pi * rho0) *\
        rc2*p3(rcut)/(8*(p2(rcut2)+p2(rc2))) *\
        (sqrt2*(A - Bmin + Bplus) - 8*rc*rcut*numpy.arctan(r/rc))
    return M_gas_below_r


class Base_interp(object):
    """ Numerical Recipes (Press+ 2007) Chapter 3 """
    def __init__(self, x, y, m):
        self.debug = False

        # Tables
        self.xx = x
        self.yy = y
        if self.debug:
            print self.xx
            print self.yy

        self.n = len(x)
        self.mm = m
        self.jsav = 0
        self.cor = False

        self.dj = numpy.minimum(1, int(numpy.power(self.n, 0.25)))
        if self.debug: print "dj =", self.dj

    def interp(self, x):
        jlo = self.hunt(x) if self.cor else self.locate(x)
        if self.debug: print "jlo =", jlo
        return self.rawinterp(jlo, x)

    def hunt(self, x):
        jl = self.jsav
        jm, jl = 0, 0
        inc = 1

        if self.n < 2 or self.mm < 2 or self.mm > self.n:
            raise ValueError("Hunt size error")

        # bool: True if ascending order of table, else False
        ascnd = (self.xx[self.n-1] >= self.xx[0])
        if self.debug: print "ascnd =", ascnd

        if jl < 0 or jl > self.n-1:  # input guess not useful, use bisection
            jl = 0
            ju = n-1
            if self.debug: print "hunt input not useful"
        else:
            if x >= self.xx[jl] and ascnd:  # hunt up
                if self.debug: print "hunting up"
                while True:
                    ju = jl + inc
                    if self.debug: print "ju =", ju
                    if ju >= self.n-1:  # outside end of table
                        ju = self.n-1
                        break
                    elif x < self.xx[ju] and ascnd:  # found bracket
                        break
                    else:
                        jl = ju
                        inc += inc
            else:  # hunt down
                if self.debug: print "hunting down"
                ju = jl
                while True:
                    jl = jl - inc
                    if jl <= 0:  # outside begin of table
                        jl = 0
                        break
                    elif x >= self.xx[jl] and ascnd:
                        break
                    else:
                        ju = jl
                        inc += inc
        # Done hunting. Now finish using bisection
        if self.debug:
            print "Now bisection"
            print "ju =", ju
            print "jl =", jl
        while ju - jl > 1:
            jm = (ju+jl >> 1)
            # print "jl =", jl
            if self.debug: print "jm =", jm
            # print "ju =", ju
            if x >= self.xx[jm] and ascnd:
                jl = jm
            else:
                ju = jm

        self.cor = True if numpy.absolute(jl - self.jsav) > self.dj else False
        self.jsav = jl
        if self.debug:
            print "cor =", self.cor
            print "jsav = ", self.jsav

        return numpy.maximum(0, numpy.minimum(self.n-self.mm, (jl-(self.mm-2)>>1)))

    def locate(self, x):
        if self.n < 2 or self.mm < 2 or self.mm > self.n:
            raise ValueError("Locate size error")

        # bool: True if ascending order of table, else False
        ascnd = (self.xx[self.n-1] >= self.xx[0])
        if self.debug:
            print "n =", self.n
            print "x[n-1] =", self.xx[self.n-1]
            print "x[0] =", self.xx[0]
            print "ascnd =", ascnd
        # lower bound
        jl = 0
        # upper bound
        ju = self.n-1
        while (ju - jl > 1):
            jm = (ju+jl) >> 1  # compute midpoint
            if self.debug:
                # print "jl = ", jl
                print "jm =", jm
                # print "ju = ", ju
            # replace lower/upper bound
            if x >= self.xx[jm] and ascnd:
                jl = jm
                if self.debug:
                    # print "setting lower bound to current bound"
                    pass
            else:
                ju=jm
                if self.debug:
                    # print "setting upper bound to current bound"
                    pass
        # decide whether to use hunt or locate next time
        self.cor = True if numpy.absolute(jl - self.jsav) > self.dj else False
        self.jsav = jl
        if self.debug:
            print "cor =", self.cor
            print "jsav = ", self.jsav

        return numpy.maximum(0, numpy.minimum(self.n-self.mm, jl-((self.mm-2)>>1)))

    def rawinterp(self, jlo, x):
        print "Not implemented: Placeholder / constructur / template..."
        pass


class Linear_interp(Base_interp):
    def __init__(self, xv, yv):
        super(Linear_interp, self).__init__(xv, yv, 2)

    def rawinterp(self, j, x):
        if self.xx[j] == self.xx[j+1]:
            print "Table is defective, but we can recover"
            return self.yy[j]
        else:
            return self.yy[j] + ((x-self.xx[j])/(self.xx[j+1]-self.xx[j]))*(self.yy[j+1]-self.yy[j])

class Spline_interp(Base_interp):
    def __init__(self, xv, yv, yp1=1.e99, ypn=1.e99):
        """
        xv : table with x values
        yv : table with y values

        if yp1 and ypn are equal to 1e99 or larger, the routine is signaled
        to set the corresponding boundary condition for a natural spline,
        with zero second derivative on that boundary; otherwise, they are
        the values of the first derivatives at the endpoints.
        """

        super(Spline_interp, self).__init__(xv, yv, 2)
        self.y2 = numpy.zeros(self.n)
        self.sety2(self.xx, self.yy, yp1, ypn)

    def sety2(self, xv, yv, yp1, ypn):
        """ This routine stores an array y2[0..n-1] with second derivatives
        of the interpolating function at the tabulated points pointed to by
        xv, using function values pointed to by yv. If yp1 and/or ypn are equal
        to 1e99 or larger, the routine is signaled to set the corresponding
        boundary condition for a natural spline, with zero second derivative on
        that boundary; otherwise, they are the values of the first derivatives
        at the endpoints. """

        u = numpy.zeros(self.n-1)

        # The lower boundary condition is set either to be natural,
            # or to have specified first derivative.
        if yp1 > 0.99e99:
            self.y2[0] = 0.0
            u[0] = 0.0
        else:
            self.y2[0] = -0.5
            u[0] = (3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1)

        # This is the decomposition loop of the tridiagonal algorithm.
        # y2 and u are used for temporary storage of the decomposed factors.
        for i in range(1, self.n-1):
            sig = (xv[i]-xv[i-1])/(xv[i+1]-xv[i-1])
            p = sig*self.y2[i-1] + 2.0
            self.y2[i] = (sig-1.0)/p
            u[i] = (yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1])
            u[i] = (6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p

        # The upper boundary condition is set either to be natural,
            # or to have a specified first derivative.
        if ypn > 0.99e99:
            qn = 0.0
            un = 0.0
        else:
            qn = 0.5
            un = (3.0/(xv[self.n-1]-xv[self.n-2]))*\
                (ypn-(yv[self.n-1]-yv[self.n-2])/(xv[self.n-1]-xv[self.n-2]))
        self.y2[self.n-1] = (un-qn*u[self.n-2])/(qn*self.y2[self.n-2]+1.0)
        # This is the backsubstitution loop for the tridiagonal algorithm.
        for k in range(self.n-2, -1, -1):
            self.y2[k] = self.y2[k]*self.y2[k+1]+u[k]

    def rawinterp(self, jl, x):
        klo = jl
        khi = jl+1

        h = self.xx[khi]-self.xx[klo]
        if h == 0.0:  # does this work in Python?
            raise ValueError("Bad input to routine splint")
        a = (self.xx[khi]-x)/h
        b = (x - self.xx[klo])/h

        if self.debug: print "klo =", klo, "\nkhi =", khi, "\nh =", h, "\na =", a, "\nb =", b

        y = a*self.yy[klo]+b*self.yy[khi] + ((a*a*a-a)*self.y2[klo]
                +(b*b*b-b)*self.y2[khi])*(h*h/6.0)

        return y


if __name__ == "__main__":
    pyplot.figure(figsize=(12,9))

    radii = numpy.arange(-10, 12, 2)
    spline = Spline_interp(radii, p2(radii))
    yvals_num = numpy.zeros(len(radii))
    yvals_int = numpy.zeros(len(radii))

    for i, r in enumerate(radii):
        yvals_num[i] = spline.interp(r)
        if i < len(radii)-1:
            h = spline.xx[i+1] - spline.xx[i]
            yvals_int[i] = -h * (0.5*spline.yy[i] + 0.5*spline.yy[i+1] - p2(h)/24 * spline.y2[i] - p2(h)/24 * spline.y2[i+1])

    pyplot.plot(radii, yvals_num, label="spline")
    pyplot.plot(radii, yvals_int, label="spline integrated")
    radii = numpy.arange(-10, 10, 0.001)
    pyplot.plot(radii, p2(radii), label="analytical")
    pyplot.plot(radii, p3(radii)/3, label="analytical integrated")

    pyplot.legend(loc="best")
    pyplot.show()
    import sys; sys.exit(0)



    # table
    radii = numpy.arange(0, 4000, 0.1)
    density = rho(radii, 1, 30.0, 1200.0, 1/1.5)

    mass_numerical = 4*numpy.pi*p2(radii)*rho(radii, 1, 30.0, 1200.0, 1/1.5)
    mass_analytical = mass(radii, 1, 30.0, 1200.0)

    spline = Spline_interp(radii, mass_numerical)
    spline.debug = False

    pyplot.figure(figsize=(12,9))
    pyplot.plot(radii, mass_analytical, label="M(<r)")

    spline.integrate(1, 4000, logbins=False)

    # radii = numpy.arange(0, 4000, 0.1)
    # linear_interp = numpy.zeros(len(radii))
    # spline_interp = numpy.zeros(len(radii))
    # for i, r in enumerate(radii):
    #     linear_interp[i] = linear.interp(r)
    #     spline_interp[i] = spline.interp(r)

    pyplot.figure(figsize=(12,9))
    pyplot.plot(xvals, yvals, label="table")
    # pyplot.plot(radii, linear_interp, label="linear interp")
    # pyplot.plot(radii, spline_interp, label="spline interp")
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.xlim(1, 4000)
    #pyplot.ylim(1e-6, 2)
    pyplot.legend(loc="best")
    pyplot.show()
