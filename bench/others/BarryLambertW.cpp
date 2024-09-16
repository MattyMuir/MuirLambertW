#include "BarryLambertW.h"

#include <cmath>
#include <iostream>

double BarryLambertW(int64_t branch, double x)
{
    double an2;
    static double an3;
    static double an4;
    static double an5;
    static double an6;
    static double c13;
    static double c23;
    static double d12;
    double delx;
    static double em;
    static double em2;
    static double em9;
    double eta;
    int i;
    static int init = 0;
    static int nbits;
    static int niter = 1;
    double reta;
    static double s2;
    static double s21;
    static double s22;
    static double s23;
    double t;
    static double tb;
    double temp;
    double temp2;
    double ts;
    double value;
    static double x0;
    static double x1;
    double xx;
    double zl;
    double zn;

    value = 0.0;

    if (init == 0)
    {
        init = 1;
        nbits = 52;
        //
        //  Various mathematical constants.
        //
        em = -exp(-1.0);
        em9 = -exp(-9.0);
        c13 = 1.0 / 3.0;
        c23 = 2.0 * c13;
        em2 = 2.0 / em;
        d12 = -em2;
        tb = pow(0.5, nbits);
        x0 = pow(tb, 1.0 / 6.0) * 0.5;
        x1 = (1.0 - 17.0 * pow(tb, 2.0 / 7.0)) * em;
        an3 = 8.0 / 3.0;
        an4 = 135.0 / 83.0;
        an5 = 166.0 / 39.0;
        an6 = 3167.0 / 3549.0;
        s2 = sqrt(2.0);
        s21 = 2.0 * s2 - 3.0;
        s22 = 4.0 - 3.0 * s2;
        s23 = s2 - 2.0;
    }

    if (x == em)
    {
        value = -1.0;
        return value;
    }
    xx = x;
    delx = xx - em;
    //
    //  Calculations for Wp.
    //
    if (branch == 0)
    {
        if (fabs(xx) <= x0)
        {
            value = xx / (1.0 + xx / (1.0 + xx
                / (2.0 + xx / (0.6 + 0.34 * xx))));
            return value;
        }
        else if (xx <= x1)
        {
            reta = sqrt(d12 * delx);
            value = reta / (1.0 + reta / (3.0 + reta / (reta
                / (an4 + reta / (reta * an6 + an5)) + an3)))
                - 1.0;
            return value;
        }
        else if (xx <= 20.0)
        {
            reta = s2 * sqrt(1.0 - xx / em);
            an2 = 4.612634277343749 * sqrt(sqrt(reta +
                1.09556884765625));
            value = reta / (1.0 + reta / (3.0 + (s21 * an2
                + s22) * reta / (s23 * (an2 + reta)))) - 1.0;
        }
        else
        {
            zl = log(xx);
            value = log(xx / log(xx
                / pow(zl, exp(-1.124491989777808 /
                    (0.4225028202459761 + zl)))));
        }
    }
    //
    //  Calculations for Wm.
    //
    else
    {
        if (xx <= x1)
        {
            reta = sqrt(d12 * delx);
            value = reta / (reta / (3.0 + reta / (reta / (an4
                + reta / (reta * an6 - an5)) - an3)) - 1.0) - 1.0;
            return value;
        }
        else if (xx <= em9)
        {
            zl = log(-xx);
            t = -1.0 - zl;
            ts = sqrt(t);
            value = zl - (2.0 * ts) / (s2 + (c13 - t
                / (270.0 + ts * 127.0471381349219)) * ts);
        }
        else
        {
            zl = log(-xx);
            eta = 2.0 - em2 * xx;
            value = log(xx / log(-xx / ((1.0
                - 0.5043921323068457 * (zl + 1.0))
                * (sqrt(eta) + eta / 3.0) + 1.0)));
        }

    }

    for (i = 1; i <= niter; i++)
    {
        zn = log(xx / value) - value;
        temp = 1.0 + value;
        temp2 = temp + c23 * zn;
        temp2 = 2.0 * temp * temp2;
        value = value * (1.0 + (zn / temp) * (temp2 - zn)
            / (temp2 - 2.0 * zn));
    }

    return value;
}

double BarryLambertW0(double x)
{
    return BarryLambertW(0, x);
}

double BarryLambertWm1(double x)
{
    return BarryLambertW(-1, x);
}