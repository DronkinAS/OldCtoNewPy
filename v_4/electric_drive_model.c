#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>


void converter_U(double A, double *U, double angle)
{
    double Ua, Ub, Uc;

    Ua = A * sin(angle);
    Ub = A * sin(angle + 4*M_PI/3);
    Uc = A * sin(angle + 2*M_PI/3);

    *U = Ua*cos(angle) + ((Ub - Uc)/sqrt(3.0))*sin(angle);
    *(U + 1) = ((Ub - Uc)/sqrt(3.0))*cos(angle) - Ua*sin(angle);
}

double calculation_Id(double L, double R, double p,
                      double U, double *I, double w)
{
    double Id, Iq, out;

    Id = *I;
    Iq = *(I + 1);
    out = (-1*R*Id)/L + p*w*Iq + U/L;

    return out;
}


double calculation_Iq(double L, double R, double p, double F,
                      double U, double *I, double w)
{
    double Id, Iq, out;

    Id = *I;
    Iq = *(I + 1);
    out = (-1*R*Iq)/L - p*w*Id + U/L - (p*w*F)/L;

    return out;
}

double calculation_Mn(double Mn_0, double Mn_1, double Mn_2,
                      double w, double kr, double nr)
{
    double Mn;

    Mn = Mn_1 + Mn_2;
    Mn = Mn/(kr*nr);
    Mn = Mn + Mn_0;

    if (w < 0)
    {
        Mn = -1*Mn;
    }

    return Mn;
}

double calculation_J(double J_0, double J_1, double J_2,
                     double kr, double nr)
{
    double J;

    J = J_1 + J_2;
    J = J/(kr*kr*nr);
    J = J + J_0;

    return J;
}

double PID(double *error, double dt)
{
    double Td = 0.8, Kp = 100.0, Ti = 0.01, signal;
    
    signal = *(error + 1) * Kp;
    signal += Td * ( *(error + 1) - *(error) ) / dt;
    signal += ( *error + *(error + 1)* dt ) / Ti;
    
    return signal;
}