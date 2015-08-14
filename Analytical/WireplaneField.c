/// \file WireplaneField.c
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "WireplaneField.h"
#include <math.h>

void wireplaneField(double r, double d, double l, double a, double* Eperp, double* Eparr) {
    a = fmod(a,d);
    if(a*a + l*l < r*r) { *Eperp = *Eparr = 0; return; }
    
    double c1 = 2*M_PI*l/d;
    double c2 = 2*M_PI*a/d;
    double denom = 1 - cos(c2)/cosh(c1);
    *Eperp = tanh(c1)/denom;
    *Eparr = sin(c2)/cosh(c1)/denom;
}

double wireplaneVOffset(double r, double d) {
    return d/M_PI * log(d/(2*M_PI*r));
}

double wireplanePotential(double r, double d, double l, double a) {
    a = fmod(a,d);
    if(a*a + l*l < r*r) return 0;
    
    double c1 = M_PI*l/d;
    if(c1 > 20) return l + wireplaneVOffset(r,d);       // asymptotic limit
    
    double V = d/M_PI * log(sinh(c1) / sinh(M_PI*r/d)); // integrate Eperp at constant a=0
    double csh = cosh(2*c1);
    double c2 = 2*M_PI*a/d;
    V += d/(2*M_PI) * log((csh - cos(c2))/(csh - 1));   // integrate Eparr at constant l
    return V;
}
