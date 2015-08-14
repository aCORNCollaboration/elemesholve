/// \file EndBesselCalc.c
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "EndBesselCalc.h"
#include <gsl/gsl_sf_bessel.h>
#include <math.h>
#include <assert.h>

void init_j0n() {
    unsigned int n;
    bessel_j0n[0] = 0;
    for(n=1; n<MAX_BESSEL_TERMS; n++) bessel_j0n[n] = gsl_sf_bessel_zero_J0(n);
}

double flatcircle_coeff(unsigned int n, double r) {
    assert(n < MAX_BESSEL_TERMS);
    double j1zn = gsl_sf_bessel_J1(bessel_j0n[n]);
    if(r == 1) return 2./(bessel_j0n[n]*j1zn);
    return 2*gsl_sf_bessel_J1(r*bessel_j0n[n])/(bessel_j0n[n]*j1zn*j1zn);
}

void addCircle(double coeffs[MAX_BESSEL_TERMS], double r, double V) {
    unsigned int n;
    for(n=1; n<MAX_BESSEL_TERMS; n++) coeffs[n] += V*flatcircle_coeff(n,r);
}

/// truncated coefficient sum mechanism
double coeffSum(const double coeffs[MAX_BESSEL_TERMS], double (*f)(unsigned int, void*), void* params, double relErr, double absErr) {
    double s = 0;
    int ntiny = 0;
    unsigned int n;
    for(n=1; n<MAX_BESSEL_TERMS; n++) {
        double ds = coeffs[n]*(*f)(n,params);
        s += ds;
        if(fabs(ds) < fabs(relErr*s) || fabs(ds) < absErr) ntiny++;
        else ntiny = 0;
        if(ntiny > 5) break;
    }
    return s;
}

double _sumBessel(unsigned int n, void* zr) {
    double z = ((double*)zr)[0];
    double r = ((double*)zr)[1];
    return exp(-bessel_j0n[n]*z)*gsl_sf_bessel_J0(bessel_j0n[n]*r);
}

double sumBessel(const double coeffs[MAX_BESSEL_TERMS], double z, double r) {
    if(r >= 1) return 0;
    double zr[2] = {z,r};
    return coeffSum(coeffs, &_sumBessel, zr, 1e-5, 1e-6);
}

double _sumBesselDerivR(unsigned int n, void* zr) {
    double z = ((double*)zr)[0];
    double r = ((double*)zr)[1];
    return -exp(-bessel_j0n[n]*z)*bessel_j0n[n]*gsl_sf_bessel_J1(bessel_j0n[n]*r);
}

double sumBesselDerivR(const double coeffs[MAX_BESSEL_TERMS], double z, double r) {
    if(r >= 1) return 0;
    double zr[2] = {z,r};
    return coeffSum(coeffs, &_sumBesselDerivR, zr, 1e-5, 1e-6);
}

double _sumBesselDerivZ(unsigned int n, void* zr) {
    double z = ((double*)zr)[0];
    double r = ((double*)zr)[1];
    return -bessel_j0n[n]*exp(-bessel_j0n[n]*z)*gsl_sf_bessel_J0(bessel_j0n[n]*r);
}

double sumBesselDerivZ(const double coeffs[MAX_BESSEL_TERMS], double z, double r) {
    if(r >= 1) return 0;
    double zr[2] = {z,r};
    return coeffSum(coeffs, &_sumBesselDerivZ, zr, 1e-5, 1e-6);
}
                                
////////////////////////////////////////////////////////////
                
void initDoubleBessel(struct DoubleBessel* BB, const double c0[MAX_BESSEL_TERMS], const double c1[MAX_BESSEL_TERMS]) {
    unsigned int n;
    for(n=1; n<MAX_BESSEL_TERMS; n++) {
        double ez = exp(-bessel_j0n[n]*BB->dz);
        double k = 1./(1-ez*ez);
        double c0n = (c0[n]-ez*c1[n])*k;
        double c1n = (-ez*c0[n]+c1[n])*k;
        BB->c0[n] = c0n;
        BB->c1[n] = c1n;
    }
}

double sumDoubleBessel(const struct DoubleBessel* BB, double z, double r) {
    return sumBessel(BB->c0, z, r) + sumBessel(BB->c1, BB->dz-z, r);
}

double sumDoubleBesselDR(const struct DoubleBessel* BB, double z, double r) {
    return sumBesselDerivR(BB->c0, z, r) + sumBesselDerivR(BB->c1, BB->dz-z, r);
}

double sumDoubleBesselDZ(const struct DoubleBessel* BB, double z, double r) {
    return sumBesselDerivZ(BB->c0, z, r) - sumBesselDerivZ(BB->c1, BB->dz-z, r);
}
