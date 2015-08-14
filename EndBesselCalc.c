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
#include <stdio.h>

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

void addCircle(double r, double V, double coeffs[MAX_BESSEL_TERMS]) {
    unsigned int n;
    for(n=1; n<MAX_BESSEL_TERMS; n++) coeffs[n] += V*flatcircle_coeff(n,r);
}

double sumBessel(double z, double r, const double coeffs[MAX_BESSEL_TERMS]) {
    if(r >= 1) return 0;
    double s = 0;
    int ntiny = 0;
    unsigned int n;
    for(n=1; n<MAX_BESSEL_TERMS; n++) {
        double ds = coeffs[n]*exp(-bessel_j0n[n]*z)*gsl_sf_bessel_J0(bessel_j0n[n]*r);
        s += ds;
        if(fabs(ds) < fabs(1e-5*s) || fabs(ds) < 1e-6) ntiny++;
        else ntiny = 0;
        if(ntiny > 5) break;
    }
    //printf("Summed to %u terms.\n", n);
    return s;
}

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
    return sumBessel(z, r, BB->c0) + sumBessel(BB->dz-z, r, BB->c1);
}

/////////////////////////

// gcc -O3 -o EndBesselCalc EndBesselCalc.c -lgsl -lblas -lm

int main(int argc, char** argv) {
    init_j0n();
    
    double c0[MAX_BESSEL_TERMS] = {0};
    double c1[MAX_BESSEL_TERMS] = {0};
    
    addCircle(1, 1.0, c0);
    
    struct DoubleBessel BB;
    BB.dz = 1;
    initDoubleBessel(&BB, c0, c1);
    
    double z;
    double r = 0.3;
    for(z=0; z<0.1; z+=0.001) printf("z=%.3f\tr=%.2f:\t%.3f\n", z, r, sumDoubleBessel(&BB, z, r));
    
    return 0;
}
