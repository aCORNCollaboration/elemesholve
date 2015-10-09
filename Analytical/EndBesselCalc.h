/// \file EndBesselCalc.h Bessel function series for boundary conditions on endcap of cylinder
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef ENDBESSELCALC_H
#define ENDBESSELCALC_H

/// maximum number of terms to compute
#define MAX_BESSEL_TERMS 1024

/// initialize J_0 zeroes array
void init_j0n();

/// calculate coefficient for a flat circular patch of radius r, unit potential on unit-radius endcap
double flatcircle_coeff(unsigned int n, double r);

/// add circular patch of specified (normalized) radius, voltage to coefficients
void addCircle(double coeffs[MAX_BESSEL_TERMS], double r, double V);

/// sum Bessel terms at z,r (scaled to radius 1)
double sumBessel(const double coeffs[MAX_BESSEL_TERMS], double z, double r);

/// radial component dPhi/dr
double sumBesselDerivR(const double coeffs[MAX_BESSEL_TERMS], double z, double r);

/// axial component dPhi/dz
double sumBesselDerivZ(const double coeffs[MAX_BESSEL_TERMS], double z, double r);

/// coefficients for boundary conditions defined on both ends of a grounded cylinder
struct DoubleBessel {
    double c0[MAX_BESSEL_TERMS];        ///< bottom end coeffs
    double c1[MAX_BESSEL_TERMS];        ///< top end coeffs
    double dz;                          ///< distance between ends (normalized to radius units)
};

/// initialize two-ended Bessel from single-end coeffs. dz should already be set.
void initDoubleBessel(struct DoubleBessel* BB, const double c0[MAX_BESSEL_TERMS], const double c1[MAX_BESSEL_TERMS]);

/// sum two-ended bessel terms (scaled to radius 1)
double sumDoubleBessel(const struct DoubleBessel* BB, double z, double r);
/// two-ended bessel r derivative (scaled to radius 1)
double sumDoubleBesselDR(const struct DoubleBessel* BB, double z, double r);
/// two-ended bessel z derivative (scaled to radius 1)
double sumDoubleBesselDZ(const struct DoubleBessel* BB, double z, double r);

#endif
