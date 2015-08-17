/// \file aCORN_EMirror_Field.h Analytical approximation for aCORN electrostatic mirror field
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef ACORN_EMIRROR_FIELD_H
#define ACORN_EMIRROR_FIELD_H

#include "EndBesselCalc.h"

/// Simplified aCORN electrostatic mirror geometry. Wires at z=0.
struct aCORN_EMirror {
    double E0;                  ///< asymptotic field [V/cm] inside mirror
    double V0;                  ///< finite-wire-size-induced potential offset [V]
    double wire_radius;         ///< grid wire radius [cm]
    double wire_spacing;        ///< spacing between wires [cm]
    double mirror_radius;       ///< radius of mirror bands
    double entrance_radius;     ///< wire grid entrance (proton direction) radius [cm]
    double exit_radius;         ///< wire grid exit radius [cm]
    double plate_radius;        ///< radius of grounded plate holding wire grid [cm]
    double bore_radius;         ///< radius of bore around assembly [cm]
    double lowerField[MAX_BESSEL_TERMS];        ///< Bessel expansion coefficients for field below grid
    struct DoubleBessel upperField;             ///< Bessel expansion coefficients for field above grid
};

/// initialize aCORN fields struct
void init_aCORN(struct aCORN_EMirror* M);

/// calculate electric field [V/cm] at specified (x,y,z) [cm] position. Wires spaced in x, centered at z=0; negative z inside mirror.
void calc_aCORN_field(struct aCORN_EMirror* M, const double x[3], double E[3]);

/// calculate electical potential [V] at specified (x,y,z) [cm] position.
double calc_aCORN_potential(struct aCORN_EMirror* M, const double x[3]);

#endif
