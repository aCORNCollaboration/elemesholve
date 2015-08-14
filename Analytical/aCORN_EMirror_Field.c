/// \file aCORN_EMirror_Field.c
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "aCORN_EMirror_Field.h"
#include "WireplaneField.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

void init_aCORN(struct aCORN_EMirror* M) { 
    assert(M);
    init_j0n(); // do this here to be sure!
    
    M->E0 = 70;
    M->wire_radius = 0.005;
    M->wire_spacing = 0.2;
    M->mirror_radius = 5.5;
    M->entrance_radius = 3.26*2.54/2.;
    M->exit_radius = 3.928;
    M->plate_radius = 6.5;
    M->bore_radius = 10.;
    
    double V0 = wireplaneVOffset(M->wire_radius, M->wire_spacing) * M->E0/2;
    printf("Initializing aCORN geometry with mirror field %g V/cm, finite wire effect %.2f V.\n", M->E0, V0);
    
    unsigned int n;
    for(n=0; n<MAX_BESSEL_TERMS; n++) M->lowerField[n] = 0;
    addCircle(M->lowerField, M->entrance_radius / M->mirror_radius, V0);
    
    double c0[MAX_BESSEL_TERMS] = {0};
    double c1[MAX_BESSEL_TERMS] = {0};
    addCircle(c0, 1, 4*V0);
    addCircle(c0, M->plate_radius / M->bore_radius, -4*V0);
    addCircle(c0, M->exit_radius / M->bore_radius, V0);
    M->upperField.dz = 1;
    initDoubleBessel(&M->upperField, c0, c1);
}

void calc_aCORN_field(struct aCORN_EMirror* M, const double x[3], double E[3]) {
    E[0] = E[1] = E[2] = 0;
    
    // wireplane field contribution
    wireplaneField(M->wire_radius,
                   M->wire_spacing,
                   x[2], x[0] + 0.5*M->wire_spacing,    // note shift so x=0 falls between wires
                   &E[2], &E[0]);
    E[2] = (1-E[2])*M->E0/2;   // add bands linear field
    E[0] *= -M->E0/2;
    
    // matching to wall boundary conditions
    double r = sqrt(x[0]*x[0] + x[1]*x[1]); // radial coordinate
    if(x[2] < 0) {
        double znorm = -x[2]/M->mirror_radius;
        double rnorm = r/M->mirror_radius;
        if(znorm < 2e-3) znorm = 2e-3;
        if(rnorm > 1e-6) {
            double Er = -sumBesselDerivR(M->lowerField, znorm, rnorm)/M->mirror_radius;
            E[0] += Er*x[0]/r;
            E[1] += Er*x[1]/r;
        }
        E[2] -= -sumBesselDerivZ(M->lowerField, znorm, rnorm)/M->mirror_radius;
    } else if(x[2] < 0.999*M->bore_radius) {
        double znorm = x[2]/M->bore_radius;
        double rnorm = r/M->bore_radius;
        if(znorm < 2e-3) znorm = 2e-3;
        if(rnorm > 1e-6) {
            double Er = -sumDoubleBesselDR(&M->upperField, znorm, rnorm)/M->bore_radius;
            E[0] += Er*x[0]/r;
            E[1] += Er*x[1]/r;
        }
        E[2] += -sumDoubleBesselDZ(&M->upperField, znorm, rnorm)/M->mirror_radius;
    }
}


