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

void init_aCORN_params(struct aCORN_EMirror* M) {
    assert(M);
    
    M->E0 = 70;
    M->wire_radius = 0.005;
    M->wire_spacing = 0.2;
    M->wire_shift = 0.;
    M->mirror_radius = 5.45;
    M->mirror_length = 23;
    M->entrance_radius = 3.26*2.54/2.;
    M->exit_radius = M->entrance_radius; //3.928;
    M->plate_radius = 6.5;
    M->bore_radius = 10.;
}

void init_aCORN_calcs(struct aCORN_EMirror* M) { 
    assert(M);
    init_j0n(); // do this here to be sure!
    
    M->V0 = wireplaneVOffset(M->wire_radius, M->wire_spacing) * M->E0/2;
    printf("Initializing aCORN geometry with mirror field %g V/cm, finite wire effect %.2f V.\n", M->E0, M->V0);
    
    double l0[MAX_BESSEL_TERMS] = {0};
    double l1[MAX_BESSEL_TERMS] = {0};
    M->lowerField.dz = M->mirror_length/M->mirror_radius;
    addCircle(l0, M->entrance_radius / M->mirror_radius, M->V0);
    initDoubleBessel(&M->lowerField, l0, l1);
    
    double c0[MAX_BESSEL_TERMS] = {0};
    double c1[MAX_BESSEL_TERMS] = {0};
    addCircle(c0, 0.9, 3*M->V0);
    addCircle(c0, M->plate_radius / M->bore_radius, -3*M->V0);
    addCircle(c0, M->exit_radius / M->bore_radius, M->V0);
    M->upperField.dz = 1;
    initDoubleBessel(&M->upperField, c0, c1);
}

void calc_aCORN_field(struct aCORN_EMirror* M, const double x[3], double E[3]) {
    E[0] = E[1] = E[2] = 0;
    
    double r = sqrt(x[0]*x[0] + x[1]*x[1]); // radial coordinate
    // special case outside mirror
    if(x[2] < 0 && r > M->mirror_radius) {
        if(r > M->bore_radius) return;
        E[2] = M->E0 * (1-log(r/M->mirror_radius)/log(M->bore_radius/M->mirror_radius));
        double Er = -x[2] * M->E0 / (r * log(M->bore_radius/M->mirror_radius));
        E[0] += Er*x[0]/r;
        E[1] += Er*x[1]/r;
        return;
    }
    
    // wireplane field contribution
    if(M->wire_radius) {
        wireplaneField(M->wire_radius,
                       M->wire_spacing,
                       x[2], x[0] + M->wire_shift*M->wire_spacing,
                       &E[2], &E[0]);
        E[0] *= -M->E0/2;
        E[2] = (1-E[2])*M->E0/2;   // correct direction; add overall linear field
    } else if(x[2] < 0) {
        E[2] = M->E0;
    }
    
    // matching to wall boundary conditions
    if(x[2] < 0) {
        double znorm = -x[2]/M->mirror_radius;
        double rnorm = r/M->mirror_radius;
        if(znorm < 2e-3) znorm = 2e-3;
        if(rnorm > 1e-6) {
            double Er = -sumDoubleBesselDR(&M->lowerField, znorm, rnorm)/M->mirror_radius;
            E[0] += Er*x[0]/r;
            E[1] += Er*x[1]/r;
        }
        E[2] -= -sumDoubleBesselDZ(&M->lowerField, znorm, rnorm)/M->mirror_radius;
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

double calc_aCORN_potential(struct aCORN_EMirror* M, const double x[3]) {
    double phi = 0;
    double r = sqrt(x[0]*x[0] + x[1]*x[1]); // radial coordinate
    
    // special case outside mirror
    if(x[2] < 0 && r > M->mirror_radius) {
        if(r>M->bore_radius) return 0;
        return -x[2] * M->E0 * (1-log(r/M->mirror_radius)/log(M->bore_radius/M->mirror_radius));
    }
    
    
    // wireplane
    if(M->wire_radius) {
        phi = wireplanePotential(M->wire_radius,
                                 M->wire_spacing,
                                 x[2], x[0] + M->wire_shift*M->wire_spacing) * M->E0/2;
    } else {
        phi = fabs(x[2])*M->E0/2;
    }
    
    // main mirror field
    phi += -x[2]*M->E0/2;
        
    // subtract this off now... added back in by Bessel expansions
    phi -= M->V0;
    
    // matching to wall boundary conditions
    if(x[2] < 0) {
        double znorm = -x[2]/M->mirror_radius;
        double rnorm = r/M->mirror_radius;
        if(znorm < 2e-3) znorm = 2e-3;
        phi += sumDoubleBessel(&M->lowerField, znorm, rnorm);
    } else if(x[2] < 0.999*M->bore_radius) {
        double znorm = x[2]/M->bore_radius;
        double rnorm = r/M->bore_radius;
        if(znorm < 2e-3) znorm = 2e-3;
        phi += sumDoubleBessel(&M->upperField, znorm, rnorm);
    }
    
    return (phi == phi)? phi : 0;
}
