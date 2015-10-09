/// \file aCORN_Example.c Example use of aCORN analytical field calculators
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

// gcc -O3 -o aCORN_Example aCORN_Example.c aCORN_EMirror_Field.c EndBesselCalc.c WireplaneField.c -lgsl -lblas -lm

#include "aCORN_EMirror_Field.h"
#include <stdio.h>

int main(int argc, char** argv) {
    
    // create struct containing info for field calculations
    struct aCORN_EMirror M;
    // initialize struct with correct geometry dimensions and pre-calculated info
    init_aCORN_params(&M);
    init_aCORN_calcs(&M);
    
    // array for position coordinate (x,y,z) [cm]
    // wires are spaced in x direction (0.2 cm apart), with bore axis along z
    // wireplane at z=0; mirror region z<0. x=0 is in gap between wires at +/-0.1cm.
    double x[3] = {2.55, 2.55, 0};
    // array to hold E field values (in V/cm)
    double E[3];
    
    // scan over some points and get field
    for(x[2]=-5; x[2]<=5; x[2]+=0.1) {
        calc_aCORN_field(&M, x, E); // calculate field E(x)
        printf("%.2f\t%.2f\t%.2f\t\t%.3f\t%.3f\t%.3f\n", x[0], x[1], x[2], E[0], E[1], E[2]);
    }
    
    return 0;
}
