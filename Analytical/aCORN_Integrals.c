/// \file aCORN_Integrals.c Integrated transverse fields in aCORN electrostatic mirror analytical model
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

// gcc -O3 -o aCORN_Integrals aCORN_Integrals.c aCORN_EMirror_Field.c EndBesselCalc.c WireplaneField.c -lgsl -lblas -lm

#include "aCORN_EMirror_Field.h"
#include <stdio.h>
#include <gsl/gsl_integration.h>

struct Eintegral_params {
    struct aCORN_EMirror* M;
    double x[3];
};

double Etransverse(double z, void* pp) {
    struct Eintegral_params* p = (struct Eintegral_params*)pp;
    double E[3];
    p->x[2] = z;
    calc_aCORN_field(p->M, p->x, E);
    return E[0];
}

int main(int argc, char** argv) {
    
    struct aCORN_EMirror M;
    init_aCORN_params(&M);
    M.wire_radius = 0.01;       // double wire radius
    init_aCORN_calcs(&M);
    M.wire_radius = 0;          // disable wire fields in subsequent calculations
    
    struct Eintegral_params EP;
    EP.M = &M;
    EP.x[0] = EP.x[1] = EP.x[2] = 0;
    
    const unsigned int INTEG_WS_SIZE = 512;
    gsl_integration_workspace* gslIntegrationWS = gsl_integration_workspace_alloc(INTEG_WS_SIZE);
    gsl_integration_cquad_workspace * gsl_cqd_ws = gsl_integration_cquad_workspace_alloc(INTEG_WS_SIZE/4);
    
    double rm,rp,e;
    size_t neval;
    double abs_err = 1e-3;
    double rel_err = 1e-3;
    
    gsl_function F;
    F.function = &Etransverse;
    F.params = &EP;
    
    for(EP.x[0]=0; EP.x[0]<=4.0001; EP.x[0] += 0.05) {
        rm = rp = 0;
        int err = gsl_integration_cquad(&F, -10, 0, abs_err, rel_err, gsl_cqd_ws, &rm, &e, &neval);
        err = gsl_integration_cquad(&F, 0, 10, abs_err, rel_err, gsl_cqd_ws, &rp, &e, &neval);
        printf("%.2f\t%.4f\t%.4f\t%.4f\n", EP.x[0], rm, rp, rm+rp);
    }
    
    gsl_integration_workspace_free(gslIntegrationWS);
    gsl_integration_cquad_workspace_free(gsl_cqd_ws);
    
    return EXIT_SUCCESS;
}
