/// \file plot_slice.cc Load and render mesh slice data
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

/*
g++ -O3 --std=c++11 -o plot_slice -I../ -I../Analytical -I${MPMUTILS}/GeneralUtils/ -L${MPMUTILS}/GeneralUtils/ \
plot_slice.cc BrianFields.cc SVGAnalyticalRenderer.cc SVGPixelRenderer.cc \
SVGSliceRenderer.cc ../CellMatrix.cc ../Analytical/aCORN_EMirror_Field.c \
../Analytical/EndBesselCalc.c ../Analytical/WireplaneField.c -lgsl -lblas -lm -lMPMGeneralUtils


./plot_slice ../../elemesholve-bld/slices_xyz.dat
rsvg-convert -f pdf -o slices_xyz_1.pdf slices_xyz_1.svg
inkscape slices_xyz_1.svg --export-pdf=slices_xyz_1.pdf

inkscape BrianField.svg --export-pdf=BrianField.pdf
*/

#include "SVGSliceRenderer.hh"
#include "SVGPixelRenderer.hh"
#include "SVGAnalyticalRenderer.hh"
#include "aCORN_EMirror_Field.h"
#include "StringManip.hh"
#include "BrianFields.hh"
#include <stdio.h>

float grid_integral(const TriCubic& G, size_t i0[3], size_t axis, size_t x0 = 0, size_t x1 = -1) {
    float s = 0;
    if(x1==-1) x1 = G.NX[axis];
    for(size_t i = x0; i < x1; i++) {
        i0[axis] = i;
        s += G.at(i0);
    }
    s /= G.sx[axis];
    return s;
}

void plot_gridslice() {
    TriCubic G[3];
    load_Brian_mirror(G);
    
    const size_t PLOTX = 0;
    const size_t PLOTY = 2;
    const size_t PLOTZ = 1;
    
    SVGPixelRenderer P;
    float x0[3];
    float x1[3];
    size_t ii[3] = {0,0,0};
    G[0].gridpos(ii,x0);
    for(int i=0; i<3; i++) ii[i] = G[0].NX[i]-1;
    G[0].gridpos(ii,x1);
    P.make_grid(x0[PLOTX], x0[PLOTY], x1[PLOTX], x1[PLOTY], G[0].NX[PLOTX], G[0].NX[PLOTY] , true);
    
    float x[3] = { 0, 0, 0 };
    for(auto& p: P.pxls) {
        x[PLOTX] = p.pos(0.5, 0);
        x[PLOTY] = p.pos(0.5, 1);
        float E[3] = { G[0](x), G[1](x), G[2](x) };
        p.z = fabs(E[PLOTX]); // transverse field in V/cm
        //p.z = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2])*0.01; // field magnitude, V/cm
    }
    
    P.autoscale = false;
    P.zAxis.range.lo[0] = 0.005;
    P.zAxis.range.hi[0] = 500;
    P.zAxis.logscale = true;
    P.draw_axis = false;
    P.outcoord_scale = 0.1;
    P.write_svg("../../elemesholve-bld/BrianField.svg");
    
    //////////////////////////////////////////////////////
    // spot check axis field. Compare to analytical model.
    struct aCORN_EMirror M;
    init_aCORN_params(&M);
    init_aCORN_calcs(&M);
    x[0] = 3.0;
    x[1] = 0;
    for(x[2] = -6; x[2] <= 6; x[2] += 0.02) {
        double E[3];
        double xx[3] = { x[0], x[1], x[2] };
        calc_aCORN_field(&M, xx, E);
        double phi = calc_aCORN_potential(&M, xx);
        float GE[3];
        for(int i=0; i<3; i++) GE[i] = G[i].eval_linear(x);
        printf("%.2f\t%.3f\t%.3f\t%.3f\t\t%.3f\t%.3f\t%.3f\t\t%.3f\n", x[2], GE[0], GE[1], GE[2], E[0], E[1], E[2], phi);
    }
    
    //////////////
    // y integrals
    ii[PLOTZ] = G[0].NX[PLOTZ]/2;
    for(size_t ix=0; ix<G[0].NX[PLOTX]; ix++) {
        ii[PLOTX] = ix;
        G[0].gridpos(ii,x0);
        double s0 = grid_integral(G[PLOTX], ii, PLOTY, 0, 78);
        double s1 = grid_integral(G[PLOTX], ii, PLOTY, 80, -1);
        printf("%.1f\t%.3f\t%.3f\t%.3f\n", x0[PLOTX], s0, s1, s0+s1);
    }
    
}

int main(int argc, char** argv) {
    
    plot_gridslice();
    return EXIT_SUCCESS;
    
    if(argc != 2) {
        printf("./plot_slice <filename>.dat");
        return EXIT_FAILURE;
    }
    string infl = argv[1];
    ifstream is(infl.c_str(),  std::ios::in | std::ios::binary);
    
    for(int i=0; i<3; i++) {
        SVGSliceRenderer SR;
        SR.outcoord_scale = 0.1;
        SR.dcmode = SVGSliceRenderer::TRANSVERSE;
        SR.autoscale = true;
        double zmin=0, zmax=0;
        
        if(i==2) {
            SR.dcmode = SVGSliceRenderer::PHI;
            //zmin = 0;
            //zmax = 3.5;
        }
        
        if(SR.dcmode == SVGSliceRenderer::MEAN_PHI) {
            zmin = 0;
            zmax = 700;
        }
        
        if(SR.dcmode == SVGSliceRenderer::TRANSVERSE) {
            SR.zAxis.logscale = true;
            SR.autoscale = false;
            zmin = 0.005;
            zmax = 500;
        }
        
        SR.zAxis.range.expand(&zmin);
        SR.zAxis.range.expand(&zmax);
        
        SR.read(is);
        SR.SH.display();
        load_aCORN_simfield(SR);
        //load_Brian_simfield(SR);
        
        if(false && SR.dcmode == SVGSliceRenderer::PHI) {
            //SR.makeMeshVis(0.001);
            SR.vis_rmax2 = 5.5*5.5;
            SR.merge_gradient_regions(0.05);
        }
        
        string outfl = dropLast(infl,".")+"_"+to_str(i)+".svg";
        SR.write_svg(outfl);
    }
    
    is.close();
    
    return EXIT_SUCCESS;
}
