/// \file plot_slice.cc Load and render mesh slice data
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

/*
g++ -O3 --std=c++11 -o plot_slice -I../ -I../Analytical -I${MPMUTILS}/GeneralUtils/ -L${MPMUTILS}/GeneralUtils/ \
plot_slice.cc BrianFields.cc SVGAnalyticalRenderer.cc SVGGradientAxis.cc SVGPixelRenderer.cc \
SVGSliceRenderer.cc ../CellMatrix.cc ../Analytical/aCORN_EMirror_Field.c \
../Analytical/EndBesselCalc.c ../Analytical/WireplaneField.c -lgsl -lblas -lm -lMPMGeneralUtils


./plot_slice ../../elemesholve-bld/slices_xyz.dat
rsvg-convert -f pdf -o slices_xyz_1.pdf slices_xyz_1.svg
inkscape slices_xyz_1.svg --export-pdf=slices_xyz_1.pdf
*/

#include "SVGSliceRenderer.hh"
#include "SVGPixelRenderer.hh"
#include "SVGAnalyticalRenderer.hh"
#include "StringManip.hh"
#include "BrianFields.hh"
#include <stdio.h>

void plot_gridslice() {
    TriCubic G[3];
    load_Brian_mirror(G);
    
    const size_t PLOTX = 1;
    const size_t PLOTY = 2;
    
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
        p.z = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);
    }
    
    P.autoscale = false;
    P.zAxis.range.lo[0] = 0.1;
    P.zAxis.range.hi[0] = 10000;
    P.zAxis.logscale = true;
    P.write_svg("../../elemesholve-bld/BrianField.svg");
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
        //load_aCORN_simfield(SR);
        load_Brian_simfield(SR);
        
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
