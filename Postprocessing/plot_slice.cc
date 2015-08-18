/// \file plot_slice.cc Load and render mesh slice data
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

/*
g++ -O3 --std=c++11 -o plot_slice -I../ -I../Analytical -I${MPMUTILS}/GeneralUtils/ -L${MPMUTILS}/GeneralUtils/ \
plot_slice.cc SVGAnalyticalRenderer.cc SVGSliceRenderer.cc ../CellMatrix.cc ../Analytical/aCORN_EMirror_Field.c \
../Analytical/EndBesselCalc.c ../Analytical/WireplaneField.c -lgsl -lblas -lm -lMPMGeneralUtils


./plot_slice ../../elemesholve-bld/slices_xyz.dat
rsvg-convert -f pdf -o slices_xyz_1.pdf slices_xyz_1.svg
inkscape slices_xyz_1.svg --export-pdf=slices_xyz_1.pdf
*/

#include "SVGSliceRenderer.hh"
#include "SVGAnalyticalRenderer.hh"
#include "StringManip.hh"
#include <stdio.h>

int main(int argc, char** argv) {
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
            zmin = 0;
            zmax = 3.5;
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
