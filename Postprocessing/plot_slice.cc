/// \file plot_slice.cc Load and render mesh slice data
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

// g++ -O3 --std=c++11 -o plot_slice -I../ -I${MPMUTILS}/GeneralUtils/ -L${MPMUTILS}/GeneralUtils/ plot_slice.cc SVGSliceRenderer.cc ../CellMatrix.cc -lMPMGeneralUtils
// ./plot_slice ../../elemesholve-bld/slices_xyz.dat
// rsvg-convert -f pdf -o slices_xyz_0.pdf slices_xyz_0.svg
// inkscape slices_xyz_2.svg --export-pdf=slices_xyz_2.pdf

#include "SVGSliceRenderer.hh"
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
        if(i<3) SR.dcmode = SVGSliceRenderer::PHI;
        //SR.zAxis.logscale = i<2;
        
        SR.read(is);
        if(SR.dcmode == SVGSliceRenderer::PHI) SR.merge_gradient_regions(0.01);
        
        string outfl = dropLast(infl,".")+"_"+to_str(i)+".svg";
        SR.write_svg(outfl);
    }
    
    is.close();
    
    return EXIT_SUCCESS;
}