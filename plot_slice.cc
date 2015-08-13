/// \file plot_slice.cc Load and render mesh slice data
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015


// g++ -O3 --std=c++11 -o plot_slice -I${MPMUTILS}/GeneralUtils/ -L${MPMUTILS}/GeneralUtils/ plot_slice.cc SVGSliceRenderer.cc CellMatrix.cc -lMPMGeneralUtils

#include "SVGSliceRenderer.hh"
#include "StringManip.hh"
#include <stdio.h>

int main(int argc, char** argv) {
    if(argc != 2) {
        printf("./plot_slice <filename>.dat");
        return EXIT_FAILURE;
    }
    
    SVGSliceRenderer SR;
    SR.dcmode = SVGSliceRenderer::PHI;
    
    string infl = argv[1];
    ifstream is(infl.c_str(),  std::ios::in | std::ios::binary);
    SR.read(is);
    is.close();
    
    string outfl = dropLast(infl,".")+".svg";
    SR.write_svg(outfl);
        
    return EXIT_SUCCESS;
}
