/// \file SVGSliceRenderer.hh SVG renderings of mesh slice data
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef SVGSLICERENDERER_HH
#define SVGSLICERENDERER_HH

#include "HDS.hh"
#include "SliceHeader.hh"
#include <string>
using std::string;
#include <cfloat>

class SVGSliceRenderer: public HalfedgeFatDS< HDS_Vertex<3,float>, HDS_Face<4,float> > {
public:
    /// Constructor
    SVGSliceRenderer() { }
    
    /// face coloring modes
    enum DisplayColorMode {
        MAG_GRAD,       ///< |E|
        DOT_AXIAL,      ///< component of E along axis vector
        TRANSVERSE,     ///< magnitude of E transverse to axis vector
        PHI             ///< gradient-shaded potential
    } dcmode = MAG_GRAD;        ///< mode for color data
    bool logscale = false;      ///< z axis log scale
    
    SliceHeader<3,double> SH;           ///< info on slice plane
    double vis_rmax2 = DBL_MAX;         ///< radius^2 of SVG visualization view
    double vis_center[2] = {0,0};       ///< center of SVG visualization view
    double axis_direction[3] = {0,0,1}; ///< direction for axial/transverse components
    double outcoord_scale = 1.0;        ///< coordinate scaling for SVG output
    bool vis_all_inside = true;         ///< whether to require all points to be in vis_rmax2, or just some
    
    /// read in data file
    virtual void read(istream& is) { SH.read(is); HalfedgeFatDS< HDS_Vertex<3,float>, HDS_Face<4,float> >::read(is); }
    /// generate SVG output
    void write_svg(const string& fname) const;
};

#endif
