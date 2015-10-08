/// \file SVGPixelRenderer.hh SVG renderings of "pixellated" grid data
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef SVGPIXELRENDERER_HH
#define SVGPIXELRENDERER_HH

#include "SVGGradientAxis.hh"

/// Struct for z-colored pixel data
class Pixel: public BBox<2,double> {
public:
    /// constructor
    Pixel() { }
    double z;           ///< z value
};

class SVGPixelRenderer {
public:
    /// Constructor
    SVGPixelRenderer() { }

    double vis_center[2] = {0,0};       ///< center of SVG visualization view
    double outcoord_scale = 1.0;        ///< coordinate scaling for SVG output
    bool autoscale = true;              ///< whether to auto-expand axis range
    SVGGradientAxis zAxis;              ///< color z axis
    
    vector<Pixel> pxls;                 ///< pixels in image
    
    /// generate pixel set to cover grid, optionally with pixel centers at grid locations
    void make_grid(double x0, double y0, double x1, double y1, size_t nx, size_t ny, bool centered = false);
    /// generate SVG output
    void write_svg(const string& fname);   
};

#endif
