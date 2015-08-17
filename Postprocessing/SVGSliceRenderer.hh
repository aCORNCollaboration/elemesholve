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
#include "SVGBuilder.hh"
#include <string>
using std::string;
#include <cfloat>

/// D-dimensional plane equation
template<size_t D, typename val_tp>
class PlaneEquation {
public:
    /// Constructor
    PlaneEquation() { }
    /// Evaluate at point
    val_tp operator()(const val_tp* x) const {
        val_tp s = 0;
        for(size_t i=0; i<D; i++) s += P[i+1]*(x[i]-x0[i]);
        return s;
    }
    
    val_tp x0[D];       ///< relative centerpoint
    val_tp P[D+1];      ///< coefficients, y = P[0] + P[i+1]*x[i]
};

/// Color axis
class SVGGradientAxis {
public:
    /// Constructor
    SVGGradientAxis();
    /// normalize to axis internal coordinates
    double axisUnits(double x) const;
    /// derivative of axis transformation
    double dAxisUnits(double x) const;
    /// finalize range; set up text
    void finalize();
    /// Determine gradient mapping given face plane equation
    string gradient_remap(const PlaneEquation<2,float>& P) const;
    
    bool logscale = false;                                      ///< log scale setting
    BBox<1,double> range = empty_double_bbox<1>();              ///< axis range
    SVG::group* axisGroup = new SVG::group();                   ///< group containing axis information
    
    color::Gradient G;                                          ///< gradient color definition
    XMLBuilder* Gaxis = new XMLBuilder("linearGradient");       ///< axis plot SVG
    SVG::lingradient* base_gradient;                            ///< gradient in SVG form
};

class SVGSliceRenderer: public HalfedgeFatDS< HDS_Vertex<3,float>, HDS_Face<4,float> > {
public:
    /// Constructor
    SVGSliceRenderer() { }
    /// Destructor
    ~SVGSliceRenderer() { if(mesh_vis) mesh_vis->release(); }
    
    /// face coloring modes
    enum DisplayColorMode {
        MAG_GRAD,       ///< |E|
        DOT_AXIAL,      ///< component of E along axis vector
        TRANSVERSE,     ///< magnitude of E transverse to axis vector
        PHI             ///< gradient-shaded potential
    } dcmode = MAG_GRAD;        ///< mode for color data
    
    SliceHeader<3,double> SH;           ///< info on slice plane
    double vis_rmax2 = DBL_MAX;         ///< radius^2 of SVG visualization view
    double vis_center[2] = {0,0};       ///< center of SVG visualization view
    double axis_direction[3] = {0,0,1}; ///< direction for axial/transverse components
    double outcoord_scale = 1.0;        ///< coordinate scaling for SVG output
    bool vis_all_inside = true;         ///< whether to require all points to be in vis_rmax2, or just some 
    bool orientation = true;            ///< CCW or CW element orientation
    SVG::group* mesh_vis = NULL;        ///< mesh lines, prior to mesh gradient recombination
    SVGGradientAxis zAxis;              ///< color z axis
    
    /// read in data file
    virtual void read(istream& is) { SH.read(is); HalfedgeFatDS< HDS_Vertex<3,float>, HDS_Face<4,float> >::read(is); }
    
    /// generate plane equation for face
    PlaneEquation<2,float> getFacePlane(index_tp f) const;
    
    /// pre-scan vertices for z range
    void prescan_phi_range();
    
    /// generate mesh image
    void makeMeshVis(double wmin = 0);
    
    /// merge regions sharing common phi gradient
    virtual void merge_gradient_regions(double tol = 0.01);
    
    /// generate SVG output
    void write_svg(const string& fname);   
};

#endif
