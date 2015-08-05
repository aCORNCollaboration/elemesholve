/// \file FEMesh3Slice.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef FEMESH3SLICE_HH
#define FEMESH3SLICE_HH

#include "MeshSlice.hh"
#include "FEMesh3.hh"
#include <string>
using std::string;
#include <cfloat>

/// Slice through a 3D finite element calculation mesh
class FEMesh3Slice: public MeshSlice {
public:
    /// Constructor
    FEMesh3Slice(const C3t3& MM, K::Plane_3 PP, const CoordinateTransform* CT = NULL): MeshSlice(MM, PP, CT) { }
    
    /// calculate interpolated values at vertices; over faces
    void calc_vtxvals(const FEMesh3& F);
    
    /// get value at specified vertex
    double get_vtxval(MS_HDS::Vertex_handle h) const;
    /// draw projected with z = vtxvals
    void draw_projection() const;
    
    /// SVG output of projection colored by |E|
    void write_svg(const string& fname, const FEMesh3& F) const;
    
    bool color_logz = true;             ///< whether to color by log scale
    double vis_rmax2 = DBL_MAX;         ///< radius^2 of SVG visualization view
    double vis_center[2] = {0,0};       ///< center of SVG visualization view
    bool vis_all_inside = true;         ///< whether to require all points to be in vis_rmax2, or just some
    
protected:
    map<MS_HDS::Vertex_handle, double> vtxvals;                 ///< function values at plane intersection vertices
    double grsq_min;                                            ///< min gradient^2
    double grsq_max;                                            ///< max gradient^2
};

#endif
