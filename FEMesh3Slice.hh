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
    FEMesh3Slice(const C3t3& MM, K::Plane_3 PP, const CoordinateTransform* CT = NULL, size_t nd = 1): MeshSlice(MM, PP, CT, nd) { }
    
    /// calculate interpolated values at vertices; over faces
    void calc_vtxvals(const FEMesh3& F);
    
    /// get value at specified vertex
    double get_vtxval(MS_HDS::Vertex_const_handle h) const;
    /// draw projected with z = vtxvals
    void draw_projection() const;
    /// dump to file as binary HDS
    void dump_HDS(ostream& o, const FEMesh3& F) const;
        
protected:
    /// calculate vertex position in plane coordinates and z height
    K::Point_3 vertex_coordinate(MS_HDS::Vertex_const_handle h) const;

    map<MS_HDS::Vertex_const_handle, double> vtxvals;           ///< function values at plane intersection vertices
};

#endif
