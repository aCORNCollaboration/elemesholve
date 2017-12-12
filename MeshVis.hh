/// \file MeshVis.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef MESHVIS_HH
#define MESHVIS_HH

#include "CGAL_Types.hh"
#include "GeomDef.hh"
#include "FEMesh.hh"
#include "Visr.hh"

#include <set>
using std::set;

/// convert Point_3 to visualizer coordinate vector
inline vsr::vec3 PtoV(const Point_3 P) { return {{P.x(), P.y(), P.z()}}; }
/// visualize edge in a 3D triangulation
void draw_Tr3Edge(const Tr3Edge& L, const CoordinateTransform* T = NULL);
/// visualize polyhedron
void draw_polyhedron(const Polyhedron& P);
/// visualize facet of 3D triangulation
void draw_facet(const Facet& F);
/// visualize cell of 3D triangulation
void draw_cell(const Cell& C);
/// visualize polyline
void draw_polyline(const Polyline_3& p);

/// Visualization class for drawing 3D triangulated mesh
class C3t3_Vis {
public:
    /// Constructor
    C3t3_Vis(const C3t3& C, const CoordinateTransform* CT = NULL, bool di = false);
    /// visualize the mesh
    void draw(const MeshBoundaryConditions<Tr::Vertex_handle>* B = NULL) const;
    /// print info about the mesh
    void print() const;
    
    bool draw_interior;         ///< whether to draw interior edges
    const CoordinateTransform* T;       ///< coordinate transform between mesh and physical coordinates
    
    set< Tr3Edge > interior;    ///< interior edges
    set< Tr3Edge > surface;     ///< edges on surface
    set< Tr3Edge > features;    ///< edges belonging to 1D features
};

#endif
