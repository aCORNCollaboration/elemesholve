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

inline vsr::vec3 PtoV(const Point_3& P) { return vsr::vec3(P.x(), P.y(), P.z()); }

void draw_Tr3Edge(const Tr3Edge& L, const CoordinateTransform* T = NULL);
void draw_polyhedron(const Polyhedron& P);
void draw_facet(const Facet& F);
void draw_cell(const Cell& C);
void draw_polyline(const Polyline_3& p);

class C3t3_Vis {
public:
    /// Constructor
    C3t3_Vis(const C3t3& C, const CoordinateTransform* CT = NULL, bool di = false);
    
    void draw(const MeshBoundaryConditions<Tr::Vertex_handle>* B = NULL) const;
    void print() const;
    
    bool draw_interior;
    const CoordinateTransform* T;
    
    set< Tr3Edge > interior;
    set< Tr3Edge > surface;
    set< Tr3Edge > features;
};

#endif
