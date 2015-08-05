/// \file FEMesh2.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef FEMESH2_HH
#define FEMESH2_HH

#include "FEMesh.hh"
#include "CGAL_Types.hh"

/// Boundary points segment on 2D mesh
class MeshBoundary {
public:
    /// Constructor
    MeshBoundary(double z = 0, bool c = true): closed(c), z0(z) { }
    
    /// Determine boundary value at given vertex
    virtual double boundaryValue(const Vertex_handle& h) const { return z0; }
    
    /// Find added vertices after refinement
    void find_refined_boundary(const CDT& cdt);
    
    /// insert polygon of constraint points
    void insert_constraint_poly(const vector<Point>& v, CDT& cdt);
    /// make circular constraint
    void make_circle(Point ccenter, double r, CDT& cdt, int npts = 25);
    
    vector<Vertex_handle> vertices;     ///< polygonal boundary path
    bool closed;                        ///< whether the boundary is closed
    double z0;                          ///< default boundary value
};

/// Boundary conditions on 2D mesh
class MeshBoundaryConditions2: public MeshBoundaryConditions<CDT::Vertex_handle> {
public:
    /// Constructor
    MeshBoundaryConditions2() { }
    
    /// Calculate boundary values from boundary segments
    void calc_bvals(const CDT& cdt);
    
    vector<MeshBoundary> bounds;        ///< boundary segments
};

/// 2-dimensional solver class
class FEMesh2: public FEMeshSolver<2,CDT::Vertex_handle> {
public:
    /// Constructor, from 2D triangulation mesh
    FEMesh2(CDT& M);
    /// Destructor
    ~FEMesh2();
    
    /// Draw mesh
    void draw();
    double draw_logz = 0;       ///< optional log-z drawing scale
protected:
    /// vertex position dump subroutine
    virtual void dump_vertex_position(const CDT::Vertex_handle v, ostream& o) const;
    /// map from triangle to CellMatrix
    map<CDT::Face_handle, size_t> trcells;
};

/// demo test of 2D field solver
void mesh_demo_2D();
    
#endif
