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
class MeshBoundaryConditions2: public MeshBoundaryConditions {
public:
    /// Constructor
    MeshBoundaryConditions2() { }
    
    /// Calculate boundary values from boundary segments
    void calc_bvals(const CDT& cdt);
    
    vector<MeshBoundary> bounds;        ///< boundary segments
};

/// 2-dimensional solver class
class FEMesh2: public FEMeshSolver<2> {
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
    virtual void dump_vertex_position(const void* v, ostream& o) const;
};

/// demo test of 2D field solver
void mesh_demo_2D();
    
#endif
