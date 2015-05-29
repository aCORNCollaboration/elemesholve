#ifndef FEMESH_HH
#define FEMESH_HH

#include "UmfSparse.hh"

#include <stdio.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include "UmfSparse.hh"

#include <iostream>
using std::cout;
#include <vector>
using std::vector;
#include <list>
using std::list;
#include <map>
using std::map;
#include <cmath>
#include <cassert>

#include "Visr.hh"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Vector_2<K> Vec2;
typedef CGAL::Direction_2<K> Direction_2;

typedef CDT::Point Point;
typedef CDT::Vertex_handle Vertex_handle;

/// Calculations related to a triangle of the mesh
class MeshTriangle {
public:
    /// Constructor from mesh face
    MeshTriangle(const CDT::Face_handle& F);
    
    Point vertices[3];
    double Delta;
    double alpha[3];
    double beta[3];
    double gamma[3];
    double kij[3][3];
};

/// Boundary points on mesh
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

class FEMesh {
public:
    /// Constructor, from mesh
    FEMesh(CDT& M);
    /// Enumerate vertices given bound vertex list; generate K, ~K matrices
    void set_boundary_points(const vector<Vertex_handle>& bpts);
    /// Set boundary conditions
    void set_boundary_values(const vector<double>& v);
    /// solve
    void solve();
    /// Draw mesh
    void draw();
    
    double draw_logz = 0;
    
protected:
    CDT& cdt;                                   ///< the mesh
    size_t nbound;                              ///< number of bound degrees of freedom
    size_t nfree;                               ///< number of free degrees of freedom
    map<Vertex_handle, int> vertex_enum;        ///< vertices, enumerated for matrix index
    
    vector<MeshTriangle> meshtri;               ///< mesh triangle calculations, in finite_faces in domain iterator order
    UmfSparse K;
    double Knorm, Kcond;
    UmfSparse tK;
    vector<double> tKh;                         ///< RHS forcing terms of system
    vector<double> bvals;                       ///< boundary vertex values
    vector<double> cvals;                       ///< free vertex values
};

#endif
