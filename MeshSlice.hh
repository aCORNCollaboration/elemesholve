/// \file MeshSlice.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef MESHSLICE_HH
#define MESHSLICE_HH

#include "GeomDef.hh"
#include "MeshVis.hh"

#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_face_base.h>

#include <set>
using std::set;

/// Custom halfedge vertex type referring back to originating line segment
template <class Refs>
class MS_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_false, K::Point_3> {
public:
    /// Constructor
    MS_vertex(Tr3Edge s = Tr3Edge(NULL,NULL)): mySeg(s) { }
    Tr3Edge mySeg;      ///< edge this vertex is located on
    double c = 0;       ///< normalized coordinate between mySeg endpoints
};

/// Custom face type with reference to C3t3 cell
template <class Refs>
struct MS_face : public CGAL::HalfedgeDS_face_base<Refs> {
    /// Constructor
    MS_face() { }    
    Tr::Cell_handle myCell;     ///< cell in 3D mesh associated with this 2D slice face
};

/// Customized face, edge, vertex items
struct MS_items : public CGAL::HalfedgeDS_items_2 {
    template <class Refs, class Traits>
    struct Face_wrapper { typedef MS_face<Refs> Face; };
    template <class Refs, class Traits>
    struct Halfedge_wrapper { typedef CGAL::HalfedgeDS_halfedge_base<Refs> Halfedge; };
    template <class Refs, class Traits>
    struct Vertex_wrapper { typedef MS_vertex<Refs> Vertex; };
};

/// A traits struct for MeshSlice halfedges
struct MS_traits {
    typedef int Point_2; ///< arbitrary point type, not used here.
};

/// Half-edge data structure type using customized types
typedef CGAL::HalfedgeDS_default<MS_traits, MS_items> MS_HDS;
/// Identifier for edge in HDS
typedef OrderedPair<MS_HDS::Vertex_handle> MSEdge;
/// convert HDS vertex to visualizer coordinate vector
inline vsr::vec3 MSVtoV(MS_HDS::Vertex_handle h) { return PtoV(h->point()); }

/// Intersection of 3D mesh with a plane
class MeshSlice {
public:
    /// Constructor
    MeshSlice(const C3t3& MM, K::Plane_3 PP, const CoordinateTransform* CT = NULL, size_t nd = 1);
    
    /// visualize
    void draw() const;
    
    MS_HDS slice;       ///< Halfedge data structure on slice plane
    int verbose = 0;    ///< debugging verbosity level

protected:
    const C3t3& M;      ///< mesh being sliced
    
    K::Plane_3 P;               ///< slice plane
    K::Vector_3 pcoords[3];     ///< plane coordinate basis vectors
    
    /// wrapper for handling polymorphic type segment/plane intersections using Boost magic
    class Intersection_handler {
    public:
        /// Constructor
        Intersection_handler(K::Plane_3 PP, MS_HDS& MM, const CoordinateTransform* CT): e(NULL,NULL), P(PP), M(MM), T(CT) { }
        
        typedef MS_HDS::Vertex_handle result_type;      ///< return type for Boost
        
        /// handle point-like intersections
        MS_HDS::Vertex_handle operator()(const K::Point_3& p);
        /// handle in-plane degenerate intersections
        MS_HDS::Vertex_handle operator()(const K::Segment_3& s);
        
        /// test whether a cell edge from vertex k to l intersects plane; check/update internal intersections list
        MS_HDS::Vertex_handle edgeIntersectsPlane(Tr::Cell_handle C, int k, int l);
        
        Tr3Edge e;                                              ///< edge being checked
        const CoordinateTransform* T;                           ///< mesh coordinates transform
        K::Plane_3 P;                                           ///< intersecting plane
        MS_HDS& M;                                              ///< data structure reference
        int degeneracy;                                         ///< degeneracy of intersection (0 = midpoint; 1 = one vertex; 2 = both vertices)
        MS_HDS::Vertex_handle other_vtx;                        ///< other vertex for degenerate intersection with line
        map<Tr3Edge, MS_HDS::Vertex_handle> intersections;      ///< intersection points, listed by incident vertices
        map<C3t3::Vertex_handle, MS_HDS::Vertex_handle> degenerate_points;  ///< vertices hit by plane
        
    } myIntersectionHandler;    ///< mesh intersection handler
    
    set<Tr::Cell_handle> toProcess;                             ///< potentially intersecting cells yet to be processed
    set<Tr::Cell_handle> foundCells;                            ///< cells already evaluated in meshing process
    
    map<MSEdge, MS_HDS::Halfedge_handle> foundEdges;            ///< index to halfedge pairs by incident vertices.
    
    size_t findndomains;                                ///< number of disconnected domains to search for
    bool is_first_intersection = false;                 ///< flag for cells not connected to others
    C3t3::Cells_in_complex_iterator searchPosition;     ///< intersection search position in cells list
    /// finds one cell intersecting plane; sets up initial half-edge
    void find_first_intersection();
    /// find all intersections with cell
    void find_intersections(const Tr::Cell_handle& C, set< MS_HDS::Vertex_handle >& ixn_vertices);
    /// generate a new half-edge pair between vertices; returns halfedge incident into vto
    MS_HDS::Halfedge_handle new_edge(MS_HDS::Vertex_handle vfrom, MS_HDS::Vertex_handle vto);
    /// get already-created "opposite" half edge if available
    MS_HDS::Halfedge_handle get_opposite(MSEdge e) const;
    /// determine correct point ordering into vector for consistency with existing halfedges; provide starting halfedge from last to first.
    MS_HDS::Halfedge_handle order_face_points(set<MS_HDS::Vertex_handle>& vxs, vector<MS_HDS::Vertex_handle>& v) const;
    /// create half-edge path connecting intersection points; return list of traversed vertices.
    vector<Tr::Vertex_handle> make_face(const Tr::Cell_handle& C, set<MS_HDS::Vertex_handle>& vxs);
};

/// dot product between 3D point (as vector) and vector
inline double pdotv(K::Point_3 p, K::Vector_3 v) { return p.x()*v.x() + p.y()*v.y() + p.z()*v.z(); }

#endif
