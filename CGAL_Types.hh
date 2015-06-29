#ifndef CGAL_TYPES_HH
#define CGAL_TYPES_HH


#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <vector>
using std::vector;
#include <list>
using std::list;
using std::pair;

/// Convenience class for ordered pairs
template<typename T>
class OrderedPair {
public:
    OrderedPair(const T& a, const T& b): first(a < b? a : b), second(a < b? b : a) { }
    bool operator<(const OrderedPair& rhs) const { return first < rhs.first || (first == rhs.first && second < rhs.second); }
    bool intersects(const OrderedPair& rhs) const { return first == rhs.first || first == rhs.second || second == rhs.first || second == rhs.second; }
    T first;
    T second;
};

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<double> K;

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

// p[olyhedra
typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Point_3 Point_3;
typedef Polyhedron::Vertex_handle Vertex_handle_3;

// implicit function isosurfaces
#include "GeomDef.hh"

#endif
