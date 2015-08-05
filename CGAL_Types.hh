/// \file CGAL_Types.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

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
    /// Constructor
    OrderedPair(const T& a, const T& b): first(a < b? a : b), second(a < b? b : a) { }
    /// comparison operator
    bool operator<(const OrderedPair& rhs) const { return first < rhs.first || (first == rhs.first && second < rhs.second); }
    /// check if this ordered pair contains an item also in another
    bool intersects(const OrderedPair& rhs) const { return first == rhs.first || first == rhs.second || second == rhs.first || second == rhs.second; }
    T first;    ///< first item in ordered pair
    T second;   ///< second item in ordered pair
};

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
/// core geometry kernel type
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<double> K;

//typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;

/// vertex for 2D triangulation
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
/// face for 2D triangulation
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
/// data structure for 2D triangulation
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
/// 2D triangulation generator
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
/// mesh criteria for 2D triangulations
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
/// 2D vector
typedef CGAL::Vector_2<K> Vec2;
/// 2D direction
typedef CGAL::Direction_2<K> Direction_2;

/// point in 2D triangulation
typedef CDT::Point Point;
/// vertex in 2D triangulation
typedef CDT::Vertex_handle Vertex_handle;

/// 3D polyhedron type
typedef CGAL::Polyhedron_3<K> Polyhedron;
/// halfedge on polyhedron
typedef Polyhedron::Halfedge_handle Halfedge_handle;
/// point in 3D polyhedron
typedef Polyhedron::Point_3 Point_3;
/// vertex in 3D polyhedron
typedef Polyhedron::Vertex_handle Vertex_handle_3;

// implicit function isosurfaces
#include "GeomDef.hh"

#endif
