/// \file GeomDef.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef GEOMDEF_HH
#define GEOMDEF_HH

#include "CGAL_Types.hh"
#include <math.h>
#include <cassert>
#include <map>
using std::map;

/// a polyline in 3D
typedef vector<K::Point_3> Polyline_3;
/// a list of 3D polylines
typedef list<Polyline_3> Polylines;

/// utility function for z-oriented polyline
void zcircle(Polylines& v, double x0, double y0, double z0, double r, int npts);

/// Base class for anisotropic meshing coordinate transforms
class CoordinateTransform {
public:
    /// Constructor
    CoordinateTransform() { }
    /// Destructor
    virtual ~CoordinateTransform() { }
    /// meshing to physical coordinates; return physical dV / mesh dV
    virtual K::Point_3 mesh2phys_J(const K::Point_3& p, double& dV) const = 0;
    /// meshing to physical coordinates; ignore Jacobean
    virtual K::Point_3 mesh2phys(const K::Point_3& p) const { double foo; return mesh2phys_J(p,foo); } 
    /// physical to meshing coordinates;
    virtual K::Point_3 phys2mesh(const K::Point_3& p) const = 0;
};

/// y coordinate "squeeze" transform
class YSqueezeTransform: public CoordinateTransform {
public:
    /// Constructor
    YSqueezeTransform(double zc, double zw, double d): z0(zc), idz2(1./(zw*zw)), dip(d) { }
    
    /// meshing to physical coordinates; return physical dV / mesh dV
    virtual K::Point_3 mesh2phys_J(const K::Point_3& p, double& dV) const { dV = 1./fsqueeze(p.z()); return K::Point_3(p.x(), p.y()*dV, p.z()); }
    /// physical to meshing coordinates; 
    virtual K::Point_3 phys2mesh(const K::Point_3& p) const { return K::Point_3(p.x(), p.y()*fsqueeze(p.z()), p.z()); }
    
    /// squeeze scaling
    virtual double fsqueeze(double z) const { return 1 - dip/(1 + (z-z0)*(z-z0)*idz2); }
    
    double z0;          ///< z at maximum squeeze
    double idz2;        ///< 1/dz^2 squeeze region length
    double dip;         ///< depth of squeeze
};

/// Base class for labelling geometry domains
class GeomDomainFunction {
public:
    /// Constructor
    GeomDomainFunction(const CoordinateTransform* CT = NULL): T(CT) { }
    
    /// labeling function
    virtual int f(double x, double y, double z) const = 0;
    /// mesh size function
    virtual double meshsize(double x, double y, double z) const = 0;
    /// linear feature preservation size function
    virtual double edgesize(double x, double y, double z) const { return meshsize(x,y,z); }
    /// recommendend "free space" mesh size
    double world_meshsize = 1.;
    
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines&) const { }
    
    CGAL::Bbox_3 myBounds;              ///< bounding box for domain
    
    const CoordinateTransform* T;       ///< meshing coordinate transform
};

#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_cv.hpp>
/// Domain function wrapper
class GeomDomainFunctionWrapper {
public:
    /// Constructor
    GeomDomainFunctionWrapper(const GeomDomainFunction* GD): G(GD) { }
    
    typedef int return_type;    ///< return type for labeling function
    /// call to labelling function for Point
    return_type operator()(const K::Point_3& p, const bool = true) const {
        if(!G->T) return G->f(p.x(), p.y(), p.z());
        K::Point_3 p2 = G->T->mesh2phys(p);
        return G->f(p2.x(), p2.y(), p2.z());
    }
    
    /// BOOST magic
    template <class T_>
    class Implicit_function_traits {
    public:
        typedef typename T_::Point Point;       ///< point type
    };
    /// more BOOST magic
    template <class RT_, class Point_>
    class Implicit_function_traits<RT_ (*)(Point_)> {
    public:
        typedef typename boost::remove_reference< typename boost::remove_cv< Point_ >::type >::type Point;      ///<
    };
        
protected:
    const GeomDomainFunction* G;        ///< domain function being wrapped
};

/// Base type for mesh domain
typedef CGAL::Labeled_mesh_domain_3<GeomDomainFunctionWrapper, K> Mesh_domain_base;
/// Mesh domain including specified features
typedef CGAL::Mesh_domain_with_polyline_features_3<Mesh_domain_base> Mesh_domain;

#ifdef NOPE_CGAL_CONCURRENT_MESH_3
typedef CGAL::Mesh_triangulation_3< Mesh_domain,
                                    CGAL::Kernel_traits<Mesh_domain>::Kernel, // Same as sequential
                                    CGAL::Parallel_tag                        // Tag to activate parallelism
                                    >::type Tr;
#else
/// 3D triangulation type                                    
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif

/// Triangulation mesh generation type
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_segment_index> C3t3;

/// 3D triangulation tetrahedral cell
typedef Tr::Cell Cell;
/// 3D triangulation triangular cell face
typedef Tr::Facet Facet;
/// Edge in 3D triangulation, specified by its vertices
typedef OrderedPair<Tr::Vertex_handle> Tr3Edge;

/// 3D meshing criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
/// 3D meshing facet criteria
typedef Mesh_criteria::Facet_criteria   Facet_criteria;
/// 3D meshing cell criteria
typedef Mesh_criteria::Cell_criteria    Cell_criteria;
/// 3D meshing edge criteria
typedef Mesh_criteria::Edge_criteria    Edge_criterea;

/// Domain mesh size wrapper; conforms to concept MeshDomainField_3
class GeomDomainMeshsizeWrapper {
public:
    /// Constructor
    GeomDomainMeshsizeWrapper(const GeomDomainFunction* GD, double lscale = 1.0): G(GD), s(lscale) { }
    
    typedef K::FT FT;                   ///< floating point type for sizing functions
    typedef Mesh_domain::Index Index;   ///< mesh domain index type
    typedef K::Point_3 Point;           ///< point type for coordinates
    typedef Edge_criterea::Point_3 Point_3;     ///< point type for edge coordinates
    
    /// mesh sizing function
    FT operator()(const K::Point_3& p, const int, const Index&) const;
    
    /// edge sizing function
    FT operator()(const Point_3& p) const;
    
    double s;   ///< scaling factor from radius field
    
protected:
    const GeomDomainFunction* G;        ///< domain providing sizing information
};

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

#endif
