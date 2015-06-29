#ifndef GEOMDEF_HH
#define GEOMDEF_HH

#include "CGAL_Types.hh"
#include <math.h>
#include <cassert>
#include <map>
using std::map;

typedef vector<K::Point_3> Polyline_3;
typedef list<Polyline_3> Polylines;

// utility function for z-oriented polyline
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
    YSqueezeTransform(double zc, double zw, double d): z0(zc), idz2(1./(zw*zw)), dip(d) { }
    
    /// meshing to physical coordinates; return physical dV / mesh dV
    virtual K::Point_3 mesh2phys_J(const K::Point_3& p, double& dV) const { dV = 1./fsqueeze(p.z()); return K::Point_3(p.x(), p.y()*dV, p.z()); }
    /// physical to meshing coordinates; 
    virtual K::Point_3 phys2mesh(const K::Point_3& p) const { return K::Point_3(p.x(), p.y()*fsqueeze(p.z()), p.z()); }
    
    /// squeeze scaling
    virtual double fsqueeze(double z) const { return 1 - dip/(1 + (z-z0)*(z-z0)*idz2); }
    
    double z0;
    double idz2;
    double dip;
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
    
    const CoordinateTransform* T;     ///< meshing coordinate transform
};

#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_cv.hpp>
/// Domain function wrapper
class GeomDomainFunctionWrapper {
public:
    GeomDomainFunctionWrapper(const GeomDomainFunction* GD): G(GD) { }
    
    /// call to labelling function for Point
    typedef int return_type;
    return_type operator()(const K::Point_3& p, const bool = true) const {
        if(!G->T) return G->f(p.x(), p.y(), p.z());
        K::Point_3 p2 = G->T->mesh2phys(p);
        return G->f(p2.x(), p2.y(), p2.z());
    }
    
    /// BOOST magic
    template <class T_>
    class Implicit_function_traits {
    public:
        typedef typename T_::Point Point;
    };
    
    template <class RT_, class Point_>
    class Implicit_function_traits<RT_ (*)(Point_)> {
    public:
        typedef typename boost::remove_reference< typename boost::remove_cv< Point_ >::type >::type Point;
    };
        
protected:
    const GeomDomainFunction* G;
};

class GeomSubVolume;

/// Geometry "world volume" (cylinder)
class GeomWorldVolume: public GeomDomainFunction {
public:
    GeomWorldVolume(const CoordinateTransform* CT = NULL);
    ~GeomWorldVolume() { }
    
    /// identify sub-volume containing point
    virtual const GeomSubVolume* findPoint(double x, double y, double z) const;
    
    /// labeling function
    virtual int f(double x, double y, double z) const;
    
    /// mesh sizing field, locating point in hierarchy
    virtual double meshsize(double x, double y, double z) const;    
    
    mutable const GeomSubVolume* prevLocated = NULL;    ///< volume responsible for previously-located point
    map<double, GeomSubVolume*> zslices;                ///< sub-volumes for slices in z /start: volume /start: volume / ... 
    vector<GeomSubVolume*> allSubVolumes;               ///< complete list of registered sub-volumes
    
    /// add internal item to Z stack
    void addLayer(GeomSubVolume* V);
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines& v) const;
    
    /// print contents summary
    void display() const;
    
    double world_dz = 1;
    double world_rr = 1;
    double world_r = 1;
};

/// Generic geometry sub-volume of world
#include <CGAL/Bbox_3.h>
class GeomSubVolume {
public:
    GeomSubVolume(GeomWorldVolume* W, GeomSubVolume* P = NULL);
    
    /// return subvolume containing point, or NULL if outside sub-volume
    virtual const GeomSubVolume* findPoint(double x, double y, double z) const = 0;
    
    /// mesh sizing field for point assumed in this volume
    virtual double meshsize(double x, double y, double z) const { return 0.1; }
    
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines&) const { }
    
    CGAL::Bbox_3 myBounds;              ///< bounding box for this item
    GeomSubVolume* parent;              ///< parent volume
    GeomWorldVolume* theWorld;          ///< world volume
    int myLabel = 0;                    ///< domain label
};

class GeomSphereFunction: public GeomSubVolume {
public:
    /// Constructor
    GeomSphereFunction(GeomWorldVolume* W, double rsquared, double xc, double yc, double zc, GeomSubVolume* P = NULL);
    /// return subvolume containing point, or NULL if outside sub-volume
    virtual const GeomSubVolume* findPoint(double x, double y, double z) const;
    
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines&) const;
    
    double rr;  ///< radius squared
    double r;   ///< radius
    double x0, y0, z0; ///< center position
};

class GeomTorusFunction: public GeomSubVolume {
public:
    /// Constructor
    GeomTorusFunction(GeomWorldVolume* W, double RR, double rr, GeomSubVolume* P = NULL);
    /// return subvolume containing point, or NULL if outside sub-volume
    virtual const GeomSubVolume* findPoint(double x, double y, double z) const;
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines&) const;
    
    double R;   ///< major radius
    double r;   ///< minor radius
};

/// Box sliced into equal divisions in x
class GeomXSliceBox: public GeomSubVolume {
public:
    /// Constructor
    GeomXSliceBox(GeomWorldVolume* W, const CGAL::Bbox_3& B, GeomSubVolume* P = NULL);
    /// return subvolume containing point, or NULL if outside sub-volume
    virtual const GeomSubVolume* findPoint(double x, double y, double z) const;
    /// mesh sizing field for point assumed in this volume
    virtual double meshsize(double x, double y, double z) const;
    
    /// add contents box
    void addContents(GeomSubVolume* V);
    /// get box number for given x
    int boxnum(double x) const;
    /// box center n for m subdivisions
    void calcCenter(int n, int m, double& x, double& y, double& z) const;
    
    double mesh_rmin = 0.1;             ///< minimum radius for adaptive mesh refinement
    
protected:
    vector<GeomSubVolume*> contents;    ///< contents for each x slice
    double x0, y0, z0;                  ///< box center
};

/// y-directed cylinder
class GeomYRod: public GeomSubVolume {
public:
    /// Constructor
    GeomYRod(GeomWorldVolume* W, double xc, double yc, double zc, double ly, double rxy, GeomSubVolume* P = NULL);
    /// return subvolume containing point, or NULL if outside sub-volume
    virtual const GeomSubVolume* findPoint(double x, double y, double z) const;
    /// mesh sizing field for point assumed in this volume
    virtual double meshsize(double x, double y, double z) const;
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines& v) const;
    
    double rr, r;               ///< radius^2, radius
    double x0, y0, z0;          ///< center coorindate
    double dy;                  ///< length in y 
};

/// ring electrode
class GeomRing: public GeomSubVolume {
public:    
    /// Constructor
    GeomRing(GeomWorldVolume* W, double rin, double rout, double zmn, double zmx, GeomSubVolume* P = NULL);
    /// return subvolume containing point, or NULL if outside sub-volume
    virtual const GeomSubVolume* findPoint(double x, double y, double z) const;
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines& v) const;
    
    double ri, ro, rri, rro;    ///< inner, outer radius and radius^2
    double z0, z1;              ///< z range
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// select Mesh_domain type
typedef CGAL::Labeled_mesh_domain_3<GeomDomainFunctionWrapper, K> Mesh_domain_base;
typedef CGAL::Mesh_domain_with_polyline_features_3<Mesh_domain_base> Mesh_domain;

#ifdef NOPE_CGAL_CONCURRENT_MESH_3
typedef CGAL::Mesh_triangulation_3< Mesh_domain,
                                    CGAL::Kernel_traits<Mesh_domain>::Kernel, // Same as sequential
                                    CGAL::Parallel_tag                        // Tag to activate parallelism
                                    >::type Tr;
#else
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_segment_index> C3t3;

typedef Tr::Cell Cell;
typedef Tr::Facet Facet;
typedef OrderedPair<Tr::Vertex_handle> Tr3Edge;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria   Facet_criteria;
typedef Mesh_criteria::Cell_criteria    Cell_criteria;
typedef Mesh_criteria::Edge_criteria    Edge_criterea;

/// Domain mesh size wrapper; conforms to concept MeshDomainField_3
class GeomDomainMeshsizeWrapper {
public:
    GeomDomainMeshsizeWrapper(const GeomDomainFunction* GD, double lscale = 1.0): G(GD), s(lscale) { }
    
    typedef K::FT FT;
    typedef Mesh_domain::Index Index;
    typedef K::Point_3 Point;
    typedef Edge_criterea::Point_3 Point_3;
    
    /// mesh sizing function
    FT operator()(const K::Point_3& p, const int, const Index&) const;
    
    /// edge sizing function
    FT operator()(const Point_3& p) const;
    
    double s;   ///< scaling factor from radius field
    
protected:
    const GeomDomainFunction* G;
};

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

#endif
