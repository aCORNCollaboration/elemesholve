/// \file ZSliceCylDomain.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef ZSLICECYLDOMAIN_HH
#define ZSLICECYLDOMAIN_HH

#include "GeomDef.hh"

class GeomSubVolume;

/// Geometry "world volume" base class
class ZSliceCylDomain: public GeomDomainFunction {
public:
    ZSliceCylDomain(const CoordinateTransform* CT = NULL);
    ~ZSliceCylDomain() { }
    
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
    GeomSubVolume(ZSliceCylDomain* W, GeomSubVolume* P = NULL);
    
    /// return subvolume containing point, or NULL if outside sub-volume
    virtual const GeomSubVolume* findPoint(double x, double y, double z) const = 0;
    
    /// mesh sizing field for point assumed in this volume
    virtual double meshsize(double x, double y, double z) const { return 0.1; }
    
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines&) const { }
    
    CGAL::Bbox_3 myBounds;              ///< bounding box for this item
    GeomSubVolume* parent;              ///< parent volume
    ZSliceCylDomain* theWorld;          ///< world volume
    int myLabel = 0;                    ///< domain label
};

class GeomSphereFunction: public GeomSubVolume {
public:
    /// Constructor
    GeomSphereFunction(ZSliceCylDomain* W, double rsquared, double xc, double yc, double zc, GeomSubVolume* P = NULL);
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
    GeomTorusFunction(ZSliceCylDomain* W, double RR, double rr, GeomSubVolume* P = NULL);
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
    GeomXSliceBox(ZSliceCylDomain* W, const CGAL::Bbox_3& B, GeomSubVolume* P = NULL);
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
    GeomYRod(ZSliceCylDomain* W, double xc, double yc, double zc, double ly, double rxy, GeomSubVolume* P = NULL);
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
    GeomRing(ZSliceCylDomain* W, double rin, double rout, double zmn, double zmx, GeomSubVolume* P = NULL);
    /// return subvolume containing point, or NULL if outside sub-volume
    virtual const GeomSubVolume* findPoint(double x, double y, double z) const;
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines& v) const;
    
    double ri, ro, rri, rro;    ///< inner, outer radius and radius^2
    double z0, z1;              ///< z range
};

#endif