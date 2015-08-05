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
    /// Constructor
    ZSliceCylDomain(const CoordinateTransform* CT = NULL);
    
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
    
    double world_dz = 1;        ///< cylinder z half-length
    double world_rr = 1;        ///< cylinder radius^2
    double world_r = 1;         ///< cylinder radius
};

#include <CGAL/Bbox_3.h>
/// Generic geometry sub-volume of world
class GeomSubVolume {
public:
    /// Constructor
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

/// Spherical sub-volume
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
    double x0;  ///< center x coordinate
    double y0;  ///< center y coordinate
    double z0;  ///< center z coorindate
};

/// Toroidal sub-volume
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
    double x0;  ///< center x coordinate
    double y0;  ///< center y coordinate
    double z0;  ///< center z coorindate
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
    
    double rr;          ///< radius^2
    double r;           ///< radius
    double x0;          ///< center x coordinate
    double y0;          ///< center y coordinate
    double z0;          ///< center z coorindate
    double dy;          ///< length in y 
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
    
    double ri;  ///< inner radius
    double ro;  ///< outer radius
    double rri; ///< inner radius^2
    double rro; ///< outer radius^2
    double z0;  ///< start of z range
    double z1;   ///< end of z range
};

#endif