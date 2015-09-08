/// \file GeomSetup.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef GEOMSETUP_HH
#define GEOMSETUP_HH

#include "GeomDef.hh"
#include "aCORN_EMirror_Geom.hh"
#include "RoundRect.hh"
#include "FEMesh.hh"

/// Base class for 3D geometry setup
class GeomSetup: public MeshBoundaryConditions<Tr::Vertex_handle> {
public:
    /// Constructor
    GeomSetup() { }
    
    GeomDomainFunction* theWorld = NULL;        ///< "world" volume
    Polylines polylines;                        ///< feature-preserving curves
    
    /// Add feature preservation lines
    void add_features(Mesh_domain& D) const;
    
    /// Calculate boundary conditions for mesh
    virtual void calc_bvals(const C3t3& M) { }
};



////////////////////////////////////////////////////////////////////

/// Concentric spheres test geometry
class SphereTestDomain: public GeomDomainFunction {
public:
    /// Constructor
    SphereTestDomain(double r1, double r2, const CoordinateTransform* CT = NULL);
    /// labeling function
    virtual int f(double x, double y, double z) const;
    /// mesh size function
    virtual double meshsize(double x, double y, double z) const;
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines& v) const;
    
    double rr1; ///< outer sphere radius^2
    double rr2; ///< inner sphere radius^2
};

/// Concentric spheres test geometry
class SphereTestGeom: public GeomSetup {
public:    
    /// Constructor
    SphereTestGeom(const CoordinateTransform* CT = NULL);
    /// Calculate boundary conditions for mesh
    virtual void calc_bvals(const C3t3& M);
    
    double rinner;      ///< radius of inner sphere
};

////////////////////////////////////////////////////////////////////


/// wires cap geometry
class WireCap: public AEM_WireCap {
public:
    /// Constructor
    WireCap() { }
    
    /// geometry interior
    bool inVolume(double x, double y, double z, double rr) const;
    /// recommended meshing radius^2
    double mesh_radius2(double x, double y, double z, double rr) const;
    /// add mesh-guiding "features" to Polylines list
    void add_features(Polylines& v, double sqz = 1.0) const;
};

/// electrostatic mirror bands geometry
class MirrorBands: public AEM_MirrorBands {
public:
    /// Constructor
    MirrorBands(double zm, bool c): AEM_MirrorBands(zm, c) { }
    
    /// recommended meshing radius^2
    double mesh_radius2(double z, double rr) const;
    /// add mesh-guiding "features" to Polylines list
    void add_features(Polylines& v) const;
};

/// grounded side support bars
class SideBars {
public:
    /// Constructor
    SideBars();
    /// geometry interior
    bool inVolume(double x, double y) const;
    /// add mesh-guiding "features" to Polylines list
    void add_features(Polylines& v, double z0) const;
    /// recommended meshing radius^2
    double mesh_radius2(double x, double y) const;
    
    double x0;          ///< center x
    RoundRect RR;       ///< bar profile
};

/// aCORN electrostatic mirror
class EMirrorWorldVolume: public GeomDomainFunction {
public:
    /// Constructor
    EMirrorWorldVolume(const CoordinateTransform* CT = NULL);
    
    /// labeling function
    virtual int f(double x, double y, double z) const;
    /// linear feature preservation size function
    virtual double edgesize(double x, double y, double z) const;
    /// mesh sizing field
    virtual double meshsize(double x, double y, double z) const;
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines& v) const;
    
    WireCap WC;         ///< mirror top cap wire grid
    MirrorBands MB;     ///< mirror bands
    SideBars SB;        ///< grounded side bars
    double world_dz;    ///< world volume z half-length
    double world_r;     ///< world volume radius
    double world_rr;    ///< world volume radius^2
};

/// aCORN NG-6 electrostatic mirror test geometry
class EMirrorGeom: public GeomSetup {
public:
    /// Constructor
    EMirrorGeom(const CoordinateTransform* CT = NULL);
    /// Calculate boundary conditions for mesh
    virtual void calc_bvals(const C3t3& M);

    EMirrorWorldVolume myWorld; ///< world volume
};

#endif
