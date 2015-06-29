// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef GEOMSETUP_HH
#define GEOMSETUP_HH

#include "GeomDef.hh"
#include "FEMesh.hh"

/// Base class for 3D geometry setup
class GeomSetup: public MeshBoundaryConditions {
public:
    /// Constructor
    GeomSetup() { }
    
    GeomDomainFunction* theWorld = NULL;        /// "world" volume
    Polylines polylines;                        /// feature-preserving curves
    
    /// Add feature preservation lines
    void add_features(Mesh_domain& D) const;
    
    /// Calculate boundary conditions for mesh
    virtual void calc_bvals(const C3t3& M) { }
};



////////////////////////////////////////////////////////////////////

/// Concentric spheres test geometry
class SphereTestDomain: public GeomDomainFunction {
public:
    SphereTestDomain(double r1, double r2): rr1(r1*r1), rr2(r2*r2) { }
    virtual int f(double x, double y, double z) const;
    virtual double meshsize(double x, double y, double z) const { return 0.01; }
    double rr1, rr2;
};

/// Concentric spheres test geometry
class SphereTestGeom: public GeomSetup {
public:    
    /// Constructor
    SphereTestGeom();
    /// Calculate boundary conditions for mesh
    virtual void calc_bvals(const C3t3& M);
    
    //virtual FT operator()(const K::Point_3& p, const int, const Index&) const;
    
    double rinner;
};

////////////////////////////////////////////////////////////////////


/// wires cap geometry
class WireCap {
public:
    WireCap();
    
    /// geometry interior
    bool inVolume(double x, double y, double z, double rr) const;
    /// recommended meshing radius
    double mesh_radius(double x, double y, double z, double rr) const;
    /// add mesh-guiding "features" to Polylines list
    void add_features(Polylines& v, double sqz = 1.0) const;
    
    const double wire_radius;           ///< grid wire radius
    const double wire_radius2;          ///< grid wire radius
    const double wire_spacing;          ///< spacing between wires
    const double entrance_radius;       ///< entrance hole radius
    const double entrance_radius2;      ///< entrance hole radius^2
    const double exit_radius;           ///< exit hole radius
    const double outer_radius;          ///< ground metal outer radius
    const double outer_radius2;         ///< ground metal outer radius^2
    const double thickness;             ///< ground metal thickness
    const double gridz;                 ///< z center of wire grid
    const double platez;                ///< bottom of grounded plate
};

class MirrorBands {
public:
    MirrorBands(double zm, bool c = false);
    
    /// recommended meshing radius
    double mesh_radius(double z, double rr) const;
    /// add mesh-guiding "features" to Polylines list
    void add_features(Polylines& v) const;
    
    /// band number coordinate; = x.5 at center of band x
    inline double band_coord(double z) const { return 0.5*rel_gap + (top_band_z - z)/band_period; }
    
    /// provide band voltage for point; DBL_MAX for non-band area
    double band_V(double z) const;
    
    const bool continuous;      ///< whether or not to discretize into bands
    double mirror_radius;       ///< radius of mirror bands
    double mirror_radius2;      ///< radius^2 of mirror bands
    double band_period;         ///< band spacing band_period
    double band_gap;            ///< gap between bands
    double rel_gap;             ///< relative gap size band_gap/band_period
    double top_band_z;          ///< z level of top of band
    double zmin;                ///< cutoff minimum z
};

/// aCORN electrostatic mirror
class EMirrorWorldVolume: public GeomDomainFunction {
public:
    EMirrorWorldVolume(const CoordinateTransform* CT = NULL);
    
    /// labeling function
    virtual int f(double x, double y, double z) const;
    /// linear feature preservation size function
    virtual double edgesize(double x, double y, double z) const;
    /// mesh sizing field
    virtual double meshsize(double x, double y, double z) const;
    /// add mesh-guiding "features" to Polylines list
    virtual void add_features(Polylines& v) const;
    
    WireCap WC;
    MirrorBands MB;
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

    EMirrorWorldVolume myWorld;
};

#endif
