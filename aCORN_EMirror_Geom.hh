/// \file aCORN_EMirror_Geom.hh classes for describing aCORN electrostatic mirror geometry
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

// length units in [cm]

/// wires cap geometry
class AEM_WireCap {
public:
    /// Constructor
    AEM_WireCap();
  
    const double wire_radius;           ///< grid wire radius
    const double wire_spacing;          ///< spacing between wires
    const double entrance_radius;       ///< entrance hole radius
    const double exit_radius;           ///< exit hole radius
    const double outer_radius;          ///< ground metal outer radius
    const double thickness;             ///< ground metal thickness
    const double gridz;                 ///< z center of wire grid
    const double platez;                ///< bottom of grounded plate

    const double wire_radius2;          ///< grid wire radius^2
    const double entrance_radius2;      ///< entrance hole radius^2
    const double outer_radius2;         ///< ground metal outer radius^2
};

/// electrostatic mirror bands geometry
class AEM_MirrorBands {
public:
    /// Constructor
    AEM_MirrorBands(double zm, bool c = false);

    /// band number coordinate; = x.5 at center of band x
    inline double band_coord(double z) const { return 0.5*rel_gap + (top_band_z - z)/band_period; }
    
    /// provide band potential [V] for point; DBL_MAX for non-band area
    double band_V(double z) const;
    
    bool continuous;                    ///< whether or not to discretize into bands
    const double mirror_radius;         ///< radius of mirror bands
    const double mirror_radius2;        ///< mirror radius squared
    const double band_period;           ///< band spacing band_period
    const double band_gap;              ///< gap between bands
    const double rel_gap;               ///< relative gap size band_gap/band_period
    const double top_band_z;            ///< z level of top of band
    const double zmin;                  ///< cutoff minimum z
    double V0_z;                        ///< z for nominal V=0
    const double E0;                    ///< average electric field, V/cm
};
