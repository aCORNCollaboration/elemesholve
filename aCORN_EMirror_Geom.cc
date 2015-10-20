/// \file aCORN_EMirror_Geom.cc classes for describing aCORN electrostatic mirror geometry
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "aCORN_EMirror_Geom.hh"
#include <cfloat>
#include <cmath>

// length units in [cm]

AEM_WireCap::AEM_WireCap():
wire_radius(0.005),     // nominal 0.005
wire_spacing(0.19),      // nominal 0.19
entrance_radius(3.26*2.54/2.),
exit_radius(3.928),
outer_radius(5.84),
thickness(1.08),
gridz(3*wire_radius),
platez(gridz-wire_radius),
wire_radius2(wire_radius*wire_radius),
entrance_radius2(entrance_radius*entrance_radius),
outer_radius2(outer_radius*outer_radius)
{ }

AEM_MirrorBands::AEM_MirrorBands(double zm, bool c):
continuous(c),
mirror_radius(5.45),
mirror_radius2(mirror_radius*mirror_radius),
band_period(0.7),        // nominal 0.7  cm
band_gap(0.03),          // nominal 0.03 cm
rel_gap(band_gap/band_period),
top_band_z(-band_gap),
V0_z(0),
zmin(zm),
E0(70)
{ }

double AEM_MirrorBands::band_V(double z) const {
    if(continuous) {
        if(z > top_band_z) return DBL_MAX;
        return E0*(V0_z - z);
    }
    
    double bz = band_coord(z);
    double nz = 0;
    bz = modf(bz, &nz);
    if(bz > 0.49*rel_gap && bz < 1-0.49*rel_gap) return nz*(E0*band_period);
    return DBL_MAX;
}
