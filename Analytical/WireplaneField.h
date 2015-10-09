/// \file WireplaneField.h Electrical field calculations for an infinite wireplane
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef WIREPLANEFIELD_H
#define WIREPLANEFIELD_H

/// field around wireplane producing unit field at infinity; output units are [distance] = [potential]/[field]
void wireplaneField(double r,           ///< wire radius
                    double d,           ///< wirespacing
                    double l,           ///< perpendicular distance from plane
                    double a,           ///< distance from wire center
                    double* Eperp,      ///< output for field component perpendicular to plane
                    double* Eparr       ///< output for field component parallel to plane
);

/// finite wire effect potential offset; output units are [distance] = [potential]/[field]
double wireplaneVOffset(double r,       ///< wire radius
                        double d        ///< wirespacing
);

/// potential around (grounded) wireplane producing unit field at infinity
double wireplanePotential(double r,     ///< wire radius
                          double d,     ///< wirespacing
                          double l,     ///< perpendicular distance from plane
                          double a      ///< distance from wire center
);

#endif
