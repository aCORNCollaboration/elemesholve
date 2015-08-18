/// \file RoundRect.hh Rounded-edge rectangle geometry calculations
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef ROUNDRECT_HH
#define ROUNDRECT_HH

#include <vector>
using std::vector;
#include <utility>
#include <cmath>

/// Rounded-edge rectangle specification
class RoundRect {
public:
    /// Constructor
    RoundRect(double w, double h, double r0): xh(w/2), yh(h/2), r(r0), rr(r0*r0) { }
    
    /// check if point inside
    bool inside(double x, double y) const {
        x = fabs(x);
        y = fabs(y);
        if(x > xh || y > yh) return false;
        if(x < xh-r || y < yh-r) return true;
        x -= xh-r;
        y -= yh-r;
        return x*x + y*y < rr;
    }

    typedef std::pair<double, double> xypt;
    
    /// get points approximating perimeter
    vector<xypt> perimeterpoints(int ncurve = 20) const {
        vector<xypt> v;
        for(int i=0; i<ncurve; i++) {
            double th = M_PI/2*(ncurve-1-i)/(ncurve-1);
            v.push_back(xypt(xh+(cos(th)-1)*r, yh+(sin(th)-1)*r));
        }
        
        for(int i=ncurve-1; i>=0; i--) v.push_back(xypt(v[i].first, -v[i].second));
        for(int i=0; i<ncurve; i++) v.push_back(xypt(-v[i].first, -v[i].second));
        for(int i=ncurve-1; i>=0; i--) v.push_back(xypt(-v[i].first, v[i].second));
        
        return v;
    }
    
    double xh;  ///< half-width in x
    double yh;  ///< half-width in y
    double r;   ///< corner radius
    double rr;  ///< corner radius^2
};

#endif
