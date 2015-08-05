/// \file GeomDef.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "GeomDef.hh"
#include <cassert>
#include <iostream>

void zcircle(Polylines& v, double x0, double y0, double z0, double r, int npts) {
    v.push_back(Polyline_3());
    for(int i = 0; i < npts; ++i) {
        double th = i*2*CGAL_PI/npts;
        v.back().push_back( K::Point_3(x0+r*cos(th), y0 + r*sin(th), z0) );
    }
    v.back().push_back(v.back().front());
}

GeomDomainMeshsizeWrapper::FT GeomDomainMeshsizeWrapper::operator()(const K::Point_3& p, const int, const Index&) const { 
    if(!G->T) return G->meshsize(p.x(), p.y(), p.z())*s;
    double dV;
    K::Point_3 p2 = G->T->mesh2phys_J(p,dV);
    double ms = G->meshsize(p2.x(), p2.y(), p2.z())*s;
    double ws = G->world_meshsize*s/pow(dV,1./3.);
    return std::min(ms,ws);
}

GeomDomainMeshsizeWrapper::FT GeomDomainMeshsizeWrapper::operator()(const Point_3& p) const { 
    if(!G->T) return G->edgesize(p.x(), p.y(), p.z());
    K::Point_3 p2 = G->T->mesh2phys(p);
    return G->edgesize(p2.x(), p2.y(), p2.z());
}
