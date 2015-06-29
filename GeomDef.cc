// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "GeomDef.hh"
#include <cassert>
#include <iostream>

/////////////////////////////

void zcircle(Polylines& v, double x0, double y0, double z0, double r, int npts) {
    v.push_back(Polyline_3());
    for(int i = 0; i < npts; ++i) {
        double th = i*2*CGAL_PI/npts;
        v.back().push_back( K::Point_3(x0+r*cos(th), y0 + r*sin(th), z0) );
    }
    v.back().push_back(v.back().front());
}

/////////////////////////////

GeomWorldVolume::GeomWorldVolume(const CoordinateTransform* CT): GeomDomainFunction(CT) {
    zslices[-world_dz] = NULL;
    zslices[world_dz] = NULL;
}

const GeomSubVolume* GeomWorldVolume::findPoint(double x, double y, double z) const {
    // try searching around previously located volume
    while(prevLocated) { 
        const GeomSubVolume* V = prevLocated->findPoint(x,y,z);
        if(V) return (prevLocated = V);
        prevLocated = prevLocated->parent;
    }
    
    if(x*x + y*y > world_rr) return NULL;       // outside radius range
        
    // find z-slice subdomain
    auto it = zslices.lower_bound(z);
    if(it == zslices.end())  return NULL; // outside of volume z range
    if(it->second) prevLocated = it->second->findPoint(x,y,z);
    return prevLocated;
}    

int GeomWorldVolume::f(double x, double y, double z) const {
    if(x*x + y*y > world_rr || fabs(z) > world_dz) return 0;
    
    const GeomSubVolume* V = findPoint(x,y,z);
    if(V) return V->myLabel;
    
    return 1;
}

double GeomWorldVolume::meshsize(double x, double y, double z) const {
    const GeomSubVolume* V = findPoint(x,y,z);
    if(V) {
        double s = V->meshsize(x,y,z);
        return s < world_meshsize? s : world_meshsize;
    }
    return world_meshsize;
}

void GeomWorldVolume::add_features(Polylines& v) const {
    for(auto it = allSubVolumes.begin(); it != allSubVolumes.end(); it++) (*it)->add_features(v);
    for(int zi = -1; zi <= 1; zi += 2) zcircle(v, 0, 0, world_dz*zi, world_r, 1000);
}

void GeomWorldVolume::addLayer(GeomSubVolume* V) {
    assert(V);
    zslices[V->myBounds.zmin()] = NULL;
    zslices[V->myBounds.zmax()] = V;
}

void GeomWorldVolume::display() const {
    std::cout << "World volume z stack divisions:\n";
    for(auto it = zslices.begin(); it != zslices.end(); it++)
        std::cout << "\t" << it->first << "\t" << it->second << "\n";
}

//////////////////////////////

GeomSubVolume::GeomSubVolume(GeomWorldVolume* W, GeomSubVolume* P):
parent(P), theWorld(W) {
    W->allSubVolumes.push_back(this);
}

/////////////////////////////

GeomSphereFunction::GeomSphereFunction(GeomWorldVolume* W, double rsquared, double xc, double yc, double zc, GeomSubVolume* P):
GeomSubVolume(W,P), rr(rsquared), r(sqrt(rsquared)), x0(xc), y0(yc), z0(zc) {
    myBounds = CGAL::Bbox_3(x0-r,y0-r,z0-r, x0+r,y0+r,z0+r);
}

const GeomSubVolume* GeomSphereFunction::findPoint(double x, double y, double z) const {
    x -= x0;
    y -= y0;
    z -= z0;
    return (x*x + y*y + z*z < rr)? this : NULL;
}

void GeomSphereFunction::add_features(Polylines& v) const {
    int npts = 50;    
    int nj = 4;
    for(int j=0; j < nj; j++) {
        double zj = -1 + (j+0.5)*2./nj;
        double rxy = r*sqrt(1-zj*zj);
        zcircle(v, x0, y0, z0+ r*zj, rxy, npts);
    }
}

//////////////////////////////

GeomTorusFunction::GeomTorusFunction(GeomWorldVolume* W, double RR, double rr, GeomSubVolume* P):
GeomSubVolume(W,P), R(RR), r(rr) {
    myBounds = CGAL::Bbox_3(-(R+r),-(R+r),-r, (R+r),(R+r),r);
}

const GeomSubVolume* GeomTorusFunction::findPoint(double x, double y, double z) const {
    double x2=x*x, y2=y*y, z2=z*z;
    double w = x2 + y2 + z2 + R*R - r*r;
    return w*w - 4*R*R*(x2+y2) < 0? this : NULL;
}

void GeomTorusFunction::add_features(Polylines& v) const {    
    int nhelix = 1;
    for(int h = 0; h < nhelix; h++) {
        double phi0 = h*2*CGAL_PI/nhelix;
        
        int nloops = int(R/r)+1;
        int ppl = 25;
        int npts = ppl*nloops;
        v.push_back(Polyline_3());
        for(int i = 0; i < npts; ++i) {
            double w = R + r*cos(phi0 + (i%ppl)*2*CGAL_PI/ppl);
            double c0 = cos(i*2*CGAL_PI/npts);
            double s0 = sin(i*2*CGAL_PI/npts);
            K::Point_3 p(w*c0, w*s0, r*sin(phi0 + (i%ppl)*2*CGAL_PI/ppl));
            v.back().push_back(p);
        }
        v.back().push_back(v.back().front()); // close the line
        
    }
}

//////////////////////////////////

GeomXSliceBox::GeomXSliceBox(GeomWorldVolume* W, const CGAL::Bbox_3& B, GeomSubVolume* P):
GeomSubVolume(W,P) {
    myLabel = 1;
    myBounds = B;
    y0 = 0.5*(myBounds.ymax() + myBounds.ymin());
    z0 = 0.5*(myBounds.zmax() + myBounds.zmin());
}

void GeomXSliceBox::addContents(GeomSubVolume* V) {
    assert(V);
    V->parent = this;
    contents.push_back(V);
}

int GeomXSliceBox::boxnum(double x) const {
    return contents.size()*(x - myBounds.xmin())/(myBounds.xmax() - myBounds.xmin());
}

const GeomSubVolume* GeomXSliceBox::findPoint(double x, double y, double z) const {
    if(!(  myBounds.xmin() < x && x < myBounds.xmax() 
        && myBounds.ymin() < y && y < myBounds.ymax() 
        && myBounds.zmin() < z && z < myBounds.zmax())) return NULL;
    
    int n = boxnum(x);
    if(n < (int)contents.size()) {
        const GeomSubVolume* V = contents[n]->findPoint(x,y,z);
        if(V) return V;
    }
    
    return this;
}

double GeomXSliceBox::meshsize(double x, double y, double z) const {
    int n = boxnum(x);
    if(n < 0 || n >= (int)contents.size()) return 1000;
    double x0i = myBounds.xmin() + (n+0.5)*(myBounds.xmax() - myBounds.xmin())/contents.size();
    double rr = (x-x0i)*(x-x0i) + (z-z0)*(z-z0);
    
    double ymn = contents[n]->myBounds.ymin();
    double ymx = contents[n]->myBounds.ymax();
    if(y < ymn) rr += (y-ymn)*(y-ymn);
    if(y > ymx) rr += (y-ymx)*(y-ymx);
    
    double r = sqrt(rr);
    if(!(mesh_rmin < r)) r = mesh_rmin;
    return r;
}

void GeomXSliceBox::calcCenter(int n, int m, double& x, double& y, double& z) const {
    x = myBounds.xmin() + (n+0.5)*(myBounds.xmax() - myBounds.xmin())/m;
    y = y0;
    z = z0;
}

////////////////////////////////////////////

GeomYRod::GeomYRod(GeomWorldVolume* W, double xc, double yc, double zc, double ly, double rxy, GeomSubVolume* P):
GeomSubVolume(W,P), rr(rxy*rxy), r(rxy), x0(xc), y0(yc), z0(zc), dy(ly) {
    myBounds = CGAL::Bbox_3(x0-r, y0-0.5*dy, z0-r,  x0+r, y0+0.5*dy, z0+r);
}

const GeomSubVolume* GeomYRod::findPoint(double x, double y, double z) const {
    if(fabs(y-y0) >= dy/2) return NULL;
    x -= x0;
    z -= z0;
    return x*x + z*z < rr ? this : NULL;
}

double GeomYRod::meshsize(double, double, double) const { return 0.5*r; }

void GeomYRod::add_features(Polylines& v) const {
    int ppl = 25;
    int nloops = int(0.5/5*dy/r)+1;
    int npts = ppl*nloops;
    v.push_back(Polyline_3());
    for(int i=0; i<=npts; i++) {
        double l = double(i)/npts;
        double th = l*nloops*2*CGAL_PI;
        K::Point_3 p(x0 + r*cos(th), y0 + (l-0.5)*dy, z0 + r*sin(th));
        v.back().push_back(p);
    }
}

////////////////////////////////////////////

GeomRing::GeomRing(GeomWorldVolume* W, double rin, double rout, double zmn, double zmx, GeomSubVolume* P):
GeomSubVolume(W,P), ri(rin), ro(rout), rri(rin*rin), rro(rout*rout), z0(zmn<zmx? zmn : zmx), z1(zmx < zmn? zmn : zmx) {
    myBounds = CGAL::Bbox_3(-ro, -ro, z0,  ro, ro, z1);
}

const GeomSubVolume* GeomRing::findPoint(double x, double y, double z) const {
    if(!(z0 < z && z < z1)) return NULL;
    double rr = x*x + y*y;
    return rri < rr && rr < rro? this : NULL;
}

void GeomRing::add_features(Polylines& v) const {
    int npts = 111;
    zcircle(v, 0, 0, z0, ri, npts);
    zcircle(v, 0, 0, z0, ro, npts);
    zcircle(v, 0, 0, z1, ri, npts);
    zcircle(v, 0, 0, z1, ro, npts);
}

////////////////////////////////////////////

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
