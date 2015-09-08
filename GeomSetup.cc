/// \file GeomSetup.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "GeomSetup.hh"
#include "Visr.hh"
#include "MeshVis.hh"
#include "CellMatrix.hh"


void GeomSetup::add_features(Mesh_domain& D) const {
    if(!theWorld || !theWorld->T) { D.add_features(polylines.begin(), polylines.end()); return; }
    
    Polylines pl2;
    for(auto it = polylines.begin(); it != polylines.end(); it++) {
        pl2.push_back(Polyline_3());
        for(auto it2 = it->begin(); it2 != it->end(); it2++)
            pl2.back().push_back(theWorld->T->phys2mesh(*it2));
    }
    D.add_features(pl2.begin(), pl2.end());
}

///////////////////////////////////////////////////////////

bool WireCap::inVolume(double x, double y, double z, double rr) const {
    if(z < platez || z > platez+thickness) return false;
    
    if(rr < entrance_radius2) {
        if(z > gridz + wire_radius) {
            if(z > platez + 0.5*thickness) {
                double r = sqrt(rr);
                double dz = z-(platez+0.75*thickness);
                double dr = r-(exit_radius + 0.25*thickness);
                if(dr > 0 || dr*dr + dz*dz < 0.25*0.25*thickness*thickness) return true;
            }
            return false;
        }
        double wx = x/wire_spacing-0.5+100;
        wx = (wx-int(wx)-0.5)*wire_spacing;
        return wx*wx + (z-gridz)*(z-gridz) < wire_radius2;
    }
    
    if(rr < outer_radius2) return true;
    
    double dr = sqrt(rr)-outer_radius;
    double dz = z-(platez+0.5*thickness);
    return dr*dr + dz*dz < 0.25*thickness*thickness;
}

double WireCap::mesh_radius2(double x, double y, double z, double rr) const {
    double mr2 = 0;
    if(rr < entrance_radius2) {
        double wx = x/wire_spacing-0.5+100;
        wx = (wx-int(wx)-0.5)*wire_spacing;
        mr2 = wx*wx + (z-gridz)*(z-gridz);
    } else {
        double dr = sqrt(rr)-entrance_radius;
        double dz = z-gridz;
        mr2 = dr*dr + dz*dz;
    }
    return std::max(mr2, wire_radius2);
}

void WireCap::add_features(Polylines& v, double sqz) const {
    //zcircle(v, 0, 0, platez, entrance_radius, 1000);
    zcircle(v, 0, 0, platez+thickness/2., entrance_radius, 1000);
    zcircle(v, 0, 0, platez+3*thickness/4., exit_radius, 1000);
    
    int ppc = 100;
    for(double x0 = -wire_spacing*int(entrance_radius/wire_spacing);
        x0 < entrance_radius; x0 += wire_spacing) {
        if(x0 <= -entrance_radius) continue;
        
        if(true) { // end circles; use with pi/6*r edge spacing -> 192s mesh 0.005 x 0.2; 730984 vertices
            for(int ymul = -1; ymul <= 1; ymul += 2) {
                Polyline_3 coil;
                for(int i=0; i<ppc; i++) {
                    double th = i*2*M_PI/ppc;
                    double x = x0 + wire_radius*cos(th);
                    double z = gridz + wire_radius*sin(th);
                    if(x < entrance_radius) {
                        double y = ymul*sqrt(entrance_radius2 - x*x);
                        coil.push_back(K::Point_3(x,y,z));
                    }
                }
                if(coil.size() >= 2) {
                    coil.push_back(coil.front());
                    v.push_back(coil);
                }
            }
        }
        
        if(false) { // helices; use with pi/6*r edge spacing -> 208s mesh 0.005 x 0.2; 783217 vertices
            //if(sqz != 1) printf("Protecting wires with squeeze %g\n", sqz);
            int ppl = 25;
            int nloops = int(0.5/sqz * 0.5 * entrance_radius/wire_radius)+1;
            int npts = ppl*nloops;
            
            Polyline_3 coil;
            for(int i=0; i<=npts; i++) {
                double l = double(i)/npts;
                double th = l*nloops*2*CGAL_PI;
                K::Point_3 p(x0 + wire_radius*cos(th), 2*(l-0.5)*entrance_radius, gridz + wire_radius*sin(th));
                if(p.x()*p.x() + p.y()*p.y() < entrance_radius2) coil.push_back(p);
            }
            if(coil.size() >= 2) v.push_back(coil);
        }
        
        if(false) { // edge lines
            int nedg = 1;
            for(int n=0; n<nedg; n++) {
                double th = n*2*M_PI/nedg;
                double x = x0 + wire_radius*cos(th);
                double z = gridz + wire_radius*sin(th);
                if(fabs(x0) >= entrance_radius) continue;
                double dy = sqrt(entrance_radius2 - x*x);
                Polyline_3 edge;
                for(int i=-1; i<=1; i+=2) {
                    double y = i*dy;
                    edge.push_back(K::Point_3(x, y, z));
                }
                v.push_back(edge);
            }
        }
        
       
        
        if(false) { // circles
            int ppl = 25;
            int nc = 2./sqz*entrance_radius/wire_radius;
            for(int c=0; c<=nc; c++) {
                Polyline_3 circ;
                double y = (c*2./nc - 1.)*entrance_radius;
                for(int i=0; i<ppl; i++) {
                    double th = i*2*CGAL_PI/ppl;
                    double x = x0 + wire_radius*cos(th);
                    if(x*x + y*y >= entrance_radius2) break;
                    double z = gridz + wire_radius*sin(th);
                    circ.push_back(K::Point_3(x,y,z));
                }
                if((int)circ.size() == ppl) {
                    circ.push_back(circ.front());
                    v.push_back(circ);
                }
            }
        }
    }
    
    return;
    
   
}

////////////////////////////////////////

double MirrorBands::mesh_radius2(double z, double rr) const {
    if(continuous) return 1e6;
    
    double znear = top_band_z;
    if(z < top_band_z) {
        double bz = band_coord(z);
        double ibz = int(bz);
        if(bz - ibz < 0.5) znear = top_band_z - ibz*band_period;
        else znear = top_band_z - (ibz+1)*band_period + band_gap;
    }
    
    double dr = sqrt(rr)-mirror_radius;
    double dz = z - znear;
    return std::max(0.49*band_gap*band_gap, dz*dz + dr*dr);
}

void MirrorBands::add_features(Polylines& v) const {
    if(continuous) {
        zcircle(v, 0, 0, top_band_z, mirror_radius, 1000);
        return;
    }
    for(double z = top_band_z; z > zmin; z -= band_period) {
        zcircle(v, 0, 0, z, mirror_radius, 1000);
        double zb = z - (band_period-band_gap);
        if(zb > zmin) zcircle(v, 0, 0, zb, mirror_radius, 1000);
    }
}

////////////////////////////////////////

SideBars::SideBars(): x0(3.25*2.54), RR(0.5*2.54, 2*2.54, 0.3) { }

bool SideBars::inVolume(double x, double y) const {
    return RR.inside(fabs(x)-x0, y);
}

void SideBars::add_features(Polylines& v, double z0) const {
    auto perim = RR.perimeterpoints();
    perim.push_back(perim[0]);
    
    Polyline_3 l;
    for(auto it = perim.begin(); it != perim.end(); it++) l.push_back(K::Point_3(x0+it->first, it->second, z0));
    v.push_back(l);
    for(auto it = l.begin(); it != l.end(); it++) *it = K::Point_3(-it->x(), it->y(), z0);
    v.push_back(l);
}

double SideBars::mesh_radius2(double x, double y) const {
    x = fabs(x) - x0 - (RR.xh-RR.r);
    y = fabs(y) - (RR.yh-RR.r);
    return std::max(x*x + y*y, RR.rr);
}

////////////////////////////////////////

EMirrorWorldVolume::EMirrorWorldVolume(const CoordinateTransform* CT):
GeomDomainFunction(CT),
MB(-10, true),
world_dz(10),
world_r(10),
world_rr(world_r*world_r) {
    world_meshsize = 1.5;
    myBounds = CGAL::Bbox_3(-world_r-1, -world_r-1, -world_dz-1,  world_r+1, world_r+1, world_dz+1);
}

int EMirrorWorldVolume::f(double x, double y, double z) const {
    if(fabs(z) > world_dz) return 0;
    double rr = x*x + y*y;
    if(rr > world_rr) return 0;
    if(SB.inVolume(x,y)) return 0;
    if(rr > MB.mirror_radius2 && z < MB.top_band_z) return 2;
    return WC.inVolume(x,y,z,rr)? 0 : 1;
}

double EMirrorWorldVolume::edgesize(double, double, double z) const {
    if(fabs(z-WC.gridz) < 1.01*WC.wire_radius) return (M_PI*WC.wire_radius)/6.;
    return 0.2;
}

void EMirrorWorldVolume::add_features(Polylines& v) const {
    for(int zi = -1; zi <= 1; zi += 2) zcircle(v, 0, 0, world_dz*zi, world_r, 1000);
    zcircle(v, 0, 0, MB.top_band_z, world_r, 1000);
    zcircle(v, 0, 0, -world_dz, MB.mirror_radius, 1000);
    zcircle(v, 0, 0, WC.platez, MB.mirror_radius, 1000);
    double sqz = 1.0;
    if(T) T->mesh2phys_J(K::Point_3(0,0,WC.gridz), sqz);
    WC.add_features(v, sqz);
    MB.add_features(v);
    for(int zi = -1; zi <= 1; zi += 2) SB.add_features(v, world_dz*zi);
    SB.add_features(v, MB.top_band_z);
}

double EMirrorWorldVolume::meshsize(double x, double y, double z) const {
    double rr = x*x + y*y;
    double rr1 = MB.mesh_radius2(z,rr);
    double rr2 = WC.mesh_radius2(x,y,z,rr);
    double rr3 = SB.mesh_radius2(x,y);
    return std::min(world_meshsize, sqrt(std::min(rr1,std::min(rr2,rr3))));
}

/////////////////////////////////////

EMirrorGeom::EMirrorGeom(const CoordinateTransform* CT): myWorld(CT)  {
    theWorld = &myWorld;
    myWorld.add_features(polylines);
}

void EMirrorGeom::calc_bvals(const C3t3& M) {
    bpts.clear();
    double rrsupports = pow(myWorld.SB.x0 - myWorld.SB.RR.xh,2);
    RoundRect RR2(2.02*myWorld.SB.RR.xh, 2.02*myWorld.SB.RR.yh, myWorld.SB.RR.r); // slightly larger rectangle for side bars
    
    for(auto it = M.facets_in_complex_begin(); it != M.facets_in_complex_end(); it++) {
        for(int i=0; i<4; i++) {
            if(i==it->second) continue;
            C3t3::Vertex_handle vi = it->first->vertex(i);
            if(bpts.count(vi)) continue;
            
            K::Point_3 p = vi->point();
            if(myWorld.T) p = myWorld.T->mesh2phys(p);
            
            double rr = p.x()*p.x() + p.y()*p.y();
            if( rr > 0.98*myWorld.world_rr                              // outer wall
                || RR2.inside(fabs(p.x())-myWorld.SB.x0, p.y())         // side bars
                || p.z() > myWorld.WC.gridz + 2*myWorld.WC.thickness) { // can top
                bpts[vi] = 0;
            } else {
                if(p.z() > myWorld.WC.gridz - 1.01*myWorld.WC.wire_radius) bpts[vi] = 0;        // grounded wires, cap
                else if(rr < 1.01*myWorld.MB.mirror_radius2) {
                    double v = myWorld.MB.band_V(p.z());
                    if(v != DBL_MAX) bpts[vi] = v;                                              // bands
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////

SphereTestDomain::SphereTestDomain(double r1, double r2, const CoordinateTransform* CT):
GeomDomainFunction(CT), rr1(r1*r1), rr2(r2*r2) {
    myBounds = CGAL::Bbox_3(-1.1*r1, -1.1*r1, -1.1*r1,  1.1*r1, 1.1*r1, 1.1*r1);
}

int SphereTestDomain::f(double x, double y, double z) const { 
    double rr = x*x + y*y + z*z; 
    return rr > rr1? 0 : rr > rr2? 1 : 0;
}

double SphereTestDomain::meshsize(double x, double y, double z) const {
    double rr = x*x + y*y + z*z;
    return (rr > rr2? sqrt(rr) : sqrt(rr2))*0.5;
}

void SphereTestDomain::add_features(Polylines& v) const {
    int nj = 4;
    double rinner = sqrt(rr2);
    for(int j=0; j < nj; j++) {
        double zj = -1 + (j+0.5)*2./nj;
        if(zj > 1) zj = 1;
        double rxy = rinner*sqrt(1-zj*zj);
        zcircle(v, 0, 0, rinner*zj, rxy, 100);
    }
}

///////////

SphereTestGeom::SphereTestGeom(const CoordinateTransform* CT): rinner(3.0) {
    double router = 7.0;
    theWorld = new SphereTestDomain(router, rinner, CT); 
    //theWorld->add_features(polylines);
    zcircle(polylines, 0, 0, 0, router, 100);
}
    
void SphereTestGeom::calc_bvals(const C3t3& M) {
    bpts.clear();
    for(auto it = M.facets_in_complex_begin(); it != M.facets_in_complex_end(); it++) {
        for(int i=0; i<4; i++) {
            if(i==it->second) continue;
            C3t3::Vertex_handle vi = it->first->vertex(i);
            K::Point_3 p = vi->point();
            if(theWorld->T) p = theWorld->T->mesh2phys(p);
            double rr = p.x()*p.x() + p.y()*p.y() + p.z()*p.z();
            bpts[vi] = (rr > rinner*rinner*1.1? 1. : 1./rinner);
        }
    }
}

//SphereTestGeom::FT SphereTestGeom::operator()(const K::Point_3& p, const int, const Index&) const {
//    FT sq_d_to_origin = CGAL::squared_distance(p, Point(CGAL::ORIGIN));
//   if(sq_d_to_origin < rinner*rinner) return 0.2*rinner; 
//    return 0.2*sqrt(sq_d_to_origin);
//}