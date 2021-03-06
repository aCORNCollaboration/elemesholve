/// \file SVGSliceRenderer.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "SVGSliceRenderer.hh"
#include "ProgressBar.hh"
#include "CellMatrix.hh"

/// extended info for face polygons to allow post-processed coloring
class SVGFacePoly: public SVG::polygon {
public:
    /// Constructor
    SVGFacePoly(const string& st = ""): SVG::polygon(st) { }
    
    /// expand dimensions to cover up gaps
    void expand(double mul = -0.01) {
        vector< pair<double,double> > offsets;
        double dx,dy;
        for(size_t i=0; i<pts.size(); i++) {
            size_t iprev = (i==0? pts.size()-1 : i-1);
            size_t inext = (i==pts.size()-1?  0 : i+1);
            dx = pts[inext].first - pts[iprev].first;
            dy = pts[inext].second - pts[iprev].second;
            offsets.push_back(pair<double, double>(dy*mul, -dx*mul));
        }
        for(size_t i=0; i<pts.size(); i++) {
            pts[i].first += offsets[i].first;
            pts[i].second += offsets[i].second;
        }
        
    }
    
    double orientation_sum() const {
        double s = 0; 
        for(size_t i=0; i<pts.size(); i++) {
            size_t iprev = (i==0? pts.size()-1 : i-1);
            s += (pts[i].first - pts[iprev].first)*(pts[i].second + pts[iprev].second);
        }
        return s;
    }
    
    vector<double> vtxz;        ///< vertex z info for coloring
    index_tp face;              ///< face index in HDS 
};

/////////////////////////////////

void SVGSliceRenderer::prescan_phi_range() {
    for(auto it = vertices.begin(); it != vertices.end(); it++) {
        if(pow(it->x[0]-vis_center[0],2) + pow(it->x[1]-vis_center[1],2) > vis_rmax2) continue;
        double z = it->x[2];
        zAxis.range.expand(&z);
    }
    cout << "Z range: " << zAxis.range.lo[0] << " -- " << zAxis.range.hi[0] << "\n";
}

void SVGSliceRenderer::makeMeshVis(double wmin) {
    if(mesh_vis) mesh_vis->release();
    mesh_vis = new SVG::group();
    mesh_vis->retain();
    mesh_vis->attrs["style"]="stroke:black;stroke-opacity:0.25;stroke-linecap:round"; 
    cout << "Drawing " << edges.size()/2 << " edges...\n";
    ProgressBar* PB = new ProgressBar(edges.size());
    int nedges = 0; 
    for(index_tp e = 0; e < edges.size(); e++) {
        PB->update(nedges++);
        index_tp eopp = edges[e].opposite;
        if(eopp <= e) continue;
        double w = 0.05*pow(faces[edges[e].face].x[0], 1./3.);
        if(w<wmin) continue;
        HDS_Vertex<3,float>& v1 = vertices[edges[e].vtx];
        HDS_Vertex<3,float>& v2 = vertices[edges[eopp].vtx];
        SVG::line* l = new SVG::line(outcoord_scale*v1.x[0],
                                     outcoord_scale*v1.x[1], 
                                     outcoord_scale*v2.x[0], 
                                     outcoord_scale*v2.x[1]);
        l->attrs["stroke-width"] = to_str(outcoord_scale*w);
        mesh_vis->addChild(l);
    }
    delete PB;
}

void SVGSliceRenderer::write_svg(const string& fname) {
    // output file
    std::ofstream o;
    o.open (fname);
    SVG::svg::make_standalone_header(o);
    
    // top-level SVG element
    XMLBuilder::indent = "\t";
    SVG::svg s;
    s.addChild(new SVG::title("elemesholve slice"));
    //s.attrs["shape-rendering"] = "crispEdges"; // "geometricPrecision";
    
    // style definitions group
    SVG::defs* d = new SVG::defs();
    s.addChild(d);
        
    d->addChild(zAxis.base_gradient);
    s.addChild(zAxis.axisGroup);
    
    // face polygons group
    SVG::group* g = new SVG::group();
    s.addChild(g);
    
    if(mesh_vis) s.addChild(mesh_vis);
    
    BBox<2,double> BB = empty_double_bbox<2>(); // image range bounding box
    
    double vtxpt[3];
    vector<SVGFacePoly*> facepoly; // face polygons for color assignment
    int nfaces = 0;
    int noriented[2] = {0,0};
    cout << "Drawing " << faces.size() << " faces...\n";
    ProgressBar* PB = new ProgressBar(2*faces.size());
    for(auto it = faces.begin(); it != faces.end(); it++) {
        PB->update(nfaces++);
        if(!it->edge) continue;
        
        // the polygon
        //double w =  0.05*pow(C.area(), 1./3.);      // stroke width scales with (cell volume)^1/3 ~ cell width
        SVGFacePoly* p = new SVGFacePoly(); //stroke_borders? "stroke-width:"+to_str(w) : "");
        p->face = nfaces-1;
        
        // collect polygon points
        index_tp start_edge = it->edge;
        index_tp current_edge = start_edge;
        size_t nOutside = 0;
        do {
            current_edge = edges[current_edge].next;
            assert(current_edge && current_edge < edges.size());
            for(int i=0; i<3; i++) vtxpt[i] = vertices[edges[current_edge].vtx].x[i];
            for(int i=0; i<2; i++) vtxpt[i] -= vis_center[i];
            if(vtxpt[0]*vtxpt[0] + vtxpt[1]*vtxpt[1] > vis_rmax2) nOutside++;
            p->addpt(vtxpt[0]*outcoord_scale, vtxpt[1]*outcoord_scale);
            p->vtxz.push_back(vtxpt[2]);
        } while(current_edge != start_edge);
        
        // skip invalid polygons outside range
        if(p->pts.size() < 3 || !((vis_all_inside && !nOutside) || (!vis_all_inside && nOutside < p->pts.size()))) {
            delete p;
            continue;
        }
        
        noriented[p->orientation_sum() > 0]++;
        
        // collect appropriate data         
        if(dcmode != PHI) {
            double z = 0;
            if(dcmode == MEAN_PHI) {
                for(auto it = p->vtxz.begin(); it != p->vtxz.end(); it++) z += *it;
                z /= p->vtxz.size();
            }
            p->vtxz.clear();
            if(dcmode == MAG_GRAD || dcmode == TRANSVERSE) {
                for(int i=0; i<3; i++) z += it->x[i+1]*it->x[i+1];
                if(dcmode == MAG_GRAD) z = sqrt(z);
            }
            if(dcmode == DOT_AXIAL || dcmode == TRANSVERSE) {
                double z2 = 0;
                for(int i=0; i<3; i++) z2 += it->x[i+1]*axis_direction[i];
                if(dcmode == TRANSVERSE) {
                    z = sqrt(z - z2*z2);
                } else z = z2;
            }

            if(autoscale) zAxis.range.expand(&z);
            p->vtxz.push_back(z);
        }
        if(autoscale) 
            for(auto it = p->vtxz.begin(); it != p->vtxz.end(); it++) 
                zAxis.range.expand(&*it);
        
        // expand image bounding box and insert polygon
        for(auto it = p->pts.begin(); it != p->pts.end(); it++) {
            vtxpt[0] = it->first;
            vtxpt[1] = it->second;
            BB.expand(vtxpt);
        }
        g->addChild(p);
        facepoly.push_back(p);
    }
    //cout << "Orientation: " << noriented[0] << " and " << noriented[1] << "\n";
    orientation = noriented[0] > noriented[1];
    
    // assign polygon colors by z
    zAxis.finalize();
    int nsubgrad = 0;
    for(auto it = facepoly.begin(); it != facepoly.end(); it++) {
        PB->update(nfaces++);
        SVGFacePoly* p = *it;
        if(!p->vtxz.size()) continue;
        p->expand(orientation? 0.02 : -0.02);
        
        // apply flat or gradient fill
        if(p->vtxz.size() == 1) {
            p->vtxz[0] = zAxis.axisUnits(p->vtxz[0]);
            p->vtxz[0] = p->vtxz[0]<0? 0 : p->vtxz[0]>1? 1 : p->vtxz[0];
            p->attrs["fill"] = "#"+color::rgb(zAxis.G.hsvcolor(p->vtxz[0])).asHexString();
        } else if(p->vtxz.size() >= 3) {
            XMLBuilder* Gi = new XMLBuilder("linearGradient");
            Gi->attrs["id"] = zAxis.base_gradient->attrs["id"]+"_"+to_str(nsubgrad++);
            PlaneEquation<2,float> FP = getFacePlane(p->face);
            for(int i=0; i<2; i++) { FP.x0[i] *= outcoord_scale; FP.P[i+1] /= outcoord_scale; }
            Gi->attrs["gradientTransform"] = zAxis.gradient_remap(FP);
            Gi->attrs["xlink:href"] = "#"+zAxis.base_gradient->attrs["id"];
            d->addChild(Gi);
            p->attrs["fill"] = "url(#" + Gi->attrs["id"] + ")";
        }
    }
    
    delete PB;
    
    // scale/move axis into position
    double yscale = BB.dl(1);
    zAxis.axisGroup->attrs["transform"] = "translate(" + to_str(BB.pos(1.1,0)) + " " + to_str(BB.pos(0.5,1) - 0.5*yscale) + ") scale(" + to_str(yscale) + ")";
    // expand to final display window
    vtxpt[0] = BB.pos(1.5,0);
    vtxpt[1] = BB.hi[1];
    BB.expand(vtxpt);
    
    s.setView(BB, 10);
    s.write(o);
    o.close();
}

struct gradient_merge_params {
    SVGSliceRenderer* S;        ///< data
    PlaneEquation<2,float> pl;  ///< reference plane equation
    index_tp f0;                ///< starting face
    double tol;                 ///< absolute error tolerance for plane equation agreement
    set<index_tp> skip_faces;   ///< previously-checked non-merge faces
};

bool gradient_merge(index_tp e, void* params) {
    gradient_merge_params& P = *(gradient_merge_params*)params;
    SVGSliceRenderer& S = *P.S;
    index_tp eopp = S.edges[e].opposite;
    index_tp f1 = S.edges[eopp].face;
    if(P.f0 >= f1 || P.skip_faces.count(f1)) return false;
    P.skip_faces.insert(f1);
    
    index_tp e1 = eopp;
    do {
        const HDS_Vertex<3,float>& v = S.vertices[S.edges[e1].vtx];
        if(fabs(P.pl(v.x) - v.x[2]) > P.tol) return false;
        e1 = S.edges[e1].next;
    } while(e1 != eopp);
    
    return true;
}

void SVGSliceRenderer::merge_gradient_regions(double tol) {
    
    if(autoscale) prescan_phi_range();
    gradient_merge_params gmp;
    gmp.S = this;
    gmp.tol = zAxis.range.dl(0) * tol;
    
    for(index_tp f = 1; f < faces.size(); f++) {
        gmp.skip_faces.clear();
        gmp.skip_faces.insert(0);
        gmp.f0 = f;
        gmp.skip_faces.insert(f);
        gmp.pl = getFacePlane(f);
        expand_merge(f, &gradient_merge, &gmp);
    }
}

float dot3(const float* x1, const double* x2) { return x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; }

PlaneEquation<2,float> SVGSliceRenderer::getFacePlane(index_tp f) const {
    PlaneEquation<2,float> P;
    const HDS_Vertex<3,float>& v = vertices[edges[faces[f].edge].vtx]; // a vertex on the plane
    P.P[0] = v.x[2];
    for(int i=0; i<2; i++) {
        P.x0[i] = v.x[i];
        P.P[i+1] = dot3(&faces[f].x[1], SH.basis[i]);
    }
    return P;
}

/*
 * bool drawEdges = !drawFaces;
 * if(drawEdges) {
 *    cout << "Drawing " << foundEdges.size() << " edges...\n";
 *    SVG::group* g = new SVG::group();
 *    //g->attrs["opacity"]="0.3";
 *    g->attrs["style"]="stroke:black;stroke-width:0.0002";
 *    s.addChild(g);
 *    for(auto it = foundEdges.begin(); it != foundEdges.end(); it++) {
 *        MS_HDS::Vertex_handle h0 = it->second->vertex();
 *        MS_HDS::Vertex_handle h1 = it->second->opposite()->vertex();
 *        double x0 = pdotv(h0->point(),pcoords[0]);
 *        double y0 = pdotv(h0->point(),pcoords[1]);
 *        double x1 = pdotv(h1->point(),pcoords[0]);
 *        double y1 = pdotv(h1->point(),pcoords[1]);
 *        g->addChild(new SVG::line(x0,y0,x1,y1));
 *        
 *        vtxpt[0] = x0;
 *        vtxpt[1] = y0;
 *        BB.expand(vtxpt);
 *        vtxpt[0] = x1;
 *        vtxpt[1] = y1;
 *        BB.expand(vtxpt);
 *    }
 * }
 */
