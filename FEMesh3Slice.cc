/// \file FEMesh3Slice.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "FEMesh3Slice.hh"
#include "SVGBuilder.hh"
#include "BBox.hh"


void FEMesh3Slice::calc_vtxvals(const FEMesh3& F) {
    vtxvals.clear();
    for(auto it = slice.vertices_begin(); it != slice.vertices_end(); it++) {
        double z = F.vertex_value(it->mySeg.first)*(1-it->c) + F.vertex_value(it->mySeg.second)*(it->c);
        vtxvals[it] = z;
    }

    // find max/min gradients
    grsq_min = DBL_MAX;
    grsq_max = 0;
    for(auto it = slice.faces_begin(); it != slice.faces_end(); it++) {
        double grsq = F.getCell(it->myCell).maggrad2();
        if(grsq < grsq_min) grsq_min = grsq;
        if(grsq > grsq_max) grsq_max = grsq;
    }
    cout << "E_min = " << sqrt(grsq_min) << "; E_max = " << sqrt(grsq_max) << "\n";
    if(grsq_min < 1e-10*grsq_max) grsq_min = 1e-10*grsq_max;
}

double FEMesh3Slice::get_vtxval(MS_HDS::Vertex_handle h) const {
    auto it = vtxvals.find(h);
    if(it == vtxvals.end()) { assert(false); return 0; }
    return it->second;
}

double pdotv(K::Point_3 p, K::Vector_3 v) {
    return p.x()*v.x() + p.y()*v.y() + p.z()*v.z();
}
    
void FEMesh3Slice::draw_projection() const {
    for(auto it = foundEdges.begin(); it != foundEdges.end(); it++) {
        MS_HDS::Vertex_handle h0 = it->second->vertex();
        MS_HDS::Vertex_handle h1 = it->second->opposite()->vertex();
        
        vsr::vec3 v0(pdotv(h0->point(),pcoords[0]), pdotv(h0->point(),pcoords[1]), get_vtxval(h0));
        vsr::vec3 v1(pdotv(h1->point(),pcoords[0]), pdotv(h1->point(),pcoords[1]), get_vtxval(h1));
        
        vsr::line(v0, v1);
    }
}

void FEMesh3Slice::write_svg(const string& fname, const FEMesh3& F) const {
    std::ofstream o;
    o.open (fname);
    SVG::make_standalone_header(o);
    
    SVG::svg s;
   
    s.addChild(new SVG::title("elemesholve slice"));
    BBox<2,double> BB = empty_double_bbox<2>();
    double vtxpt[2];
    
    bool drawFaces = (grsq_max != grsq_min);
    if(drawFaces) {
        cout << "Drawing " << slice.size_of_faces() << " faces...\n";
        color::Gradient G;
        for(int i=0; i<=12; i++) {
            double l = double(i)/12;
            G.addStop(l, color::hsv((1-l)*1.5*M_PI,1,1));
        }
        
        SVG::defs* d = new SVG::defs();
        s.addChild(d);
        SVG::lingradient* lg = new SVG::lingradient(G, "zaxis", 0, 0, 0, 1);
        d->addChild(lg);
        
        SVG::group* g = new SVG::group();
        g->attrs["style"]="stroke:black;stroke-opacity:0.15;stroke-linejoin:round";
        s.addChild(g);
        
        for(auto it = slice.faces_begin(); it != slice.faces_end(); it++) {
            MS_HDS::Halfedge_const_handle h0 = it->halfedge();
            MS_HDS::Halfedge_const_handle h1 = h0;
            const FEMesh3::CM& C = F.getCell(it->myCell);
            
            double z = C.maggrad2();
            if(color_logz) {
                if(z < grsq_min) z = 0;
                else z = log(z/grsq_max)/log(grsq_max/grsq_min) + 1.;
            }
            else z = sqrt(z/grsq_max);
            
            double w =  0.05*pow(C.area(), 1./3.);      // stroke width scales with (cell volume)^1/3 ~ cell width
            SVG::polygon* p = new SVG::polygon("fill:#"+color::rgb(G.hsvcolor(z)).asHexString()+";stroke-width:"+to_str(w));
            size_t nOutside = 0;
            do {
                h1 = h1->next();
                K::Point_3 pvtx = h1->vertex()->point();
                for(int i=0; i<2; i++) vtxpt[i] = pdotv(pvtx,pcoords[i]) - vis_center[i];
                if(vtxpt[0]*vtxpt[0] + vtxpt[1]*vtxpt[1] > vis_rmax2) nOutside++;
                p->addpt(vtxpt[0], vtxpt[1]);
                if(p->pts.size() > 4) break;
            } while(h1 != h0);
            
            if(p->pts.size() <= 4 && ((vis_all_inside && !nOutside) || (!vis_all_inside && nOutside < p->pts.size()))) {
                for(auto it = p->pts.begin(); it != p->pts.end(); it++) {
                    vtxpt[0] = it->first;
                    vtxpt[1] = it->second;
                    BB.expand(vtxpt);
                }
                g->addChild(p);
            } else { delete p; }
        }
        
        // gradient scale bar
        SVG::rect* r = new SVG::rect(BB.pos(1.1,0), BB.lo[1], 0.1*BB.dl(0), BB.dl(1));
        s.addChild(r);
        r->attrs["fill"] = lg->idstr();
        vtxpt[0] = BB.pos(1.2,0);
        vtxpt[1] = BB.hi[1];
        BB.expand(vtxpt);
    }
    
    bool drawEdges = !drawFaces;
    if(drawEdges) {
        cout << "Drawing " << foundEdges.size() << " edges...\n";
        SVG::group* g = new SVG::group();
        //g->attrs["opacity"]="0.3";
        g->attrs["style"]="stroke:black;stroke-width:0.0002";
        s.addChild(g);
        for(auto it = foundEdges.begin(); it != foundEdges.end(); it++) {
            MS_HDS::Vertex_handle h0 = it->second->vertex();
            MS_HDS::Vertex_handle h1 = it->second->opposite()->vertex();
            double x0 = pdotv(h0->point(),pcoords[0]);
            double y0 = pdotv(h0->point(),pcoords[1]);
            double x1 = pdotv(h1->point(),pcoords[0]);
            double y1 = pdotv(h1->point(),pcoords[1]);
            g->addChild(new SVG::line(x0,y0,x1,y1));
            
            vtxpt[0] = x0;
            vtxpt[1] = y0;
            BB.expand(vtxpt);
            vtxpt[0] = x1;
            vtxpt[1] = y1;
            BB.expand(vtxpt);
        }
    }
    
    s.setView(BB, 10);
    s.write(o);
    //rsvg-convert -f pdf -o slice_y.pdf slice_y.svg
    o.close();
}
