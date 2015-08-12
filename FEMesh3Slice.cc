/// \file FEMesh3Slice.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "FEMesh3Slice.hh"
#include "SVGBuilder.hh"
#include "ProgressBar.hh"
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

double FEMesh3Slice::get_vtxval(MS_HDS::Vertex_const_handle h) const {
    auto it = vtxvals.find(h);
    if(it == vtxvals.end()) { assert(false); return 0; }
    return it->second;
}

K::Point_3 FEMesh3Slice::vertex_coordinate(MS_HDS::Vertex_const_handle h) const {
    return K::Point_3(pdotv(h->point(),pcoords[0]), pdotv(h->point(),pcoords[1]), get_vtxval(h));
}

void FEMesh3Slice::draw_projection() const {
    for(auto it = foundEdges.begin(); it != foundEdges.end(); it++) {
        MS_HDS::Vertex_handle h0 = it->second->vertex();
        MS_HDS::Vertex_handle h1 = it->second->opposite()->vertex();
        vsr::line(PtoV(vertex_coordinate(h0)), PtoV(vertex_coordinate(h1)));
    }
}

/////////////////////

/// plane equations and gradient calculations for polygon representing a plane
class PolygonPlane {
public:
    /// Constructor
    PolygonPlane(const SVG::polygon* P) {
        assert(P);
        if(P->pts.size() < 3) return;
        
        // TODO select best set of 3 vertices from larger number
        
        for(int i=0; i<3; i++) {
            V.v[i][0] = P->pts[i].first;
            V.v[i][1] = P->pts[i].second;
        }
        V.calc_vmid();
        C.calculate(V.v);
    }
    /// Determine gradient mapping for corner values
    string gradient_remap(const vector<double>& cornervals) {
        C.set_solution(cornervals.data());
        double gx = C.psolved[1];
        double gy = C.psolved[2];
        double th = atan2(gy,gx)*180./M_PI;   // rotation angle [degrees]
        double mg2 = gx*gx + gy*gy; // magnitude of gradient
        //return "translate(-2,0) scale(2.0) rotate(45)";
        string txstr = "translate("+to_str(V.vmid[0])+","+to_str(V.vmid[1])+") ";
        txstr += "rotate("+to_str(th)+") ";
        txstr += "scale("+to_str(1./sqrt(mg2))+") ";
        txstr += "translate("+to_str(-C.psolved[0])+",0)";
        return txstr;
    }
    
    CellVertices<2,double> V;    ///< representative triangle
    CellMatrix<2,double> C;      ///< plane equations for the polygon
};

/// extended info for face polygons to allow post-processed coloring
class SVGFacePoly: public SVG::polygon {
public:
    /// Constructor
    SVGFacePoly(const string& st = ""): SVG::polygon(st) { }
    vector<double> vtxz; ///< vertex z info for coloring
};

/// Color axis
class SVGGradientAxis {
public:
    /// Constructor
    SVGGradientAxis(SVG::lingradient* lg) {
        // derived axis gradient
        Gaxis->attrs["id"] = "Gaxis";
        Gaxis->attrs["gradientTransform"] = "rotate(-90) translate(-1 0)";
        Gaxis->attrs["xlink:href"] = "#"+lg->attrs["id"];
        Gaxis->attrs["gradientUnits"] = "objectBoundingBox";
        axisGroup->addChild(Gaxis);
        
        // gradient rectangle
        SVG::rect* r = new SVG::rect(0, 0, 0.1, 1);
        r->attrs["style"] = "fill:url(#" + Gaxis->attrs["id"] + ");stroke:black;stroke-width:0.002";
        axisGroup->addChild(r);
        
        axisGroup->attrs["font-size"] = "0.07";
    }
    /// normalize to axis internal coordinates
    double axisUnits(double x) const { return (x-range.lo[0])/range.dl(0); }
    
    /// finalize range; set up text
    void finalize() {
        SVG::text* uLabel = new SVG::text(to_str(range.hi[0]), 0.107, 0.);
        axisGroup->addChild(uLabel);
        SVG::text* lLabel = new SVG::text(to_str(range.lo[0]), 0.107, 0.93);
        axisGroup->addChild(lLabel);
    }
    
    BBox<1,double> range = empty_double_bbox<1>();              ///< axis range
    SVG::group* axisGroup = new SVG::group();                   ///< group containing axis information
    XMLBuilder* Gaxis = new XMLBuilder("linearGradient");       ///< gradient
};

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
        ProgressBar* PB = new ProgressBar(slice.size_of_faces(), slice.size_of_faces()/20);
        
        // style definitions, with axis gradient
        SVG::defs* d = new SVG::defs();
        s.addChild(d);
        
        // rainbow gradient scale bar
        color::Gradient G;
        for(int i=0; i<=12; i++) {
            double l = double(i)/12;
            G.addStop(l, color::hsv((1-l)*1.5*M_PI,1,1));
        }
        SVG::lingradient* lg = new SVG::lingradient(G, "zaxis", 0, 0, 1, 0);
        lg->attrs["gradientUnits"] = "userSpaceOnUse";
        d->addChild(lg);
        SVGGradientAxis zAxis(lg);
        s.addChild(zAxis.axisGroup);
        
        // face polygons group
        SVG::group* g = new SVG::group();
        bool stroke_borders = false;
        if(stroke_borders) g->attrs["style"]="stroke:black;stroke-opacity:0.15;stroke-linejoin:round";
        s.addChild(g);
        
        vector<SVGFacePoly*> facepoly; // face polygons for color assignment
        vector<double> facez; // z values for (flat) faces
        int nfaces = 0;
        for(auto it = slice.faces_begin(); it != slice.faces_end(); it++) {
            PB->update(nfaces++);
            MS_HDS::Halfedge_const_handle h0 = it->halfedge();
            MS_HDS::Halfedge_const_handle h1 = h0;
            const FEMesh3::CM& C = F.getCell(it->myCell);
            
            // the polygon
            double w =  0.05*pow(C.area(), 1./3.);      // stroke width scales with (cell volume)^1/3 ~ cell width
            SVGFacePoly* p = new SVGFacePoly(stroke_borders? "stroke-width:"+to_str(w) : "");
            
            // collect polygon points
            size_t nOutside = 0;
            do {
                h1 = h1->next();
                K::Point_3 pvtx = h1->vertex()->point();
                for(int i=0; i<2; i++) vtxpt[i] = pdotv(pvtx,pcoords[i]) - vis_center[i];
                if(vtxpt[0]*vtxpt[0] + vtxpt[1]*vtxpt[1] > vis_rmax2) nOutside++;
                p->addpt(vtxpt[0], vtxpt[1]);
                p->vtxz.push_back(1000*get_vtxval(h1->vertex()));
                if(p->pts.size() > 4) break;
            } while(h1 != h0);
            
            // skip invalid polygons outside range
            if(!(3 <= p->pts.size() &&  p->pts.size() <= 4) || !((vis_all_inside && !nOutside) || (!vis_all_inside && nOutside < p->pts.size()))) {
                delete p;
                continue;
            }
            
            // set flat-fill style z
            if(dcmode == MAG_GRAD || dcmode == LOG_MAG_GRAD) {
                p->vtxz.clear();
                double z = C.maggrad2();
                if(dcmode == LOG_MAG_GRAD) z = log(z < grsq_min? grsq_min : z)/2.;
                else z = sqrt(z);
                zAxis.range.expand(&z);
                p->vtxz.push_back(z);
            } else if(dcmode == PHI) {
                for(auto it = p->vtxz.begin(); it != p->vtxz.end(); it++) zAxis.range.expand(&*it);
            }
            
            // expand image bounding box and insert polygon
            for(auto it = p->pts.begin(); it != p->pts.end(); it++) {
                vtxpt[0] = it->first;
                vtxpt[1] = it->second;
                BB.expand(vtxpt);
            }
            
            g->addChild(p);
            facepoly.push_back(p);
        }
        delete PB;
        
        // assign polygon colors by z
        int nsubgrad = 0;
        for(auto it = facepoly.begin(); it != facepoly.end(); it++) {
            SVGFacePoly* p = *it;
            if(!p->vtxz.size()) continue;
            
            // normalize and evaluate z info
            double zsum = 0;
            double zzsum = 0;
            for(auto itz = p->vtxz.begin(); itz != p->vtxz.end(); itz++) {
                *itz = zAxis.axisUnits(*itz);
                zsum += *itz;
                zzsum += (*itz)*(*itz);
            }
            zsum /= p->vtxz.size();
            zzsum /= p->vtxz.size();
            if(p->vtxz.size() > 1 && sqrt(zzsum - zsum*zsum) < 1e-4) {
                p->vtxz.resize(1);
                p->vtxz[0] = zsum;
            }
            // apply flat or gradient fill
            if(p->vtxz.size() == 1) {
                if(p->attrs.count("style")) p->attrs["style"] += ";";
                p->attrs["style"] += "fill:#"+color::rgb(G.hsvcolor(p->vtxz[0])).asHexString();
            } else if(p->vtxz.size() >= 3) {
                PolygonPlane PP(p);
                XMLBuilder* Gi = new XMLBuilder("linearGradient");
                Gi->attrs["id"] = lg->attrs["id"]+"_"+to_str(nsubgrad++);
                Gi->attrs["gradientTransform"] = PP.gradient_remap(p->vtxz);
                Gi->attrs["xlink:href"] = "#"+lg->attrs["id"];
                SVG::lingradient* lg = new SVG::lingradient(G, "zaxis", 0, 0, 1, 0);
                d->addChild(Gi);
                p->attrs["fill"] = "url(#" + Gi->attrs["id"] + ")";
            }
        }
        
        // scale/move axis into position
        zAxis.finalize();
        double yscale = BB.dl(1);
        zAxis.axisGroup->attrs["transform"] = "translate(" + to_str(BB.pos(1.1,0)) + " " + to_str(BB.pos(0.5,1) - 0.5*yscale) + ") scale(" + to_str(yscale) + ")";
        // expand to final display window
        vtxpt[0] = BB.pos(1.4,0);
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
