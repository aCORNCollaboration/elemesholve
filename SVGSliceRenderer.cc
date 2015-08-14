/// \file SVGSliceRenderer.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "SVGSliceRenderer.hh"
#include "SVGBuilder.hh"
#include "ProgressBar.hh"
#include "CellMatrix.hh"

/// plane equations and gradient calculations for polygon representing a plane
class PolygonPlane {
public:
    /// Constructor
    PolygonPlane(const SVG::polygon* P) {
        assert(P);
        if(P->pts.size() < 3) return;
        
        if(P->pts.size() == 4) { // try to choose best (largest) triangle for plane calcs
            size_t ibest = 0;
            double area_best = 0;
            for(size_t j=0; j<4; j++) {
                for(size_t i=0; i<4; i++) if(i != j) vtxnums.push_back(i);
                for(int i=0; i<3; i++) {
                    V.v[i][0] = P->pts[vtxnums[i]].first;
                    V.v[i][1] = P->pts[vtxnums[i]].second;
                }
                V.calc_vmid();
                C.calculate(V.v);
                vtxnums.clear();
                double area = C.area();
                if(area > area_best) { area_best = area; ibest = j; }
            }
            for(size_t i=0; i<4; i++) if(i != ibest) vtxnums.push_back(i);
        } else {
            for(size_t i=0; i<3; i++) vtxnums.push_back(i);
        }
        
        for(int i=0; i<3; i++) {
            V.v[i][0] = P->pts[vtxnums[i]].first;
            V.v[i][1] = P->pts[vtxnums[i]].second;
        }
        V.calc_vmid();
        C.calculate(V.v);
    }
    /// Determine gradient mapping for corner values
    string gradient_remap(const vector<double>& cornervals) {
        double vtxvals[3];
        for(int i=0; i<3; i++) vtxvals[i] = cornervals[vtxnums[i]];
        C.set_solution(vtxvals);
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
    
    vector<size_t> vtxnums;      ///< vertex numbers selected for plane definition
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
        //Gaxis->attrs["gradientUnits"] = "objectBoundingBox"; // should need this? but causes problems?
        axisGroup->addChild(Gaxis);
        
        // gradient rectangle
        SVG::rect* r = new SVG::rect(0, 0, 0.1, 1);
        r->attrs["style"] = "fill:url(#" + Gaxis->attrs["id"] + ");stroke:black;stroke-width:0.002";
        axisGroup->addChild(r);
        
        axisGroup->attrs["font-size"] = "0.07";
    }
    /// normalize to axis internal coordinates
    double axisUnits(double x) const { 
        if(logscale) return log(x/range.lo[0]) / log(range.hi[0]/range.lo[0]);
        return (x-range.lo[0])/range.dl(0);
    }
    
    /// finalize range; set up text
    void finalize() {
        if(logscale) {
            if(range.lo[0] < 1e-6*range.hi[0]) range.lo[0] = 1e-6*range.hi[0];
        }
        SVG::text* uLabel = new SVG::text(to_str(range.hi[0]), 0.107, 0.06);
        uLabel->attrs["dominant-baseline"]="middle";
        axisGroup->addChild(uLabel);
        SVG::text* lLabel = new SVG::text(to_str(range.lo[0]), 0.107, 0.995);
        lLabel->attrs["dominant-baseline"]="middle";
        axisGroup->addChild(lLabel);
    }
    
    bool logscale = false;                                      ///< log scale setting
    BBox<1,double> range = empty_double_bbox<1>();              ///< axis range
    SVG::group* axisGroup = new SVG::group();                   ///< group containing axis information
    XMLBuilder* Gaxis = new XMLBuilder("linearGradient");       ///< gradient
};

void SVGSliceRenderer::write_svg(const string& fname) const {
    // output file
    std::ofstream o;
    o.open (fname);
    SVG::make_standalone_header(o);
    
    // top-level SVG element
    XMLBuilder::indent = "\t";
    SVG::svg s;
    s.addChild(new SVG::title("elemesholve slice"));
    
    // style definitions group
    SVG::defs* d = new SVG::defs();
    s.addChild(d);
    
    // rainbow gradient scale bar
    color::Gradient G;
    int ngradstops = 6;
    for(int i=0; i<ngradstops; i++) {
        double l = double(i)/(ngradstops-1);
        G.addStop(l, color::hsv((1-l)*1.5*M_PI,1,1));
    }
    SVG::lingradient* lg = new SVG::lingradient(G, "zaxis", 0, 0, 1, 0);
    lg->attrs["gradientUnits"] = "userSpaceOnUse";
    d->addChild(lg);
    SVGGradientAxis zAxis(lg);
    zAxis.logscale = logscale;
    s.addChild(zAxis.axisGroup);
    
    // face polygons group
    SVG::group* g = new SVG::group();
    //bool stroke_borders = false;
    //if(stroke_borders) g->attrs["style"]="stroke:black;stroke-opacity:0.15;stroke-linejoin:round";
    s.addChild(g);
    
    BBox<2,double> BB = empty_double_bbox<2>(); // image range bounding box
    
    double vtxpt[3];
    vector<SVGFacePoly*> facepoly; // face polygons for color assignment
    int nfaces = 0;
    ProgressBar* PB = new ProgressBar(2*faces.size(), faces.size()/20);
    cout << "Drawing " << faces.size() << " faces...\n";
    for(auto it = faces.begin(); it != faces.end(); it++) {
        PB->update(nfaces++);
        //const FEMesh3::CM& C = F.getCell(it->myCell);
        
        // the polygon
        //double w =  0.05*pow(C.area(), 1./3.);      // stroke width scales with (cell volume)^1/3 ~ cell width
        SVGFacePoly* p = new SVGFacePoly(); //stroke_borders? "stroke-width:"+to_str(w) : "");
        
        // collect polygon points
        index_tp start_edge = it->edge;
        index_tp current_edge = start_edge;
        size_t nOutside = 0;
        do {
            current_edge = edges[current_edge].next;
            assert(current_edge < edges.size());
            for(int i=0; i<3; i++) vtxpt[i] = vertices[edges[current_edge].vtx].x[i];
            for(int i=0; i<2; i++) vtxpt[i] -= vis_center[i];
            if(vtxpt[0]*vtxpt[0] + vtxpt[1]*vtxpt[1] > vis_rmax2) nOutside++;
            p->addpt(vtxpt[0]*outcoord_scale, vtxpt[1]*outcoord_scale);
            p->vtxz.push_back(vtxpt[2]);
            if(p->pts.size() > 4) break;
        } while(current_edge != start_edge);
        
        // skip invalid polygons outside range
        if(!(3 <= p->pts.size() &&  p->pts.size() <= 4) || !((vis_all_inside && !nOutside) || (!vis_all_inside && nOutside < p->pts.size()))) {
            delete p;
            continue;
        }
        
        // collect appropriate data
        if(dcmode == PHI) for(auto it = p->vtxz.begin(); it != p->vtxz.end(); it++) zAxis.range.expand(&*it);
        else {
            p->vtxz.clear();
            double z = 0;
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
            zAxis.range.expand(&z);
            p->vtxz.push_back(z);
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
       
    // assign polygon colors by z
    zAxis.finalize();
    int nsubgrad = 0;
    for(auto it = facepoly.begin(); it != facepoly.end(); it++) {
        PB->update(nfaces++);
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
            p->vtxz[0] = p->vtxz[0]<0? 0 : p->vtxz[0]>1? 1 : p->vtxz[0];
            p->attrs["fill"] = "#"+color::rgb(G.hsvcolor(p->vtxz[0])).asHexString();
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
    //rsvg-convert -f pdf -o slice_y.pdf slice_y.svg
    //inkscape slice_y.svg --export-pdf=slice_y.pdf
    o.close();
}




/*
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
*/
