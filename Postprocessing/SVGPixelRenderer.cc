/// \file SVGPixelRenderer.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "SVGPixelRenderer.hh"
#include "ProgressBar.hh"
#include <stdio.h>

void SVGPixelRenderer::make_grid(double x0, double y0, double x1, double y1, size_t nx, size_t ny, bool centered) {
    if(centered) {
        double dx = (x1-x0)/(nx-1);
        x0 -= 0.5*dx;
        x1 += 0.5*dx;
        double dy = (y1-y0)/(ny-1);
        y0 -= 0.5*dy;
        y1 += 0.5*dy;
    }
    
    for(size_t ix = 0; ix < nx; ix++) {
        for(size_t iy = 0; iy < ny; iy++) {
            Pixel p;
            p.lo[0] = x0 + ix*(x1-x0)/nx;
            p.hi[0] = x0 + (ix+1)*(x1-x0)/nx;
            p.lo[1] = y0 + iy*(y1-y0)/ny;
            p.hi[1] = y0 + (iy+1)*(y1-y0)/ny;
            pxls.push_back(p);
        }
    }
}

void SVGPixelRenderer::write_svg(const string& fname) {
    // output file
    std::ofstream o;
    o.open (fname);
    SVG::svg::make_standalone_header(o);
    
    // top-level SVG element
    XMLBuilder::indent = "\t";
    SVG::svg s;
    s.addChild(new SVG::title("elemesholve slice"));
    s.attrs["shape-rendering"] = "crispEdges"; // "geometricPrecision";
    
    // style definitions group
    SVG::defs* d = new SVG::defs();
    s.addChild(d);
    d->addChild(zAxis.base_gradient);
    if(draw_axis) s.addChild(zAxis.axisGroup);
    assert(zAxis.axisGroup);
    
    // pixels group
    SVG::group* g = new SVG::group();
    s.addChild(g);
    
    BBox<2,double> BB = empty_double_bbox<2>(); // image range bounding box
    
    // generate pixels
    vector<SVG::rect*> rects;
    for(auto const& p: pxls) {
        if(autoscale) zAxis.range.expand(&p.z);
        double x0[2];
        double x1[2];
        for(int i=0; i<2; i++) {
            x0[i] = p.pos(-0.05,i);
            x0[i] -= vis_center[i];
            if(i) x0[i] *= -1;
            x0[i] *= outcoord_scale;
            
            x1[i] = p.pos(1.05,i);
            x1[i] -= vis_center[i];
            if(i) x1[i] *= -1;
            x1[i] *= outcoord_scale;
        }
        BB.expand(x0);
        BB.expand(x1);
        SVG::rect* r = new SVG::rect(x0[0], x0[1], x1[0]-x0[0], x1[1]-x0[1]);
        g->addChild(r);
        rects.push_back(r);
    }
    
    // assign polygon colors by z
    zAxis.finalize();
    printf("Z axis range %g -- %g\n", zAxis.range.lo[0], zAxis.range.hi[0]);
    for(size_t i=0; i<pxls.size(); i++) {
        double z = zAxis.axisUnits(pxls[i].z);
        z = z > 0 ? (z < 1? z : 1) : 0;
        rects[i]->attrs["fill"] = "#"+color::rgb(zAxis.G.hsvcolor(z)).asHexString();
    }
            
    if(draw_axis) {
        // scale/move axis into position
        double yscale = BB.dl(1);
        zAxis.axisGroup->attrs["transform"] = "translate(" + to_str(BB.pos(1.1,0)) + " " + to_str(BB.pos(0.5,1) - 0.5*yscale) + ") scale(" + to_str(yscale) + ")";
        // expand to final display window
        double v[2] = { BB.pos(1.7,0), BB.hi[1] };
        BB.expand(v);
    }
    
    s.setView(BB, 10);
    s.write(o);
    o.close();
    
    printf("Output svg to '%s'\n", fname.c_str());
}
