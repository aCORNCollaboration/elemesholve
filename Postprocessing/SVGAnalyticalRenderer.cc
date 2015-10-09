/// \file SVGAnalyticalRenderer.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "SVGAnalyticalRenderer.hh"
#include "aCORN_EMirror_Field.h"
#include "BrianFields.hh"
#include "ProgressBar.hh"
#include <algorithm>

void load_aCORN_simfield(SVGSliceRenderer& S) {
    
    struct aCORN_EMirror M;
    init_aCORN_params(&M);
    M.wire_radius = 0.01;
    init_aCORN_calcs(&M);
    double x[3];
    
    // new vertex values
    cout << "Calculating " << S.vertices.size() << " potential values...\n";
    ProgressBar* PB = new ProgressBar(S.vertices.size());
    int nvert = 0;
    for(auto it = S.vertices.begin(); it != S.vertices.end(); it++) {
        PB->update(nvert++);
        for(int i=0; i<3; i++) x[i] = S.SH.x[i] + it->x[0]*S.SH.basis[0][i] + it->x[1]*S.SH.basis[1][i];
        x[2] -= 3*M.wire_radius;
        it->x[2] = calc_aCORN_potential(&M, x);
        //if(!(nvert % 1000)) cout << "x = " << x[0] << ", " << x[1] << ", " << x[2] << "\tphi = " << it->x[2] << "\n";
    }
    delete PB;
    
    // new gradients 
    cout << "Calculating " << S.faces.size() << " field points...\n";
    ProgressBar* PB2 = new ProgressBar(S.faces.size());
    double E[3];
    int nfaces = 0;
    for(auto it = S.faces.begin(); it != S.faces.end(); it++) {
        PB2->update(nfaces++);
        
        if(it == S.faces.begin()) continue;
        
        // calculate face midpoint
        for(int i=0; i<3; i++) x[i] = 0;
        index_tp e0 = it->edge;
        index_tp e = e0;
        int nvtx = 0;
        do {
            SVGSliceRenderer::Vertex& v = S.vertices[S.edges[e].vtx];
            for(int i=0; i<3; i++) x[i] += S.SH.x[i] + v.x[0]*S.SH.basis[0][i] + v.x[1]*S.SH.basis[1][i];
            e = S.edges[e].next;
            nvtx++;
        } while (e != e0);
        for(int i=0; i<3; i++) x[i] /= nvtx;
        
        x[2] -= 3*M.wire_radius;
        calc_aCORN_field(&M, x, E);
        for(int i=0; i<3; i++) it->x[i+1] = -E[i]; // note minus sign, since we want grad phi = -E
        //if(!(nfaces % 1000)) cout << "x = " << x[0] << ", " << x[1] << ", " << x[2] << "\tE = " << E[0] << ", " << E[1] << ", " << E[2] << "\n";
    }
    delete PB2;
    
}

void load_Brian_simfield(SVGSliceRenderer& S) {
    TriCubic G[3];
    load_Brian_mirror(G);
    
    float x[3];
    // new gradients 
    cout << "Calculating " << S.faces.size() << " field points...\n";
    ProgressBar* PB = new ProgressBar(S.faces.size());
    int nfaces = 0;
    for(auto it = S.faces.begin(); it != S.faces.end(); it++) {
        PB->update(nfaces++);
        
        if(it == S.faces.begin()) continue;
        
        // calculate face midpoint
        for(int i=0; i<3; i++) x[i] = 0;
        index_tp e0 = it->edge;
        index_tp e = e0;
        int nvtx = 0;
        do {
            SVGSliceRenderer::Vertex& v = S.vertices[S.edges[e].vtx];
            for(int i=0; i<3; i++) x[i] += S.SH.x[i] + v.x[0]*S.SH.basis[0][i] + v.x[1]*S.SH.basis[1][i];
            e = S.edges[e].next;
            nvtx++;
        } while (e != e0);
        for(int i=0; i<3; i++) x[i] /= nvtx;
        
        for(int i=0; i<3; i++) it->x[i+1] = G[i](x);
    }
    delete PB;
}
