/// \file SVGAnalyticalRenderer.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "SVGAnalyticalRenderer.hh"
#include "aCORN_EMirror_Field.h"
#include "ProgressBar.hh"

void load_aCORN_simfield(SVGSliceRenderer& S) {
    
    struct aCORN_EMirror M;
    init_aCORN(&M);
    double x[3];
    double E[3];
    
    // new vertex values
    cout << "Calculating " << S.vertices.size() << " potential values...\n";
    ProgressBar* PB = new ProgressBar(S.vertices.size(), S.vertices.size()/20);
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
    ProgressBar* PB2 = new ProgressBar(S.faces.size(), S.faces.size()/20);
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
        for(int i=0; i<3; i++) it->x[i+1] = -E[i];
    }
    delete PB2;
    
}

