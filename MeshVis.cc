// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "MeshVis.hh"
#include <iostream>
using std::cout;

void draw_Tr3Edge(const Tr3Edge& L, const CoordinateTransform* T) {
    const K::Point_3& p1 = L.first->point();
    const K::Point_3& p2 = L.second->point();
    if(!T) vsr::line(PtoV(p1) , PtoV(p2));
    else vsr::line(PtoV(T->mesh2phys(p1)) , PtoV(T->mesh2phys(p2)));
}

void draw_polyhedron(const Polyhedron& P) {
    for(auto it = P.edges_begin(); it != P.edges_end(); it++) { 
        vsr::vec3 p1 = PtoV(it->vertex()->point());
        vsr::vec3 p2 = PtoV(it->opposite()->vertex()->point());
        vsr::line(p1 , p2);
    }
}

void draw_facet(const Facet& F) {
    for(int i=0; i<3; i++) {
        if(i==F.second) continue;
        for(int j=i+1; j<4; j++) {
            if(j==F.second) continue;
            vsr::vec3 p1 = PtoV(F.first->vertex(i)->point());
            vsr::vec3 p2 = PtoV(F.first->vertex(j)->point());
            vsr::line(p1 , p2);
        }
    }
}

void draw_cell(const Cell& C) {
    for(int i=0; i<3; i++) {
        for(int j=i+1; j<4; j++) {
            vsr::vec3 p1 = PtoV(C.vertex(i)->point());
            vsr::vec3 p2 = PtoV(C.vertex(j)->point());
            vsr::line(p1 , p2);
        }
    }
}

void draw_polyline(const Polyline_3& p) {
    vsr::startLines();
    for(auto it = p.begin(); it != p.end(); it++) vsr::vertex(PtoV(*it));
    vsr::endLines();
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////

C3t3_Vis::C3t3_Vis(const C3t3& C, const CoordinateTransform* CT, bool di): draw_interior(di), T(CT) {
    for(auto it = C.facets_in_complex_begin(); it != C.facets_in_complex_end(); it++) {
        for(int i=0; i<3; i++) {
            if(i==it->second) continue;
            for(int j=i+1; j<4; j++) {
                if(j==it->second) continue;
                C3t3::Vertex_handle vi = it->first->vertex(i);
                C3t3::Vertex_handle vj = it->first->vertex(j);
                if(C.is_in_complex(vi,vj)) features.insert(Tr3Edge(vi,vj));
                else surface.insert(Tr3Edge(vi,vj));
            }
        }
    }
    
    if(!draw_interior) return;
    for(auto it = C.cells_in_complex_begin(); it != C.cells_in_complex_end(); it++) {
        for(int i=0; i<3; i++) {
            for(int j=i+1; j<4; j++) {
                Tr3Edge L(it->vertex(i), it->vertex(j));
                if(!surface.count(L)) interior.insert(L);
            }
        }
    }
}

void C3t3_Vis::draw(const MeshBoundaryConditions* B) const {
    cout << "Drawing mesh";
    if(T) cout << " with transform";
    if(B) cout << " with boundary points";
    cout << "...\n";
        
    if(draw_interior) {
        vsr::setColor(0,0,1,0.2);
        for(auto it = interior.begin(); it != interior.end(); it++) draw_Tr3Edge(*it, T);
    }
    
    if(B && B->bpts.size()) {
        for(auto it = surface.begin(); it != surface.end(); it++) {
            if(B->bpts.count(&*(it->first)) && B->bpts.count(&*(it->second))) vsr::setColor(1,0,0,0.5);
            else vsr::setColor(0,0,0,0.1);
            draw_Tr3Edge(*it, T);
        }
    } else {
        vsr::setColor(1,0,0,0.5);
        for(auto it = surface.begin(); it != surface.end(); it++) draw_Tr3Edge(*it, T);
    }
    
    vsr::setColor(0,0,0);
    for(auto it = features.begin(); it != features.end(); it++) draw_Tr3Edge(*it, T);
}

void C3t3_Vis::print() const {
    printf("C3t3 with %lu feature segments, %lu surface segments, and %lu interior segments.\n", features.size(), surface.size(), interior.size());
}
