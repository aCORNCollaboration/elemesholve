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
#include "HDS.hh"
#include "SliceHeader.hh"

void FEMesh3Slice::calc_vtxvals(const FEMesh3& F) {
    vtxvals.clear();
    for(auto it = slice.vertices_begin(); it != slice.vertices_end(); it++) {
        double z = F.vertex_value(it->mySeg.first)*(1-it->c) + F.vertex_value(it->mySeg.second)*(it->c);
        vtxvals[it] = z;
    }
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

void FEMesh3Slice::dump_HDS(ostream& o, const FEMesh3& F) const {
    
    SliceHeader<3,double> SH;
    for(int i=0; i<3; i++) {
        SH.x[i] = P.point()[i];
        for(int j=0; j<3; j++) SH.basis[i][j] = pcoords[i][j];
    }
    SH.write(o);
    
    HalfedgeDS_Builder<HDS_Vertex<3,float>, HDS_Face<4,float>, MS_HDS::Vertex_const_handle, MS_HDS::Halfedge_const_handle> B;
    
    HDS_Vertex<3,float> v;
    for(auto it = slice.vertices_begin(); it != slice.vertices_end(); it++) {
        K::Point_3 p = vertex_coordinate(it);
        v.x[0] = p.x();
        v.x[1] = p.y();
        v.x[2] = p.z();
        B.add_vertex(v,it);
    }
    
    B.enumerate_null_edge(NULL);
    for(auto it = slice.halfedges_begin(); it != slice.halfedges_end(); it++) B.enumerate_edge(it);
    for(auto it = slice.halfedges_begin(); it != slice.halfedges_end(); it++) B.setup_edge(it, it->next(), it->opposite(), it->vertex());
    
    HDS_Face<4,float> f;
    for(auto it = slice.faces_begin(); it != slice.faces_end(); it++) {
        const FEMesh3::CM& C = F.getCell(it->myCell);
        f.x[0] = C.area();
        for(int i=1; i<4; i++) f.x[i] = C.psolved[i];
        B.add_face(f,it->halfedge());
    }
    
    B.HDS.write(o);
}
