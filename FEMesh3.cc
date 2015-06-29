// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "FEMesh3.hh"
#include "UmfSparse.hh"
#include "EigenSparse.hh"

FEMesh3::FEMesh3(C3t3& M, const CoordinateTransform* T): 
FEMeshSolver(new EigenSparse(), new EigenSparse()) {
    // calculate face geometry factors
    for(auto it = M.cells_in_complex_begin(); it != M.cells_in_complex_end(); it++) {
        CellMatrix<3>& C = cells[&*it];
        for(int i=0; i<4; i++) {
            C.v_ID[i] = (void*)&*it->vertex(i);
            vertex_enum[C.v_ID[i]] = 0;
            K::Point_3 p = it->vertex(i)->point();
            if(T) p = T->mesh2phys(p);
            C.v[i][0] = p.x();
            C.v[i][1] = p.y();
            C.v[i][2] = p.z();
        }
        C.calculate();
    }
    cout << "FEMesh3 calculator on " << vertex_enum.size() << " vertices and " << cells.size() << " cells.\n";
}

FEMesh3::~FEMesh3() {
    if(K) delete K;
    if(tK) delete tK;
}

void FEMesh3::dump_vertex_position(const void* v, ostream& o) const {
    Tr::Vertex* vv = (Tr::Vertex*)v;
    o << vv->point().x() << "\t" << vv->point().y() << "\t" << vv->point().z();
}
