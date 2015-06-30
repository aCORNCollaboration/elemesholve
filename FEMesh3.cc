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
    CellVertices<3,double> CV;
    for(auto it = M.cells_in_complex_begin(); it != M.cells_in_complex_end(); it++) {
        trcells[it] = cells.size();
        cells.push_back(CM());
        CM& C = cells.back();
        for(int i=0; i<4; i++) {
            C.v_ID[i] = it->vertex(i);
            vertex_enum[C.v_ID[i]] = 0;
            K::Point_3 p = it->vertex(i)->point();
            if(T) p = T->mesh2phys(p);
            CV.v[i][0] = p.x();
            CV.v[i][1] = p.y();
            CV.v[i][2] = p.z();
        }
        CV.calc_vmid();
        C.calculate(CV.v);
    }
    cout << "FEMesh3 calculator on " << vertex_enum.size() << " vertices and " << cells.size() << " cells.\n";
}

FEMesh3::~FEMesh3() {
    if(K) delete K;
    if(tK) delete tK;
}

const typename FEMesh3::CM& FEMesh3::getCell(const Tr::Cell_handle& C) const {
    auto it = trcells.find(C);
    assert(it != trcells.end());
    return cells[it->second];
}

void FEMesh3::dump_vertex_position(const Tr::Vertex_handle v, ostream& o) const {
    o << v->point().x() << "\t" << v->point().y() << "\t" << v->point().z();
}

void FEMesh3::write(ostream& o) const {
    cout << "Outputting solved mesh...\n";
    // vertices
    int32_t n = vxnums().size();
    o.write((char*)&n, sizeof(n));
    map<Tr::Vertex_handle,int32_t> vtxenum;
    float vx[4];
    n = 0;
    for(auto it = vxnums().begin(); it !=  vxnums().end(); it++) {
        vtxenum[it->first] = n++;
        K::Point_3 p = it->first->point();
        vx[0] = p.x();
        vx[1] = p.y();
        vx[2] = p.z();
        vx[3] = vertex_value(it->first);
        o.write((char*)vx, 4*sizeof(vx[0]));
    }
    // cells
    n = getCells().size();
    o.write((char*)&n, sizeof(n));
    int32_t vxi[4];
    for(auto it = getCells().begin(); it !=  getCells().end(); it++) {
        for(int i=0; i<4; i++) vxi[i] = vtxenum[it->v_ID[i]];
        o.write((char*)vxi, 4*sizeof(vxi[0]));
    }
    cout << "Done.\n";
}