// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef SIMPLEXMESH_HH
#define SIMPLEXMESH_HH

#include "CellMatrix.hh"
#include <vector>
using std::vector;
#include <algorithm>
#include <map>
using std::map;
#include <iostream>
using std::ostream;
using std::istream;
using std::pair;
using std::cout;

/// Mesh vertex data structure with file I/O
template<size_t D>
class MeshVertex {
public:
    /// Constructor
    MeshVertex() { }
    
    double x[D];                ///< vertex position
    double v;                   ///< vertex value
    
    /// load from file
    void read(istream& is) {
        is.read((char*)x, D*sizeof(x[0]));
        is.read((char*)&v, sizeof(v));
    }
};

///  Mesh cell with adjacency information
template<size_t D>
class MeshCell: public CellMatrixV<D,int64_t> {
public:
    /// Constructor
    MeshCell() { }
    size_t c_ID[D+1];   ///< adjacent cells
    double vmid[D];     ///< midcentroid position
    
    /// load from file
    void read(istream& is) { is.read((char*)CellMatrixV<D,int64_t>::v_ID, (D+1)*sizeof(CellMatrixV<D,int64_t>::v_ID[0])); }
};

/// Triangular/tetrahedral mesh data structure with simple point location and file I/O
template<size_t D>
class SimplexMesh {
public:
    /// Constructor
    SimplexMesh() { }
    
    /// load from file
    void read(istream& is);
    
    /// locate cell containing point
    size_t locate_cell(const double* x, size_t start) const;  
    
    vector< MeshVertex<D> > vertices;   ///< mesh vertices
    vector< MeshCell<D> > cells;        ///< mesh cells
    int verbose = 1;    ///< debugging verbosity level
    
protected:
    struct face_ID {
        size_t vx[D];
        bool operator<(const face_ID& r) const {
            for(int i=0; i<D; i++) {
                if(vx[i] > r.vx[i]) return false;
                if(vx[i] < r.vx[i]) return true;
            }
            return false;
        }
    };                                                  ///< cell face ID by vertex numbers
    map<face_ID, pair<size_t, size_t> > adj_cell;       ///< cell adjacency table to (c_ID, vertex number)
};

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

template<size_t D>
void SimplexMesh<D>::read(istream& is) {
    // load vertices
    int64_t nvertices;
    is.read((char*)&nvertices, sizeof(nvertices));
    if(verbose) cout << "Loading " << nvertices << " vertices..." << std::endl;
    vertices.resize(nvertices);
    for(int64_t i=0; i<nvertices; i++) vertices[i].read(is);
    
    // load cells; re-calculate and re-construct adjacency information
    CellVertices<D> CV;
    typename SimplexMesh<D>::face_ID cellfaces[D+1];
    int64_t ncells;
    is.read((char*)&ncells, sizeof(ncells));
    if(verbose) cout << "Loading " << ncells << " cells..." << std::endl;
    cells.resize(ncells);
    size_t npaired = 0;
    for(int64_t c=0; c<ncells; c++) {       // cell cell_ID
        MeshCell<D>& C = cells[c];
        C.read(is);
        for(size_t i=0; i<D+1; i++) {       // vertex number in cell
            size_t vn = C.v_ID[i];          // global vertex v_ID
            assert(vn < vertices.size());   // check valid v_ID
            for(size_t j=0; j<D; j++) CV.v[i][j] = vertices[vn].x[j];       // collect vertex positions
            for(size_t j=0; j<D+1; j++) {   // record opposing cell face vertices
                if(j==i) continue;
                cellfaces[i].vx[j>i?j-1:j] = C.v_ID[j];
            }
            // assemble adjacency information
            std::sort(cellfaces[i].vx, cellfaces[i].vx+D);      // canonical order face vertices
            auto it = adj_cell.find(cellfaces[i]);      // check if this matches previously identified face
            if(it == adj_cell.end()) {
                adj_cell[cellfaces[i]] = pair<size_t, size_t>(c,i); // save face for future pairing
                C.c_ID[i] = -1;     // mark as unpaired face
            } else {
                C.c_ID[i] = it->second.first;           // set adjacent c_ID
                assert(C.c_ID[i] < c);                  // check valid paired cell
                MeshCell<D>& C2 = cells[C.c_ID[i]];     // the adjacent cell
                C2.c_ID[it->second.second] = c;         // set adjacent c_ID
                npaired++;
                adj_cell.erase(it);
            }
        }
        CV.calc_vmid();          // midcentroid
        for(size_t j=0; j<D; j++) C.vmid[j] = CV.vmid[j];
        C.calculate(CV.v);      // cell matrix
    }
    
    if(verbose) cout << "Loading complete; " << npaired << " paired and " << adj_cell.size() << " unpaired faces." << std::endl;
}

template<size_t D>
size_t SimplexMesh<D>::locate_cell(const double* x, size_t start) const {
    return 0;   // TODO
}

#endif
