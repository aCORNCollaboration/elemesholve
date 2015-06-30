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
#include <ctime>

/// Mesh vertex data structure with file I/O
template<size_t D, typename val_tp>
class MeshVertex {
public:
    /// Constructor
    MeshVertex() { }
    
    val_tp x[D];                ///< vertex position
    val_tp v;                   ///< vertex value
    
    /// load from file
    void read(istream& is) {
        is.read((char*)x, D*sizeof(x[0]));
        is.read((char*)&v, sizeof(v));
    }
};

///  Mesh cell with adjacency information
template<size_t D, typename val_tp>
class MeshCell: public CellMatrixV<D,val_tp,int32_t> {
public:
    /// Constructor
    MeshCell() { }
    int32_t c_ID[D+1];  ///< adjacent cells
    val_tp vmid[D];     ///< midcentroid position
    
    /// load from file
    void read(istream& is) { is.read((char*)CellMatrixV<D,val_tp,int32_t>::v_ID, (D+1)*sizeof(CellMatrixV<D,val_tp,int32_t>::v_ID[0])); }
};

/// Triangular/tetrahedral mesh data structure with simple point location and file I/O
template<size_t D, typename val_tp>
class SimplexMesh {
public:
    /// Constructor
    SimplexMesh() { }
    
    /// load from file
    void read(istream& is);
    
    /// locate cell containing point
    size_t locate_cell(const val_tp* x, size_t start) const;  
    
    vector< MeshVertex<D,val_tp> > vertices;   ///< mesh vertices
    vector< MeshCell<D,val_tp> > cells;        ///< mesh cells
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

template<size_t D, typename val_tp>
void SimplexMesh<D, val_tp>::read(istream& is) {
    
    if(verbose > 1) cout << "Vertex = " << sizeof(MeshVertex<D, val_tp>) << " bytes; cell = " << sizeof(MeshCell<D, val_tp>) << " bytes.\n";
    clock_t startTime = clock();
    
    // load vertices
    int32_t nvertices;
    is.read((char*)&nvertices, sizeof(nvertices));
    if(verbose) cout << "Loading " << nvertices << " vertices..." << std::endl;
    vertices.resize(nvertices);
    for(int32_t i=0; i<nvertices; i++) vertices[i].read(is);
    
    // load cells; re-calculate and re-construct adjacency information
    CellVertices<D,val_tp> CV;
    typename SimplexMesh<D, val_tp>::face_ID cellfaces[D+1];
    int32_t ncells;
    is.read((char*)&ncells, sizeof(ncells));
    if(verbose) cout << "Loading " << ncells << " cells..." << std::endl;
    cells.resize(ncells);
    size_t npaired = 0;
    for(int32_t c=0; c<ncells; c++) {       // cell cell_ID
        MeshCell<D, val_tp>& C = cells[c];
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
                MeshCell<D, val_tp>& C2 = cells[C.c_ID[i]];     // the adjacent cell
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
    
    clock_t endTime = clock();
    if(verbose > 1) cout << "Data loaded in " << (endTime - startTime)/float(CLOCKS_PER_SEC) << " seconds.\n";
}

template<size_t D, typename val_tp>
size_t SimplexMesh<D, val_tp>::locate_cell(const val_tp* x, size_t start) const {
    return 0;   // TODO
}

#endif
