/// \file SimplexMesh.hh
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
#include <cfloat>

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
    size_t locate_cell(const val_tp x[D], int32_t start = -1, unsigned int max_retries = 20) const;  
    /// mesh walk step, searching for position B; calculates barycentric coordinates up to first negative with crossable face
    int32_t walk_step(int32_t c0, const val_tp p[D], val_tp bcoords[D+1]) const;
    
    vector<int32_t> start_cells;        ///< initial guess cells for point location
    
    /// shorthand for cell type
    typedef MeshCell<D,val_tp> cell_tp;
    /// shorthand for vertex type
    typedef MeshVertex<D,val_tp> vtx_tp;
    
    vector<vtx_tp> vertices;    ///< mesh vertices
    vector<cell_tp> cells;      ///< mesh cells
    int verbose = 1;            ///< debugging verbosity level
    
protected:
    /// cell face ID by vertex numbers
    struct face_ID {
        size_t vx[D];   ///< face vertices
        /// sorting operator
        bool operator<(const face_ID& r) const {
            for(int i=0; i<D; i++) {
                if(vx[i] > r.vx[i]) return false;
                if(vx[i] < r.vx[i]) return true;
            }
            return false;
        }
    };                                                  
    map<face_ID, pair<size_t, size_t> > adj_cell;       ///< cell adjacency table to (c_ID, vertex number)
};



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



template<size_t D, typename val_tp>
void SimplexMesh<D, val_tp>::read(istream& is) {
    
    if(verbose > 1) cout << "Vertex = " << sizeof(vtx_tp) << " bytes; cell = " << sizeof(cell_tp) << " bytes.\n";
    clock_t startTime = clock();
    
    // load vertices
    int32_t nvertices;
    is.read((char*)&nvertices, sizeof(nvertices));
    if(verbose) cout << "Loading " << nvertices << " vertices..." << std::endl;
    vertices.resize(nvertices);
    for(int32_t i=0; i<nvertices; i++) vertices[i].read(is);
    
    // load cells; re-calculate and re-construct adjacency information
    CellVertices<D,val_tp> CV;
    face_ID cellfaces[D+1];
    int32_t ncells;
    is.read((char*)&ncells, sizeof(ncells));
    if(verbose) cout << "Loading " << ncells << " cells..." << std::endl;
    cells.resize(ncells);
    size_t npaired = 0;
    for(int32_t c=0; c<ncells; c++) {       // cell cell_ID
        cell_tp& C = cells[c];
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
                C.c_ID[i] = it->second.first;   // set adjacent c_ID
                assert(C.c_ID[i] < c);          // check valid paired cell
                cell_tp& C2 = cells[C.c_ID[i]]; // the adjacent cell
                C2.c_ID[it->second.second] = c; // set adjacent c_ID
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
int32_t SimplexMesh<D, val_tp>::walk_step(int32_t c0, const val_tp p[D], val_tp b[D+1]) const {
    const cell_tp& C = cells[c0];
    // cell-centered coordinates
    val_tp pc[D];
    for(size_t i=0; i<D; i++) pc[i] = p[i] - C.vmid[i];
    
    float r2_min = FLT_MAX;
    int32_t cbest = -1;
    bool inside = true;
    for(size_t i=0; i<D+1; i++) { // coordinates from each vertex
        b[i] = phi_i<D,val_tp>(C, i, pc);
        if(b[i] < 0) {
            inside = false;
            if(C.c_ID[i] != -1) {
                const cell_tp& C2 = cells[C.c_ID[i]];
                float r2 = 0;
                for(int j=0; j<D; j++) r2 += (p[j]-C2.vmid[j])*(p[j]-C2.vmid[j]);
                if(r2 < r2_min) {
                    r2_min = r2;
                    cbest = C.c_ID[i];
                }
            }
        }
        if(b[i] > 1) inside = false;
    }
    return inside? c0 : cbest;
}

/*
template<size_t D, typename val_tp>
int32_t SimplexMesh<D, val_tp>::walk_step(int32_t c0, const val_tp p[D], val_tp b[D+1]) const {
    const cell_tp& C = cells[c0];
    // cell-centered coordinates
    val_tp pc[D];
    for(size_t i=0; i<D; i++) pc[i] = p[i] - C.vmid[i];
    
    bool inside = true;
    for(size_t i=0; i<D+1; i++) { // coordinates from each vertex
        b[i] = phi_i<D,val_tp>(C, i, pc);
        if(b[i] < 0) {
            inside = false;
            if(C.c_ID[i] != -1) return C.c_ID[i];
        }
        if(b[i] > 1) inside = false;
    }
    return inside? c0 : -1;
}
*/

template<size_t D, typename val_tp>
size_t SimplexMesh<D, val_tp>::locate_cell(const val_tp x[D], int32_t start, unsigned int max_retries) const {
    int ntries = 0;
    unsigned int n_retries = 0;
    
    if(verbose > 2) {
        cout << "Searching for point";
        for(int i=0; i<D; i++) cout << "\t" << x[i]; 
        cout << "\n";
    }
    
    while(n_retries < start_cells.size() || n_retries < max_retries) {
        if(start == -1) {
            if(n_retries < start_cells.size()) start = start_cells[n_retries];
            else start = rand()%cells.size();
        }
        n_retries++;

        val_tp b[D+1];
        while(start != -1) {
            ntries++;
            if(verbose > 3) {
                const cell_tp& C = cells[start];
                cout << "\t" << start;
                for(int i=0; i<D; i++) cout << "\t" << x[i]-C.vmid[i];
                cout << "\n";
            }
            int32_t old_cell = start;
            start = walk_step(old_cell, x, b);
            if(start == old_cell) {
                if(verbose > 2) cout << "\tlocated in " << ntries << " steps!\n";
                return old_cell;
            }
        }
        start = -1;
        if(verbose > 2) cout << "\tSearch failed! Starting over!\n";
    }
    
    if(verbose > 2) cout << "\tWalk search teminated without result!\n";
    return -1;
}

#endif
