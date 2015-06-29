// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef FEMESH_HH
#define FEMESH_HH

#include "SparseMatrix.hh"
#include "CellMatrix.hh"
#include <iostream>
using std::cout;
using std::ostream;
#include <map>
using std::map;

/// base class for specifying mesh boundary conditions
class MeshBoundaryConditions {
public:
    /// Constructor
    MeshBoundaryConditions() { }
    
    map<void*, double> bpts;    ///< boundary points and values by vertex ID
};

template<size_t D>
class FEMeshSolver {
public:    
    /// specify list of vertex ID's to fix as boundary points; set up corresponding matrices.
    void set_boundary_points(const MeshBoundaryConditions& M);
    /// set boundary values vector, in same order as set_boundary_points(...)
    void set_boundary_values(const MeshBoundaryConditions& M);
    /// calculate solution
    void solve();
    /// get solved vertex value
    double vertex_value(void* v) const;
    /// dump solved vertex values to file
    void dump_vertices(ostream& o) const;
    /// dump solved values at cell centers to file
    void dump_centers(ostream& o) const;
    /// find CellMatrix entry
    const CellMatrix<D>& getCell(void* C) const;
    
protected:
    /// Constructor protected; make subclass for filling out cells array.
    FEMeshSolver(SparseMatrixBase* myK, SparseMatrixBase* mytK): K(myK), tK(mytK) { }
    
    /// vertex position dump subroutine
    virtual void dump_vertex_position(const void* v, ostream& o) const = 0;
    
    map<void*, CellMatrix<D> > cells;           ///< mesh cell geometry calculations, indexed by underlying cell pointer
    size_t nbound;                              ///< number of bound degrees of freedom
    size_t nfree;                               ///< number of free degrees of freedom
    
    map<void*, int> vertex_enum;                ///< vertices, enumerated for matrix index
   
    SparseMatrixBase* K;                            ///< free vertices stiffness matrix
    double Knorm, Kcond;                        ///< norm and condition for solution system
    SparseMatrixBase* tK;                           ///< ~K boundary DF matrix
    vector<double> tKh;                         ///< RHS forcing terms of system
    
    vector<double> bvals;                       ///< boundary vertex values
    vector<double> cvals;                       ///< free vertex values
};



/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////



template<size_t D>
void FEMeshSolver<D>::set_boundary_values(const MeshBoundaryConditions& M) {
    cout << "Setting boundary values...\n";
    if(!nbound) { cout << "\ttoo easy without boundary points!\n"; return; }
    assert(M.bpts.size() == nbound);
    bvals.clear();
    for(auto itb = M.bpts.begin(); itb != M.bpts.end(); itb++) bvals.push_back(itb->second);
    tK->mul(tKh, bvals);
    for(auto it = tKh.begin(); it != tKh.end(); it++) *it *= -1;
}

template<size_t D>
void FEMeshSolver<D>::solve() {
    cout << "Solving system...\n";
    if(nbound && nfree) {
        assert(tKh.size() == nfree);
        K->solve(cvals, tKh);
        assert(cvals.size() == nfree);
    } else {
        cvals.clear();
        cvals.resize(nfree);
    }
    cout << "Setting solution...\n";
    double vv[D+1];
    for(auto it = cells.begin(); it != cells.end(); it++) {
        for(size_t i = 0; i < D+1; i++) vv[i] = vertex_value(it->second.v_ID[i]);
        it->second.set_solution(vv);
    }
    cout << "Done.\n";
}

template<size_t D>
double FEMeshSolver<D>::vertex_value(void* v) const {
    if(cvals.size() != nfree ||  bvals.size() != nbound) return 0;
    auto it = vertex_enum.find(v);
    assert(it != vertex_enum.end());
    if(it==vertex_enum.end()) return 0;
    int vn = it->second;
    if(vn >= 0) {
        assert(vn < (int)cvals.size());
        return cvals[vn];
    } else {
        assert(-vn-1 < (int)bvals.size());
        return bvals[-vn-1]; 
    }
}


template<size_t D>
void FEMeshSolver<D>::set_boundary_points(const MeshBoundaryConditions& M) {
    cout << "Setting " << M.bpts.size() << " fixed boundary points...\n";
    
    // clear prior enumeration
    for(auto it = vertex_enum.begin(); it != vertex_enum.end(); it++) it->second = 0;
    // enumerate boundary points
    nbound = 0;
    for(auto itb = M.bpts.begin(); itb != M.bpts.end(); itb++) {
        auto it = vertex_enum.find(itb->first);
        assert(it != vertex_enum.end()); // make sure this is a known vertex
        it->second = -int(++nbound);
    }
    // enumerate remaining free points
    nfree = 0;
    for(auto it = vertex_enum.begin(); it != vertex_enum.end(); it++) if(it->second == 0) it->second = nfree++;
    cout << nbound << " bound + " << nfree << " free vertices = " << vertex_enum.size() << " total.\n";
    
    tK->resize(nfree, nbound);
    K->resize(nfree, nfree);
    
    //for(int i=0; i<nfree; i++) (*tK)(i,0) = 0; // pre-fill tK column to assure correct size... TODO why does this sometimes fail without?
    
    // generate interaction matrices
    cout << "Filling interaction matrices over " << cells.size() << " cells...\n";
    for(auto fit = cells.begin(); fit != cells.end(); ++fit) {
        // matrix DF numbers for vertices in cell
        int vnum[D+1];
        for(size_t i=0; i<D+1; i++) {
            assert(vertex_enum.count(fit->second.v_ID[i]));
            vnum[i] = vertex_enum[fit->second.v_ID[i]];
        }
        // sum matrix terms
        for(size_t i=0; i<D+1; i++) { // each vertex i
            if(vnum[i] >= 0) (*K)(vnum[i], vnum[i]) += fit->second.k_ij(i,i);
            for(size_t j=i+1; j<D+1; j++) { // each other vertex j in cell with i
                double kij = fit->second.k_ij(i,j);
                if(vnum[j] < 0) {
                    if(vnum[i] >= 0) (*tK)(vnum[i], -vnum[j]-1) += kij;
                } else if(vnum[i] >= 0) {
                    (*K)(vnum[i], vnum[j]) += kij;
                    (*K)(vnum[j], vnum[i]) += kij;
                } else (*tK)(vnum[j], -vnum[i]-1) += kij;
            }
        }
    }
    
    K->finalize();
    tK->finalize();
    
    K->display();
    tK->display();
    
    if(!(nbound && nfree)) return;
    
    cout << "Preparing matrices...\n";
    K->setupSolver();
}

template<size_t D>
void FEMeshSolver<D>::dump_vertices(ostream& o) const {
    cout << "Dumping " << vertex_enum.size() << " vertices to file...\n";
    for(auto it = vertex_enum.begin(); it != vertex_enum.end(); it++) {
        dump_vertex_position(it->first, o);
        o << "\t";
        if(it->second >= 0) o << cvals[it->second] << "\t" << "0";
        else o << bvals[-it->second-1] << "\t" << "1";
        o << "\n";
    }
}

template<size_t D>
void FEMeshSolver<D>::dump_centers(ostream& o) const {
    for(auto it = cells.begin(); it != cells.end(); it++) {
        for(size_t i = 0; i < D; i++) o << it->second.vmid[i] << "\t";
        o << it->second.psolved[0] << "\n";
    }
}

template<size_t D>
const CellMatrix<D>& FEMeshSolver<D>::getCell(void* C) const {
    auto it = cells.find(C);
    assert(it != cells.end());
    return it->second;
}

#endif
