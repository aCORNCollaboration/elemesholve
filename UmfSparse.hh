/// \file UmfSparse.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef UMFSPARSE_HH
#define UMFSPARSE_HH

#include "SparseMatrix.hh"

/// Convenience interface for umfpack sparse matrices
class UmfSparse: public SparseMatrixBase {
public:
    /// Constructor
    UmfSparse(): SparseMatrixBase() { }
    
    /// Destructor
    ~UmfSparse();
    
    /// Element access operator
    virtual double& operator()(unsigned int i, unsigned int j);
    /// Const access operator
    virtual double operator()(unsigned int i, unsigned int j) const;
    /// set size for matrix
    virtual void resize(unsigned int m, unsigned int n) { }
    
    /// Sort data by row
    void sort();
    /// Prepare solver
    virtual void setupSolver();
    /// Solve Ax = b. Resizes x as needed.
    virtual void solve(vector<double>& x, const vector<double>& b);
    
    /// Multiply b = Ax; resizes x as needed.
    virtual void mul(vector<double>& b, const vector<double>& x);
    
    virtual void display() const { printf("UmfSparse (%i x %i) matrix with %zu entries.\n", mmax, nmax, Ai.size()); }
    
protected:
    
    unsigned int mmax = 0;      ///< maximum row index entered
    unsigned int nmax = 0;      ///< maximum column index entered
    
    vector<int> Ai;             ///< row indices
    vector<int> Aj;             ///< column indices
    vector<double> Ax;          ///< value at row/column
    
    vector<int> Ap;             ///< sorted row data pointers
    bool is_sorted = false;     ///< whether data arrays are in sorted order
    
    map< pair<int, int>, int > index;   ///< index of entry locations in Ax data vector
    void* Symbolic = NULL;      ///< UMFPACK "symbolic" solver component
    void* Numeric = NULL;       ///< UMFPACK "numeric" solver component
};

#endif
