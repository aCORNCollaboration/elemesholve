#ifndef EIGENSPARSE_HH
#define EIGENSPARSE_HH

#include "SparseMatrix.hh"
#include <Eigen/SparseCore>
#include <Eigen/CholmodSupport>
#include <Eigen/IterativeLinearSolvers>

/// Convenience interface for Eigen sparse matrices
class EigenSparse: public MapSparse {
public:
    /// Constructor
    EigenSparse() { }
    
    /// set size for matrix
    virtual void resize(unsigned int m, unsigned int n) { M.resize(m,n); }
    /// Finalize matrix (fill sparse entries)
    virtual void finalize();
    /// Prepare solver
    virtual void setupSolver();
    /// Solve Ax = b. Resizes x as needed.
    virtual void solve(vector<double>& x, const vector<double>& b);
    
    /// Multiply b = Ax; resizes x as needed.
    virtual void mul(vector<double>& b, const vector<double>& x);
    
    virtual void display() const;
    
protected:
    typedef Eigen::SparseMatrix<double> SparseMat;
    typedef Eigen::CholmodSupernodalLLT<SparseMat, Eigen::Upper>  Solver;
    //typedef Eigen::BiCGSTAB<SparseMat>  Solver;
    SparseMat M;
    Solver slvr;
};


#endif
