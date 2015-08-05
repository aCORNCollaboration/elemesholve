/// \file SparseMatrix.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef SparseMatrixBase_HH
#define SparseMatrixBase_HH

#include <vector>
using std::vector;
#include <stdio.h>

#include <map>
using std::map;
using std::pair;

/// Virtual base class for sparse matrix interface
class SparseMatrixBase {
public:
    /// Constructor for empty (resizeable) matrix
    SparseMatrixBase() { }
    
    /// Destructor
    virtual ~SparseMatrixBase() { }
    
    /// Element access operator
    virtual double& operator()(unsigned int i, unsigned int j) = 0;
    /// Finalize matrix (fill sparse entries)
    virtual void finalize() { }
    
    /// set size for matrix
    virtual void resize(unsigned int m, unsigned int n) = 0;
    
    /// Prepare solver on filled matrix
    virtual void setupSolver() = 0;
    /// Solve Ax = b. Resizes x as needed.
    virtual void solve(vector<double>& x, const vector<double>& b) = 0;
    
    /// Multiply b = Ax; resizes x as needed.
    virtual void mul(vector<double>& b, const vector<double>& x) = 0;
    
    virtual void display() const { printf("Sparse matrix.\n"); }
};

/// Base class for simple map-based sparse
class MapSparse: public SparseMatrixBase {
public:
    /// Constructor
    MapSparse() { }
    
    double& operator()(unsigned int i, unsigned int j) { return vals[pair<int,int>(i,j)]; }
     
protected:
    map< pair<int, int>, double > vals;
};


#endif
