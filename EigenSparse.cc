#include "EigenSparse.hh"
#include <cassert>

void EigenSparse::display() const { printf("EigenSparse (%i x %i) matrix.\n", M.rows(), M.cols()); }
    
void EigenSparse::setupSolver() {
    assert(M.rows() == M.cols());
    printf("Computing solver...\n");
    slvr.compute(M);
    assert(slvr.info() == Eigen::Success);
}

void EigenSparse::finalize() {
    printf("Loading %zu matrix entries...\n", vals.size());
    M.reserve(Eigen::VectorXi::Constant(M.cols(), 20));
    for(auto it = vals.begin(); it != vals.end(); it++)  M.insert(it->first.first, it->first.second) = it->second;
    vals.clear();
}

void EigenSparse::solve(vector<double>& x, const vector<double>& b) {
    assert((int)b.size() == M.rows());
    Eigen::VectorXd bb(b.size());
    for(size_t i=0; i<b.size(); i++) bb[i] = b[i];
    
    Eigen::VectorXd xx = slvr.solve(bb);
    assert(slvr.info() == Eigen::Success);
    x.resize(M.rows());
    for(size_t i=0; i<x.size(); i++) x[i] = xx[i];
}
    
void EigenSparse::mul(vector<double>& b, const vector<double>& x) {
    assert((int)x.size() == M.cols());
    Eigen::VectorXd xx(x.size());
    for(size_t i=0; i<x.size(); i++) xx[i] = x[i];
    
    Eigen::VectorXd bb = M * xx;
    b.resize(M.rows());
    for(size_t i=0; i<b.size(); i++) b[i] = bb[i];
}
