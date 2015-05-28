#include "UmfSparse.hh"
#include <suitesparse/umfpack.h>
#include <cassert>

UmfSparse::~UmfSparse() {
    if(Numeric) umfpack_di_free_numeric (&Numeric);
}

double& UmfSparse::operator()(int i, int j) {
    auto p = pair<int,int>(i,j);
    auto it = index.find(p);
    if(it != index.end()) return Ax[it->second];
    is_sorted = false;
    index[p] = Ax.size();
    Ai.push_back(i);
    Aj.push_back(j);
    Ax.push_back(0);
    if(i >= mmax) mmax = i+1;
    if(j >= nmax) nmax = j+1;
    return Ax.back();
}

double UmfSparse::operator()(int i, int j) const {
    auto it = index.find(pair<int,int>(i,j));
    if(it != index.end()) return Ax[it->second];
    return 0;
}

void UmfSparse::setupSolver() {
    assert(!Symbolic && !Numeric);
    sort();
    umfpack_di_symbolic(mmax, nmax, Ap.data(), Aj.data(), Ax.data(), &Symbolic, NULL, NULL);
    assert(Symbolic);
    umfpack_di_numeric(Ap.data(), Aj.data(), Ax.data(), Symbolic, &Numeric, NULL, NULL);
    assert(Numeric);
    umfpack_di_free_symbolic (&Symbolic);
    Symbolic = NULL;
}

void UmfSparse::solve(vector<double>& x, const vector<double>& b) {
    assert(Numeric);
    assert(b.size() == mmax);
    x.resize(nmax);
    umfpack_di_solve(UMFPACK_A, Ap.data(), Aj.data(), Ax.data(), x.data(), b.data(), Numeric, NULL, NULL);
}

void UmfSparse::mul(vector<double>& b, const vector<double>& x) {
    assert(x.size() >= nmax);
    b.clear();
    b.resize(mmax);
    for(size_t i=0; i < Ai.size(); i++) b[Ai[i]] += Ax[i]*x[Aj[i]];
}

void UmfSparse::sort() {
    if(is_sorted) return;
    is_sorted = true;
    
    Ap.resize(mmax+1);
    Ap[0] = 0;
    vector<int> newAi(Ai.size()), newAj(Aj.size());
    vector<double> newAx(Ax.size());
    
    int i = 0;
    int prevrow = 0;
    for(auto it = index.begin(); it != index.end(); it++) {
        newAi[i] = Ai[it->second];
        newAj[i] = Aj[it->second];
        newAx[i] = Ax[it->second];
        while(newAi[i] > prevrow) Ap[++prevrow] = i;
        it->second = i++;
    }
    while(prevrow < mmax) Ap[++prevrow] = i;
    
    Ai = newAi;
    Aj = newAj;
    Ax = newAx;
}
