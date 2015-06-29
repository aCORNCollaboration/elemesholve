// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef CELL_MATRIX_HH
#define CELL_MATRIX_HH

#include <stdio.h>
#include <cassert>
#include <cmath>

template<size_t D>
inline int factorial() { return D*factorial<D-1>(); }
template< >
inline int factorial<0>() { return 1; }

/// Special matrix for D-dimensional simplex cells plane equations.
/// For D=3, M =
/// 1 x1 y1 z1
/// 1 x2 y2 z2
/// 1 x3 y3 z3
/// 1 x4 y4 z4
template<size_t D>
class CellMatrix {
public:
    /// Constructor
    CellMatrix() { }
    /// Calculate values
    void calculate();
    /// print to stdout
    void display() const;
    /// value at position x for vertex function i
    double phi_i(size_t i, const double* x, bool relcenter = false) const;
    /// value at position x for vertex values y (assumes set_solution called)
    double phi(const double* x, bool relcenter = false) const;
    /// | grad phi |^2 (assumes set_solution called)
    double maggrad2() const;
    /// cell area
    double area() const { return fabs(det/factorial<D>()); }
    /// integral over one phi plane
    double int_phi() const { return fabs(det/factorial<D+1>()); }
    /// k^{\nu}_{ij} = \int_{cell} \nabla \phi_i \cdot \nabla \phi_j dx
    double k_ij(size_t i, size_t j) const;
    /// set "solved" vertex values
    void set_solution(const double* y);
    
    void* v_ID[D+1];    ///< global vertex ID numbers
    double v[D+1][D];   ///< vertex coordinates, shifted relative to centroid by calc_vmid()
    double vmid[D];     ///< cell centroid
    double det;         ///< vertex matrix determinant
    double p[D+1][D+1]; ///< plane equation coefficients matrix = M^{-1}
    double psolved[D+1];///< total plane equation for "solved" vertex values
    
protected:
    /// Calculate mid vertex; convert v to offsets relative to mid
    void calc_vmid();
};



/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////


template<size_t D>
double CellMatrix<D>::phi_i(size_t i, const double* x, bool relcenter) const {
    double y = p[0][i];
    if(relcenter) for(size_t j=0; j<D; j++) y += p[j+1][i]*x[j];
    else for(size_t j=0; j<D; j++) y += p[j+1][i]*(x[j]-vmid[j]);
    return y;
}

template<size_t D>
double CellMatrix<D>::phi(const double* x, bool relcenter) const {
    double z = psolved[0];
    if(relcenter) for(size_t i=0; i<D; i++) z += psolved[i+1]*x[i];
    else for(size_t i=0; i<D; i++) z += psolved[i+1]*(x[i]-vmid[i]);
    return z;
}

template<size_t D>
double CellMatrix<D>::maggrad2() const {
    double g = 0;
    for(size_t i=0; i<D; i++) g += psolved[i+1]*psolved[i+1];
    return g;
}

template<size_t D>
double CellMatrix<D>::k_ij(size_t i, size_t j) const {
    double g = 0;
    for(size_t k = 0; k < D; k++) g += p[k+1][i]*p[k+1][j];
    return area()*g;
}

template<size_t D>
void CellMatrix<D>::display() const {
    for(size_t i=0; i<D+1; i++) {
        printf("|\t1");
        for(size_t j=0; j<D; j++) printf("\t%g", v[i][j]);
        printf("\t|\t\t|");
        for(size_t j=0; j<D+1; j++) printf("\t%g", p[i][j]);
        printf("\t|\n");
    }
    printf("det(M) = %g;\t k_ij =\n", det);
    
    for(size_t i=0; i<D; i++) {
        printf("|\t");
        for(size_t j=0; j<D; j++) printf("%g\t", k_ij(i,j));
        printf("|\n");
    }
    
    printf("{ ");
    for(size_t i=0; i<D+1; i++) {
        printf("{ 1");
        for(size_t j=0; j<D; j++) printf(", %g", v[i][j]);
        printf(" }, ");
    }
    printf("}\n");
}

template<size_t D>
void CellMatrix<D>::calc_vmid() {
    for(size_t i=0; i<D; i++) {
        vmid[i] = 0;
        for(size_t j=0; j<D+1; j++) vmid[i] += v[j][i];
        vmid[i] *= 1./(D+1);
        for(size_t j=0; j<D+1; j++) v[j][i] -= vmid[i];
    }
}

template<size_t D>
void CellMatrix<D>::set_solution(const double* y) {
    for(size_t i=0; i<D+1; i++) {
        psolved[i] = 0;
        for(size_t j=0; j<D+1; j++) psolved[i] += y[j]*p[i][j];
    }
}



template< >
void CellMatrix<1>::calculate();
template< >
void CellMatrix<2>::calculate();
template< >
void CellMatrix<3>::calculate();

#endif
