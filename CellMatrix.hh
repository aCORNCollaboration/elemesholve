/// \file CellMatrix.hh D-dimensional simplex (D=1,2,3) plane equations and linear interpolation
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

/// compile-time factorial definition
template<size_t D>
inline int factorial() { return D*factorial<D-1>(); }
/// compile-time factorial definition, special 0 case
template< >
inline int factorial<0>() { return 1; }

/// Simplex cell vertex coordinates, shifted relative to centroid
template<size_t D, typename val_tp>
class CellVertices {
public:
    /// Constructor
    CellVertices() { }
    
    val_tp v[D+1][D];   ///< vertex coordinates, shifted relative to centroid by calc_vmid()
    val_tp vmid[D];     ///< cell centroid
    
    /// Calculate mid vertex; convert v to offsets relative to mid
    void calc_vmid();
};

/// Special matrix calculations \f$\det(M)\f$, \f$M^{-1}\f$ for D-dimensional simplex cells plane equations.
/** For D=3, M =
    1 x1 y1 z1
    1 x2 y2 z2
    1 x3 y3 z3
    1 x4 y4 z4 */
template<size_t D, typename val_tp>
class CellMatrix {
public:
    /// Constructor
    CellMatrix() { }
    /// Calculate matrix for vertices
    void calculate(const val_tp v[D+1][D]);
    /// print to stdout
    void display() const;
    /// \f$| \nabla \phi |^2\f$ (assumes set_solution called)
    val_tp maggrad2() const;
    /// cell area
    val_tp area() const { return fabs(det/factorial<D>()); }
    /// integral over one phi plane
    val_tp int_phi() const { return fabs(det/factorial<D+1>()); }
    /// \f$k^{\nu}_{ij} = \int_{cell} \nabla \phi_i \cdot \nabla \phi_j dx\f$
    val_tp k_ij(size_t i, size_t j) const;
    /// set "solved" vertex values
    void set_solution(const val_tp* y);
    
    val_tp det;         ///< vertex matrix determinant
    val_tp p[D+1][D+1]; ///< plane equation coefficients matrix \f$M^{-1}\f$
    val_tp psolved[D+1];///< total plane equation for "solved" vertex values
};

/// Cell matrix class with vertex IDs
template<size_t D, typename val_tp, typename vtx_id>
class CellMatrixV: public CellMatrix<D, val_tp> {
public:
    /// Constructor
    CellMatrixV() { }
    vtx_id v_ID[D+1];    ///< vertex ID numbers
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template<size_t D, typename val_tp>
void CellVertices<D, val_tp>::calc_vmid() {
    for(size_t i=0; i<D; i++) {
        vmid[i] = 0;
        for(size_t j=0; j<D+1; j++) vmid[i] += v[j][i];
        vmid[i] *= 1./(D+1);
        for(size_t j=0; j<D+1; j++) v[j][i] -= vmid[i];
    }
}

/// calculate cell pyramid basis function \f$\phi_i(x)\f$
template<size_t D, typename val_tp>
val_tp phi_i(const CellMatrix<D, val_tp>& C, size_t i, const val_tp* x, const val_tp* vmid = NULL) {
    val_tp y = C.p[0][i];
    if(!vmid) for(size_t j=0; j<D; j++) y += C.p[j+1][i]*x[j];
    else for(size_t j=0; j<D; j++) y += C.p[j+1][i]*(x[j]-vmid[j]);
    return y;
}

/// calculate cell function value \f$\phi(x)\f$
template<size_t D, typename val_tp>
val_tp phi_tot(const CellMatrix<D, val_tp>& C, const val_tp* x, const val_tp* vmid = NULL) {
    val_tp z = C.psolved[0];
    if(!vmid) for(size_t i=0; i<D; i++) z += C.psolved[i+1]*x[i];
    else for(size_t i=0; i<D; i++) z += C.psolved[i+1]*(x[i]-vmid[i]);
    return z;
}

template<size_t D, typename val_tp>
val_tp CellMatrix<D, val_tp>::maggrad2() const {
    val_tp g = 0;
    for(size_t i=0; i<D; i++) g += psolved[i+1]*psolved[i+1];
    return g;
}

template<size_t D, typename val_tp>
val_tp CellMatrix<D, val_tp>::k_ij(size_t i, size_t j) const {
    val_tp g = 0;
    for(size_t k = 0; k < D; k++) g += p[k+1][i]*p[k+1][j];
    return area()*g;
}

template<size_t D, typename val_tp>
void CellMatrix<D, val_tp>::display() const {
    for(size_t i=0; i<D+1; i++) {
        //printf("|\t1");
        //for(size_t j=0; j<D; j++) printf("\t%g", v[i][j]);
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
    /*
    printf("{ ");
    for(size_t i=0; i<D+1; i++) {
        printf("{ 1");
        for(size_t j=0; j<D; j++) printf(", %g", v[i][j]);
        printf(" }, ");
    }
    printf("}\n");
    */
}

template<size_t D, typename val_tp>
void CellMatrix<D, val_tp>::set_solution(const val_tp* y) {
    for(size_t i=0; i<D+1; i++) {
        psolved[i] = 0;
        for(size_t j=0; j<D+1; j++) psolved[i] += y[j]*p[i][j];
    }
}


//////////////////////////////////

/// calculation routine for 1D float-valued CellMatrix
template< >
void CellMatrix<1,float>::calculate(const float v[2][1]);
/// calculation routine for 2D float-valued CellMatrix
template< >
void CellMatrix<2,float>::calculate(const float v[3][2]);
/// calculation routine for 3D float-valued CellMatrix
template< >
void CellMatrix<3,float>::calculate(const float v[4][3]);

/// calculation routine for 1D double-valued CellMatrix
template< >
void CellMatrix<1,double>::calculate(const double v[2][1]);
/// calculation routine for 2D double-valued CellMatrix
template< >
void CellMatrix<2,double>::calculate(const double v[3][2]);
/// calculation routine for 3D double-valued CellMatrix
template< >
void CellMatrix<3,double>::calculate(const double v[4][3]);

#endif
