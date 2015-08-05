/// \file CellMatrix.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "CellMatrix.hh"

template<typename val_tp>
val_tp calc_minv_1(const val_tp v[2][1], val_tp p[2][2]) {
    val_tp det = v[1][0] - v[0][0];
    p[0][0] = v[1][0]/det;
    p[0][1] = v[0][0]/det;
    p[1][0] = -1/det;
    p[1][1] = 1/det;
    return det;
}

template<typename val_tp>
val_tp calc_minv_2(const val_tp v[3][2], val_tp p[3][3]) {
    const size_t P2of3[3][2][2] = { { {1,2}, {2,1} }, { {2,0}, {0,2} }, { {0,1}, {1,0} } };
    int mul = 1;
    val_tp det = 0;
    for(size_t i=0; i<3; i++) {
        p[0][i] = 0;
        for(size_t n=0; n<2; n++) {
            p[0][i] += mul*v[P2of3[i][0][n]][0]*v[P2of3[i][1][n]][1];
            mul *= -1;
        }
        det += p[0][i];
    }
    
    for(size_t i=0; i<3; i++) {
        p[0][i] /= det;
        p[1][i] = (v[P2of3[i][0][0]][1] - v[P2of3[i][1][0]][1])/det;
        p[2][i] = (v[P2of3[i][1][0]][0] - v[P2of3[i][0][0]][0])/det;
    }
    return det;
}

template<typename val_tp>
val_tp calc_minv_3(const val_tp v[4][3], val_tp p[4][4]) {
    // index permutation sets for unique 3 of 4 indices excluding 1. 3 even permutations followed by 3 odd.
    const size_t P3of4[4][3][6] = {
        { {1,2,3, 3,2,1},
          {2,3,1, 2,1,3},
          {3,1,2, 1,3,2} },
        
        { {0,2,3, 3,2,0},
          {2,3,0, 2,0,3},
          {3,0,2, 0,3,2} },
        
        { {0,1,3, 3,1,0},
          {1,3,0, 1,0,3},
          {3,0,1, 0,3,1} },
        
        { {0,1,2, 2,1,0},
          {1,2,0, 1,0,2},
          {2,0,1, 0,2,1} } };
        
    // index permutation sets for unique 2 of 4 excluding 1; alternating even/odd.
    const size_t P2of4[4][2][6] = {
        { {1,2, 2,3, 3,1},
          {2,1, 3,2, 1,3} },
        
        { {0,3, 2,0, 3,2},
          {3,0, 0,2, 2,3} },
        
        { {0,1, 1,3, 3,0},
          {1,0, 3,1, 0,3} },
        
        { {0,2, 1,0, 2,1},
          {2,0, 0,1, 1,2} } };
        
    val_tp det = 0;
    int mul = 1;
    for(size_t i=0; i<4; i++) {
        p[0][i] = 0;
        for(size_t n=0; n<6; n++) {
            p[0][i] += mul * v[P3of4[i][0][n]][0] * v[P3of4[i][1][n]][1] * v[P3of4[i][2][n]][2];
            if(n==2 || n==5) mul *= -1;
        }
        mul *= -1;
        det += p[0][i];
    }
    
    assert(mul == 1);
    for(size_t i=0; i<4; i++) {
        p[0][i] /= det;
        for(size_t j=0; j<3; j++) p[j+1][i] = 0;
        for(size_t n=0; n<6; n++) {
            p[1][i] += -mul * v[P2of4[i][0][n]][1] * v[P2of4[i][1][n]][2];
            p[2][i] += -mul * v[P2of4[i][0][n]][2] * v[P2of4[i][1][n]][0];
            p[3][i] += -mul * v[P2of4[i][0][n]][0] * v[P2of4[i][1][n]][1];
            mul *= -1;
        }
        for(size_t j=0; j<3; j++) p[j+1][i] /= det;
    }
    
    return det;
}

//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////

template< >
void CellMatrix<1,double>::calculate(const double v[2][1]) { det = calc_minv_1<double>(v,p); }
template< >
void CellMatrix<2,double>::calculate(const double v[3][2]) { det = calc_minv_2<double>(v,p); }
template< >
void CellMatrix<3,double>::calculate(const double v[4][3]) { det = calc_minv_3<double>(v,p); }

template< >
void CellMatrix<1,float>::calculate(const float v[2][1]) { det = calc_minv_1<float>(v,p); }
template< >
void CellMatrix<2,float>::calculate(const float v[3][2]) { det = calc_minv_2<float>(v,p); }
template< >
void CellMatrix<3,float>::calculate(const float v[4][3]) { det = calc_minv_3<float>(v,p); }
