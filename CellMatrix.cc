// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "CellMatrix.hh"


template< >
void CellMatrix<1>::calculate() {
    calc_vmid();
    
    det = v[1][0] - v[0][0];
    p[0][0] = v[1][0]/det;
    p[0][1] = v[0][0]/det;
    p[1][0] = -1/det;
    p[1][1] = 1/det;
}

template< >
void CellMatrix<2>::calculate() {  
    calc_vmid();

    const size_t P2of3[3][2][2] = { { {1,2}, {2,1} }, { {2,0}, {0,2} }, { {0,1}, {1,0} } };
    int mul = 1;
    det = 0;
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
}


template< >
void CellMatrix<3>::calculate() {
    calc_vmid();
    
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
        
        det = 0;
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
}
