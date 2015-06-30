#include "SimplexMesh.hh"
#include <fstream>
using std::ifstream;
#include <unistd.h>

// g++ -O3 --std=c++11 -o meshtrack -Wextra -Wpedantic -Woverloaded-virtual meshtrack.cc CellMatrix.cc
int main(int, char**) {
    
    SimplexMesh<3,float> M;
    M.verbose = 2;
    ifstream is("../elemesholve-bld/mesh.dat",  std::ios::in | std::ios::binary);
    M.read(is);
    is.close();
    
    return 0;
}
