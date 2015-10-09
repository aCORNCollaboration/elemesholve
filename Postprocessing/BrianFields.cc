/// \file BrianFields.cc

/*
    g++ -O3 --std=c++11 -o BrianFields -I${MPMUTILS}/GeneralUtils/ -L${MPMUTILS}/GeneralUtils/ BrianFields.cc -lMPMGeneralUtils
*/

#include "BrianFields.hh"

#include "ProgressBar.hh"
#include "PathUtils.hh"
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
using std::string;

void load_Brian_mirror(TriCubic G[3], bool reload) {
    printf("Loading electric field interpolator...\n");
    string baseFile = "/home/mpmendenhall/Documents/aCORN/Reference/CalculatedFieldmaps/NG6M13z20-40Low";
    
    if(!reload && fileExists(baseFile+".dat")) {
        std::ifstream fin;
        fin.open((baseFile+".dat").c_str());
        for(size_t i=0; i<3; i++) G[i].read(fin);
        printf("\tDone.\n");
        return;
    }
    
    const double dz = -32.2;
    float r0[3] = {-6,-6, 20+dz};
    float r1[3] = { 6, 6, 40+dz};
    size_t dims[3] = {121, 121, 201};
    for(int i=0; i<3; i++) {
        G[i].setDimensions(dims);
        G[i].setUserRange(r0, r1);
    }
    
    size_t c[3] = {0,0,0};
    
    std::ifstream fin;
    fin.open((baseFile+".txt").c_str(), std::ifstream::in);
    string line;
    ProgressBar* PB = new ProgressBar(dims[0]*dims[1]*dims[2]);
    int nread = 0;
    while(true) {
        std::getline(fin,line);
        if(line[0] == '%') continue;
        
        float x[3];
        float E[3];
        std::istringstream in(line);
        in >> x[0] >> x[1] >> x[2] >> E[0] >> E[1] >> E[2];
        
        for(int i=0; i<3; i++) G[i].set(c, 0.01*E[i]);
        PB->update(nread++);
        if(!increment_counter<3>(c,dims)) break;
    }
    delete PB;
    
    std::ofstream fout;
    fout.open((baseFile+".dat").c_str());
    for(size_t i=0; i<3; i++) {
        printf("Saving interpolator %zu...\n",i);
        G[i].write(fout);
    }
}
