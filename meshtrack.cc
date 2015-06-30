#include "SimplexMesh.hh"
#include <fstream>
using std::ifstream;
#include <unistd.h>

#include "Visr.hh"

// g++ -O3 --std=c++11 -o meshtrack -DWITH_OPENGL=1 -I${MPMUTILS}/Visualization/ -I${MPMUTILS}/Matrix/ -L${MPMUTILS}/Visualization/ meshtrack.cc CellMatrix.cc -lMPMVis -lGL -lglut -lpthread


void* mainThread(void*) {
    
    SimplexMesh<3,float> M;
    M.verbose = 3;
    ifstream is("../elemesholve-bld/mesh.dat",  std::ios::in | std::ios::binary);
    M.read(is);
    is.close();
    
    float x[3] = {0,0,0};
    M.start_cells.push_back(M.locate_cell(x));
    x[0] = 9.5; x[1] = 0; x[2] = 9.5;
    M.start_cells.push_back(M.locate_cell(x));
    x[2] = -9.5;
    M.start_cells.push_back(M.locate_cell(x));
    x[0] = -9.5;
    M.start_cells.push_back(M.locate_cell(x));
    x[2] = 9.5;
    M.start_cells.push_back(M.locate_cell(x));
    for(int i=0; i<3; i++) x[i] = 0;
    
    float dx[3] = {0,0,0};
    float b[4];
    int nlocated = 0;
    int32_t prev_located = -1;
    for(int i=0; i<1000; i++) {
        do {
            //for(int j=0; j<3; j++) x[j] = 20.0*(double(rand())/double(RAND_MAX)-0.5);
            for(int j=0; j<3; j++) dx[j] = 1.0*(double(rand())/double(RAND_MAX)-0.5);
        } while ((x[0]+dx[0])*(x[0]+dx[0]) + (x[1]+dx[1])*(x[1]+dx[1]) > 99. || fabs(x[2]+dx[2]) > 9.5);
        for(int j=0; j<3; j++) x[j] += dx[j];
        
        //const SimplexMesh<3,float>::cell_tp& C0 = M.cells[rand()%M.cells.size()];
        //for(int j=0; j<3; j++) x[j] = C0.vmid[j];
        
        vsr::startRecording();
        vsr::setColor(0,0,1,0.2);
        vsr::startLines();
        int32_t start = prev_located == -1?  M.start_cells[0] : prev_located;
        vsr::vec3 prev_vtx;
        
        int ntries = 0;
        while(start != -1) {
            const SimplexMesh<3,float>::cell_tp& C = M.cells[start];
            prev_vtx = vsr::vec3(C.vmid[0], C.vmid[1], C.vmid[2]);
            vsr::vertex(prev_vtx);
            int32_t old_cell = start;
            start = M.walk_step(old_cell, x, b);
            if(start == old_cell) break;
            if(start == -1 && ++ntries < M.start_cells.size()) start = M.start_cells[ntries];
        }
        vsr::endLines();
       
        
        prev_located = start;
        if(start != -1) nlocated++;
        else {
            vsr::setColor(1,0,0);
            vsr::line(prev_vtx, vsr::vec3(x[0], x[1], x[2]));
        }
        
        vsr::stopRecording();
    }
    cout << "Located " << nlocated << " points.\n";
    
    vsr::pause();
    
    vsr::set_kill();
    return NULL;
}

int main(int, char**) {
    
    vsr::initWindow("meshtrack", 0.1);
    pthread_t thread;
    pthread_create( &thread, NULL, &mainThread, NULL);
    vsr::doGlutLoop();
    
    return 0;
}
