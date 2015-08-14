/// \file meshtrack.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

// g++ -O3 --std=c++11 -o meshtrack -DWITH_OPENGL=1 -I../ -I${MPMUTILS}/Visualization/ -I${MPMUTILS}/Matrix/ -L${MPMUTILS}/Visualization/ meshtrack.cc ../CellMatrix.cc -lMPMVis -lGL -lglut -lpthread

#include "SimplexMesh.hh"
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <unistd.h>

#include "Visr.hh"

void scan_line(SimplexMesh<3,float>& M, const float x0[3], const float x1[3], unsigned int n, ostream& o) {
    int32_t c = M.locate_cell(x0);
    float x[3];
    for(unsigned int i=0; i<n; i++) {
        float l = float(i)/(n-1);
        for(int j=0; j<3; j++) x[j] = x0[j]*(1-l) + x1[j]*l;
        c = M.locate_cell(x,c);
        if(c==-1) {
            cout << "Tracking failed at " << x[0] << ", " << x[1] << ", " << x[2] << " !\n";
            return;
        }
        const SimplexMesh<3,float>::cell_tp& C = M.cells[c];
        float phi = phi_tot(C, x, C.vmid);
        o << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << phi;
        for(int j=0; j<3; j++) o << "\t" << C.psolved[j+1];
        o << "\n";
    }
}

struct xycoord { float x, y; };
void radial_scanpoints(float r, float dr, vector<xycoord>& v) {
    int npts = int(2*M_PI*r/dr)+1;
    for(int i=0; i<npts; i++) {
        float th = i*2*M_PI/npts;
        xycoord xy;
        xy.x = r*cos(th);
        xy.y = r*sin(th);
        v.push_back(xy);
    }
}

void* mainThread(void*) {
    
    SimplexMesh<3,float> M;
    M.verbose = 3;
    ifstream is("../../elemesholve-bld/mesh.dat",  std::ios::in | std::ios::binary);
    M.read(is);
    is.close();
    
    float x[3] = {0,0,0};
    M.start_cells.push_back(M.locate_cell(x));
    x[2] = -9.5;
    M.start_cells.push_back(M.locate_cell(x));
    x[2] = 9.5;
    M.start_cells.push_back(M.locate_cell(x));
    
    M.verbose = 0;
    ofstream o("../../elemesholve-bld/scan.txt");
    float x0[3] = { 0, 0, -8 };
    float x1[3] = { 0, 0, 8 };
    vector<xycoord> v;
    int nptsr = 11;
    for(int nr = 0; nr < nptsr; nr++) radial_scanpoints(nr*4.0/(nptsr-1), 4.0/(nptsr-1), v);
    for(auto it = v.begin(); it != v.end(); it++) {
        x0[0] = x1[0] = it->x;
        x0[1] = x1[1] = it->y;
        scan_line(M, x0, x1, 501, o);
    }
    cout << "Done!\n";
    vsr::set_kill();
    return NULL;
    
    /*
    x[0] = 9.5; x[1] = 0; x[2] = 9.5;
    M.start_cells.push_back(M.locate_cell(x));
    x[2] = -9.5;
    M.start_cells.push_back(M.locate_cell(x));
    x[0] = -9.5;
    M.start_cells.push_back(M.locate_cell(x));
    x[2] = 9.5;
    M.start_cells.push_back(M.locate_cell(x));
    for(int i=0; i<3; i++) x[i] = 0;
    */
    
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
