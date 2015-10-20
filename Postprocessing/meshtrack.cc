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

struct xycoord { float x, y; };
class fieldpoint {
public:
    fieldpoint() { }
    
    float l;            ///< position parameter
    float E[3];         ///< field
    float phi;          ///< potential
    
    void display() const { cout << "E = (" << E[0] << ", " << E[1] << ", " << E[2] << "); phi = " << phi << "\n"; } 
};

/// Data from scan along mesh line
class ScanLine {
public:
    /// Constructor
    ScanLine() { }
    
    float x0[3];        ///< start point
    float x1[3];        ///< end point
    
    float x[3];         ///< position temporary variable
    /// calculate position for 0 <= l <= 1 between endpoints
    void calcx(float l) { for(int i=0; i<3; i++) x[i] = x0[i] + l*(x1[i]-x0[i]); }
    /// length of scan
    float length() const { float l = 0; for(int i=0; i<3; i++) l += pow(x1[i]-x0[i],2); return sqrt(l); }
    
    size_t n;     ///< number of scan points
    vector<fieldpoint> scandat; ///< data from scan
    
    /// get scan data from mesh; return success
    bool scan(SimplexMesh<3,float>& M) {
        int32_t c = M.locate_cell(x0);
        for(unsigned int i=0; i<n; i++) {
            fieldpoint f;
            f.l = float(i)/(n-1);
            calcx(f.l);
            c = M.locate_cell(x,c);
            if(c==-1) {
                cout << "Tracking failed at " << x[0] << ", " << x[1] << ", " << x[2] << " !\n";
                return false;
            }
            const SimplexMesh<3,float>::cell_tp& C = M.cells[c];
            f.phi = phi_tot(C, x, C.vmid);
            for(int j=0; j<3; j++) f.E[j] = -C.psolved[j+1];
            scandat.push_back(f);
        }
        return true;
    }
    
    /// dump scan line to file
    void dump_scan(ostream& o) {
        for(auto const& s: scandat) {
            calcx(s.l);
            o << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << s.phi;
            for(int j=0; j<3; j++) o << "\t" << s.E[j];
            o << "\n";
        }
    }
    
    /// calculate transverse field integral over scan
    float Et_integral() {
        float I = 0;
        for(auto const& s: scandat) {
            calcx(s.l);
            float r = sqrt(x[0]*x[0] + x[1]*x[1]);
            if(r < 1e-6) continue;
            I += (s.E[0]*x[0] + s.E[1]*x[1])/r;
        }
        return I*length()/scandat.size();
    }
    
    /// Show point
    bool showPoint(SimplexMesh<3,float>& M, double l = 0.5) {
        calcx(l);
        int c = M.locate_cell(x);
        if(c==-1) {
            cout << "Tracking failed at " << x[0] << ", " << x[1] << ", " << x[2] << " !\n";
            return false;
        }
        const SimplexMesh<3,float>::cell_tp& C = M.cells[c];
        fieldpoint f;
        f.phi = phi_tot(C, x, C.vmid);
        for(int j=0; j<3; j++) f.E[j] = -C.psolved[j+1];
        f.display();
        return true;
    } 
};

/// circle of scan points spaced approximately dr at radius r
void radial_scanpoints(float r, float dr, vector<xycoord>& v) {
    int npts = int(2*M_PI*r/dr)+1;
    for(int i=0; i<npts; i++) {
        float th = npts==1? rand()*2*M_PI/RAND_MAX : i*2*M_PI/npts;
        xycoord xy;
        xy.x = r*cos(th);
        xy.y = r*sin(th);
        v.push_back(xy);
    }
}

void radial_scans() {
    SimplexMesh<3,float> M;
    ifstream is("../../elemesholve-bld/mesh.dat",  std::ios::in | std::ios::binary);
    M.read(is);
    is.close();
    
    float x[3] = {0,0,0};
    M.start_cells.push_back(M.locate_cell(x));
    x[2] = -9.5;
    M.start_cells.push_back(M.locate_cell(x));
    x[2] = 9.5;
    M.start_cells.push_back(M.locate_cell(x));
    
    vector<xycoord> v;
    int nptsr = 401;
    for(int nr = 0; nr < nptsr; nr++) radial_scanpoints(nr*4.0/(nptsr-1), 10000, v); //radial_scanpoints(nr*4.0/(nptsr-1), 4.0/(nptsr-1), v);
    
    //ofstream o("../../elemesholve-bld/scan.txt");
    ofstream oEt("../../elemesholve-bld/Et_scan.txt");
    double gridz = 0.015;
    for(auto& p: v) {
        ScanLine L;
        float r = sqrt(p.x*p.x + p.y*p.y);
        L.x0[0] = L.x1[0] = p.x;
        L.x0[1] = L.x1[1] = p.y;
        L.n = 400;
        
        cout << "\nr = " << r << "\n";
        
        L.x0[2] = gridz-7.907;
        L.x1[2] = gridz-0.107;
        float Et_inner = L.scan(M)? L.Et_integral() : -1000;
        L.scandat.clear();
        L.showPoint(M);
        
        L.x0[2] = gridz+0.093;
        L.x1[2] = gridz+12.193;
        float Et_outer = L.scan(M)? L.Et_integral() : -1000;
        L.showPoint(M);
        
        if(Et_inner != -1000 && Et_outer != -1000) {
            oEt << r << "\t" << Et_inner << "\t" << Et_outer << "\t" << Et_inner + Et_outer << "\n";
        }
        
        //L.dump_scan(o);
    }
    cout << "Done!\n";
}

void location_test() {
    SimplexMesh<3,float> M;
    M.verbose = 3;
    ifstream is("../../elemesholve-bld/mesh.dat",  std::ios::in | std::ios::binary);
    M.read(is);
    is.close();
    /*
     x [0] = 9.5; x*[1] = 0; x[2] = 9.5;
     M.start_cells.push_back(M.locate_cell(x));
     x[2] = -9.5;
     M.start_cells.push_back(M.locate_cell(x));
     x[0] = -9.5;
     M.start_cells.push_back(M.locate_cell(x));
     x[2] = 9.5;
     M.start_cells.push_back(M.locate_cell(x));
     for(int i=0; i<3; i++) x[i] = 0;
     */
    
    float x[3] = {0,0,0};
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
}




void* mainThread(void*) {
    radial_scans();
    //location_test();
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
