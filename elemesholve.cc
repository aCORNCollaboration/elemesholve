/// \file elemesholve.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "FEMesh2.hh"
#include "Visr.hh"
#include "CalculationProcess.hh"
#include <CGAL/config.h>

#include <tbb/task_scheduler_init.h>
#include <ctime>
#include <fstream>
using std::ofstream;

void demo2D() {
    CDT cdt;
    double vscale = 100;
    
    list<Point> list_of_seeds;
    
    MeshBoundaryConditions2 bounds;
    
    // outer square boundary
    vector<Point> vsquare;
    double hw = 0.95;
    vsquare.push_back(Point(-hw,-hw));
    vsquare.push_back(Point(-hw,hw));
    vsquare.push_back(Point(hw,hw));
    vsquare.push_back(Point(hw,-hw));
    bounds.bounds.push_back(MeshBoundary(1e-9, true));
    bounds.bounds.back().insert_constraint_poly(vsquare, cdt);
    bounds.bounds.back().closed = false;
    bounds.bounds.clear();
    
    // ground top, HV bottom
    /*
    vector<Point> toppts;
    for(int j=0; j<10; j++) toppts.push_back(Point(hw*0.01 + j*0.99*hw/9, 0.9*hw));
    bounds.push_back(MeshBoundary(1e-9, false));
    bounds.back().insert_constraint_poly(toppts, cdt);
    toppts.clear();
    for(int j=0; j<10; j++) toppts.push_back(Point(hw*0.01 + j*0.99*hw/9, -0.99*hw));
    bounds.push_back(MeshBoundary(9*0.1*vscale, false));
    bounds.back().insert_constraint_poly(toppts, cdt);
    */
    
    // bands
    int nbands = 9;
    double lband = 0.08;
    double sband = 0.10;
    for(int i=0; i<nbands; i++) {
        for(int sgn = -1; sgn <= 1; sgn += 2) {
            int nbndpts = 10;
            vector<Point> bandpoints;
            for(int j=0; j<nbndpts; j++) bandpoints.push_back(Point(sgn*0.7, -((i+0.5)*sband + (j+0.5)*lband/(nbndpts-1))));  
            bounds.bounds.push_back(MeshBoundary(i*0.1*vscale, false));
            bounds.bounds.back().insert_constraint_poly(bandpoints, cdt);
        }
    }
    
    // wires
    int nwires = 40;
    for(int i=0; i<nwires; i++) {
        double l = -0.8 + 1.6*double(i)/(nwires -1);
        Point ccenter(l, 0);
        list_of_seeds.push_back(ccenter);
        bounds.bounds.push_back(MeshBoundary(0.*vscale));
        bounds.bounds.back().make_circle(ccenter, 0.002, cdt);
    }
    
    // generate mesh
    cout << "Number of vertices before meshing: " << cdt.number_of_vertices() << std::endl;
    cout << "Meshing the triangulation..." << std::endl;
    CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(), Criteria(0.125, 0.2));
    cout << "Number of vertices after meshing: " << cdt.number_of_vertices() << std::endl;
    
    bounds.calc_bvals(cdt);
    
    FEMesh2 M(cdt);
    M.draw_logz = 0.2;
    M.set_boundary_points(bounds);
    M.set_boundary_values(bounds);
    
    vsr::startRecording(true);
    vsr::clearWindow();
    M.draw();
    vsr::stopRecording();
    
    M.solve();
    
    vsr::startRecording(true);
    vsr::clearWindow();
    M.draw();
    vsr::stopRecording();
    vsr::pause();
}

void meshgen_test() {

    CalculationProcess CP;
    
    printf("Setting %zu features.\n", CP.myGeom->polylines.size());
    vsr::startRecording(true);
    vsr::clearWindow();
    for(auto it = CP.myGeom->polylines.begin(); it != CP.myGeom->polylines.end(); it++) draw_polyline(*it);
    vsr::stopRecording();
    
    clock_t startTime = clock();
    CP.gen_mesh();
    CP.refine_mesh();
    clock_t endTime = clock();
    cout << "Mesh generation completed in " << (endTime - startTime)/float(CLOCKS_PER_SEC) << " seconds.\n";
    
    //std::ofstream f("mesh.dat");
    //CGAL::set_binary_mode(f);
    //f << c3t3.triangulation();
    //cout << "Triangulation mesh to file 'mesh.dat'\n";
    
    // Find boundary points; impose boundary conditions
    CP.setup_solver();
    assert(CP.M);
    
    // Slice!
    K::Plane_3 SPx(K::Point_3(0, 0, 0),   K::Vector_3(-1,0,0));
    K::Plane_3 SPy(K::Point_3(0, 0, 0),   K::Vector_3(0,1,0));
    K::Plane_3 SPz(K::Point_3(0, 0, 0.15), K::Vector_3(0,0,1));
    FEMesh3Slice MSx(CP.c3t3, SPx, CP.myCT);
    FEMesh3Slice MSy(CP.c3t3, SPy, CP.myCT);
    FEMesh3Slice MSz(CP.c3t3, SPz, CP.myCT, 2);
    
    // Visualize
    C3t3_Vis V(CP.c3t3, CP.myCT);
    V.print();

    vsr::startRecording(true);
    vsr::clearWindow();
    vsr::setColor(0,0,1);
    MSy.draw();
    V.draw(CP.myGeom);
    vsr::stopRecording();
    //vsr::pause();
    
    startTime = clock();
    CP.M->solve();
    endTime = clock();
    cout << "Matrix solution calculated in " << (endTime - startTime)/float(CLOCKS_PER_SEC) << " seconds.\n";
    
    ofstream slicedump;
    slicedump.open("slices_xyz.dat", std::ios::out | std::ios::binary);
    
    MSx.calc_vtxvals(*CP.M);
    MSx.dump_HDS(slicedump,*CP.M);
    //
    MSy.calc_vtxvals(*CP.M);
    MSy.dump_HDS(slicedump,*CP.M);
    //
    MSz.calc_vtxvals(*CP.M);
    MSz.dump_HDS(slicedump,*CP.M);
    
    slicedump.close();
    
    ofstream meshdump;
    meshdump.open("mesh.dat", std::ios::out | std::ios::binary);
    CP.M->write(meshdump);
    meshdump.close();
    
    vsr::startRecording(true);
    vsr::clearWindow();
    vsr::setColor(0,0,1);
    MSy.draw_projection();
    vsr::stopRecording();
    vsr::pause();
}


void* mainThread(void*) {

    //mesh_demo_2D();
    
    //int n = tbb::task_scheduler_init::default_num_threads();
    //std::cout << "TBB defaulting to " << n << " threads.\n\n";
    //tbb::task_scheduler_init init(4);
    
    meshgen_test();
    
    vsr::set_kill();
    return NULL;
}

int main(int, char**) {
    std::cout << "My CGAL library is " <<  CGAL_VERSION_NR << " (1MMmmb1000)\n"; 
    std::cout << "\twhere MM is the major number release, mm is the minor number release, and " << "b is the bug fixing number release.\n\n";
   
    vsr::initWindow("elemesholve", 0.1);
    pthread_t thread;
    pthread_create( &thread, NULL, &mainThread, NULL);
    vsr::doGlutLoop();
}
