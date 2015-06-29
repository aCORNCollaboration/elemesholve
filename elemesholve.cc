#include "FEMesh2.hh"
#include "FEMesh3.hh"
#include "Visr.hh"
#include "MeshVis.hh"
#include "GeomSetup.hh"
#include "MeshSlice.hh"
#include "FEMesh3Slice.hh"
#include <CGAL/config.h>
#include <CGAL/triangulate_polyhedron.h>

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

    //SphereTestGeom G;
    YSqueezeTransform YSQ(0, 10., 0.7);
    EMirrorGeom G(&YSQ);
    YSQ.z0 = G.myWorld.WC.gridz;
    
    GeomDomainFunctionWrapper GW(G.theWorld);
    GeomDomainMeshsizeWrapper RadiusMesh(G.theWorld, 0.5);
    
    Mesh_domain domain(GW, K::Sphere_3(CGAL::ORIGIN, 250), 1e-6);
    G.add_features(domain);

    Edge_criterea edge_criterea(RadiusMesh);
    Facet_criteria facet_criteria(15,           // angle bound
                                  RadiusMesh,   // radius bound field
                                  3e-2,         // distance bound
                                  CGAL::FACET_VERTICES_ON_SURFACE      // facet topology requirement
                                 );
    Cell_criteria cell_criteria(5,              // radius-edge ratio
                                RadiusMesh      // sizing field
                               );
    Mesh_criteria criteria(edge_criterea, facet_criteria, cell_criteria);
    
    // Mesh generation
    clock_t startTime = clock();
    printf("Generating mesh...\n");
    //C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());
    printf("Perturbing mesh...\n");
    CGAL::perturb_mesh_3(c3t3, domain, time_limit = 5);
    printf("Exuding mesh...\n");
    CGAL::exude_mesh_3(c3t3,12); 
    clock_t endTime = clock();
    cout << "Mesh generation completed in " << (endTime - startTime)/float(CLOCKS_PER_SEC) << " seconds.\n";
    
    // Find boundary points
    G.calc_bvals(c3t3);
    
    // Slice!
    K::Plane_3 SPx(K::Point_3(0,0,0), K::Vector_3(-1,0,0));
    K::Plane_3 SPy(K::Point_3(0,0,0), K::Vector_3(0,1,0));
    FEMesh3Slice MSx(c3t3, SPx, &YSQ);
    FEMesh3Slice MSy(c3t3, SPy, &YSQ);
    
    // Visualize
    C3t3_Vis V(c3t3, &YSQ);
    V.print();

    vsr::startRecording(true);
    vsr::clearWindow();
    vsr::setColor(0,0,1);
    MSy.draw();
    V.draw(&G);
    vsr::stopRecording();
    
    //vsr::pause();
    
    // electrostatic calculation
    FEMesh3 M(c3t3, &YSQ);
    M.set_boundary_points(G);
    M.set_boundary_values(G);
    
    startTime = clock();
    M.solve();
    endTime = clock();
    cout << "Matrix solution calculated in " << (endTime - startTime)/float(CLOCKS_PER_SEC) << " seconds.\n";
    
    /*
    ofstream vtxdump;
    vtxdump.open ("vtxdump.txt");
    M.dump_vertices(vtxdump);
    //M.dump_centers(vtxdump);
    vtxdump.close();
    */
    
    MSx.calc_vtxvals(M);
    MSx.write_svg("slice_x.svg",M);
    MSy.calc_vtxvals(M);
    MSy.write_svg("slice_y.svg",M);
    
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
