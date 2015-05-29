#include "FEMesh.hh"
#include "Visr.hh"

void demo2D() {
    CDT cdt;
    double vscale = 100;
    
    list<Point> list_of_seeds;
    vector<MeshBoundary> bounds;
    
    // outer square boundary
    vector<Point> vsquare;
    double hw = 0.95;
    vsquare.push_back(Point(-hw,-hw));
    vsquare.push_back(Point(-hw,hw));
    vsquare.push_back(Point(hw,hw));
    vsquare.push_back(Point(hw,-hw));
    bounds.push_back(MeshBoundary(1e-9, true));
    bounds.back().insert_constraint_poly(vsquare, cdt);
    bounds.back().closed = false;
    bounds.clear();
    
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
            bounds.push_back(MeshBoundary(i*0.1*vscale, false));
            bounds.back().insert_constraint_poly(bandpoints, cdt);
        }
    }
    
    // wires
    int nwires = 40;
    for(int i=0; i<nwires; i++) {
        double l = -0.8 + 1.6*double(i)/(nwires -1);
        Point ccenter(l, 0);
        list_of_seeds.push_back(ccenter);
        bounds.push_back(MeshBoundary(0.*vscale));
        bounds.back().make_circle(ccenter, 0.002, cdt);
    }
    
    // generate mesh
    cout << "Number of vertices before meshing: " << cdt.number_of_vertices() << std::endl;
    cout << "Meshing the triangulation..." << std::endl;
    CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(), Criteria(0.125, 0.2));
    cout << "Number of vertices after meshing: " << cdt.number_of_vertices() << std::endl;
    
    vector<Vertex_handle> boundary_verts;
    vector<double> bval;
    for(auto it = bounds.begin(); it != bounds.end(); it++) {
        it->find_refined_boundary(cdt);
        boundary_verts.insert(boundary_verts.end(), it->vertices.begin(), it->vertices.end());
        for(auto it2 = it->vertices.begin(); it2 != it->vertices.end(); it2++) bval.push_back(it->boundaryValue(*it2));
    }
    
    FEMesh M(cdt);
    M.draw_logz = 0.2;
    M.set_boundary_points(boundary_verts);
    M.set_boundary_values(bval);
    
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


void* mainThread(void*) {
    demo2D();
    vsr::set_kill();
    return NULL;
}

int main(int, char**) {
    vsr::initWindow("elemesholve");
    pthread_t thread;
    pthread_create( &thread, NULL, &mainThread, NULL);
    vsr::doGlutLoop();
}
