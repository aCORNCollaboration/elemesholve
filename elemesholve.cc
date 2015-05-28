#include "FEMesh.hh"
#include "Visr.hh"

void* mainThread(void* args) {
    CDT cdt;
    
    printf("Hello, world.\n");
    
    list<Point> list_of_seeds;
    vector<MeshBoundary> bounds;
    
    
    // outer square boundary
    vector<Point> vsquare;
    double hw = 0.95;
    vsquare.push_back(Point(-hw,-hw));
    vsquare.push_back(Point(hw,-hw));
    vsquare.push_back(Point(hw,hw));
    vsquare.push_back(Point(-hw,hw));
    bounds.push_back(MeshBoundary(0.1));
    bounds.back().insert_constraint_poly(vsquare, cdt);
    
    
    // two internal wires
    Point ccenter(0.5, -0.3);
    list_of_seeds.push_back(ccenter);
    bounds.push_back(MeshBoundary(0.8));
    bounds.back().make_circle(ccenter, 0.2, cdt);
    
    Point ccenter2(-0.7, 0.1);
    list_of_seeds.push_back(ccenter2);
    bounds.push_back(MeshBoundary(-0.3));
    bounds.back().make_circle(ccenter2, 0.002, cdt);
    
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
    M.set_boundary_points(boundary_verts);
    M.set_boundary_values(bval);
    
    vsr::startRecording(true);
    vsr::clearWindow();
    M.draw();
    vsr::stopRecording();
    vsr::pause();
    
    M.solve();
    
    vsr::startRecording(true);
    vsr::clearWindow();
    M.draw();
    vsr::stopRecording();
    vsr::pause();
    
    vsr::set_kill();
    return NULL;
}

int main(int, char**) {
    vsr::initWindow("elemesholve");
    pthread_t thread;
    pthread_create( &thread, NULL, &mainThread, NULL);
    vsr::doGlutLoop();
}
