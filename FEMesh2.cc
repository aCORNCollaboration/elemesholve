// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "FEMesh2.hh"
#include "UmfSparse.hh"
#include "Visr.hh"

////////////////////////////////////////////////////////

void MeshBoundary::insert_constraint_poly(const vector<Point>& v, CDT& cdt) {
    for(size_t i = 0; i < v.size(); i++) {
        vertices.push_back(cdt.insert(v[i]));
        if(i) cdt.insert_constraint(vertices[i-1], vertices[i]);
    }
    if(closed && vertices.size() >= 2) cdt.insert_constraint(vertices[0], vertices.back());
}

void MeshBoundary::make_circle(Point ccenter, double r, CDT& cdt, int npts) {
    vector<Point> vcircle;
    closed = true;
    for(int i = 0; i < npts; i++) {
        double th = double(i)/npts*2*M_PI;
        vcircle.push_back(Point(ccenter.x()+r*cos(th), ccenter.y()+r*sin(th)));
    }
    insert_constraint_poly(vcircle, cdt);
}

bool samedirection(const Vec2& a, const Vec2& b, double epsilon = 0.0001) {
    double n = sqrt((a.x()*a.x()+a.y()*a.y())*(b.x()*b.x()+b.y()*b.y()));
    double dd = (a.x()*b.x() + a.y()*b.y())/n;
    return fabs(dd-1) < epsilon;
}

void MeshBoundary::find_refined_boundary(const CDT& cdt) {
    if(vertices.size() < 2) {
        printf("Warning: original boundary degenerate. Skipping.\n");
        return;
    }
    
    vector<Vertex_handle> newvert;
    newvert.push_back(vertices[0]);
    Vec2 dv;
    int origvert = 0;
    Vertex_handle nextvert = vertices[0];
    do {
        if(nextvert == vertices[origvert]) {
            origvert++;
            if(!closed && nextvert == vertices.back()) break;
            dv = vertices[origvert%vertices.size()]->point() - nextvert->point();
        }
        
        auto vtxcirc = cdt.incident_vertices(nextvert);
        int ntries = 0;
        do {
            Vec2 dx = vtxcirc->point() - nextvert->point();
            if(samedirection(dx,dv)) {
                nextvert = vtxcirc;
                if(nextvert != vertices[0]) newvert.push_back(nextvert);
                break;
            }
            vtxcirc++;
        } while(ntries++ < 100);
        if(ntries == 100) { printf("Error! Unable to find next vertex!\n"); break; }
    } while( nextvert != vertices[0] );
    
    vertices = newvert;
}

////////////////////////////////////////////////////////

void MeshBoundaryConditions2::calc_bvals(const CDT& cdt) {
    bpts.clear();
    for(auto it = bounds.begin(); it != bounds.end(); it++) {
        it->find_refined_boundary(cdt);
        for(auto itv = it->vertices.begin(); itv != it->vertices.end(); itv++) {
            bpts[*itv] = it->boundaryValue(*itv) ;
        }
    }
}

/////////////////////////////////////////////////////////

FEMesh2::FEMesh2(CDT& M): FEMeshSolver(new UmfSparse(), new UmfSparse()) {
    // calculate face geometry factors
    CellVertices<2> CV;
    for(auto fit = M.finite_faces_begin(); fit != M.finite_faces_end(); ++fit) {
        if(!fit->is_in_domain()) continue;
        trcells[fit] = cells.size();
        cells.push_back(CM());
        CM& C = cells.back();
        for(int i=0; i<3; i++) {
            C.v_ID[i] = fit->vertex(i);
            vertex_enum[C.v_ID[i]] = 0;
            CV.v[i][0] = fit->vertex(i)->point().x();
            CV.v[i][1] = fit->vertex(i)->point().y();
        }
        CV.calc_vmid();
        C.calculate(CV.v);
    }
    cout << "FEMesh2 calculator on " << vertex_enum.size() << " vertices and " << cells.size() << " faces.\n";
}

FEMesh2::~FEMesh2() {
    delete K;
    delete tK;
}

void FEMesh2::draw() {
    cout << "Drawing mesh...\n";
    for(auto it = cells.begin(); it != cells.end(); it++) {
        vsr::startLines();
        for(int i=0; i<4; i++) {
            double z = vertex_value(it->v_ID[i%3]);
            if(draw_logz) z = z>0? draw_logz*log(z) : 0;
            assert(false); // TODO get vmid back
            //vsr::vertex(vsr::vec3(it->v[i%3][0]+it->vmid[0], it->v[i%3][1]+it->vmid[1], z));
        }
        vsr::endLines();
    }
}

void FEMesh2::dump_vertex_position(const CDT::Vertex_handle v, ostream& o) const {
    o << v->point().x() << "\t" << v->point().y();
}

/////////////////////////////////////////////////////////////

void mesh_demo_2D() {
    
    CDT cdt;
    
    list<Point> list_of_seeds;
    MeshBoundaryConditions2 bounds;
    
    // outer square boundary
    vector<Point> vsquare;
    double hw = 0.95;
    vsquare.push_back(Point(-hw,-hw));
    vsquare.push_back(Point(-hw,hw));
    vsquare.push_back(Point(hw,hw));
    vsquare.push_back(Point(hw,-hw));
    bounds.bounds.push_back(MeshBoundary(0.0, true));
    bounds.bounds.back().insert_constraint_poly(vsquare, cdt);
    
    // two circles
    Point ccenter(-0.3, 0.5);
    list_of_seeds.push_back(ccenter);
    bounds.bounds.push_back(MeshBoundary(1.0));
    bounds.bounds.back().make_circle(ccenter, 0.2, cdt);
    
    Point ccenter2(0.5, -0.3);
    list_of_seeds.push_back(ccenter2);
    bounds.bounds.push_back(MeshBoundary(-1.0));
    bounds.bounds.back().make_circle(ccenter2, 0.002, cdt);
    
    // generate mesh
    cout << "Number of vertices before meshing: " << cdt.number_of_vertices() << std::endl;
    cout << "Meshing the triangulation..." << std::endl;
    CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(), Criteria(0.125, 0.2));
    cout << "Number of vertices after meshing: " << cdt.number_of_vertices() << std::endl;
    
    bounds.calc_bvals(cdt);
    
    FEMesh2 M(cdt);
    M.set_boundary_points(bounds);
    M.set_boundary_values(bounds);
    
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
}