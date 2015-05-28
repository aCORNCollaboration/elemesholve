#include <stdio.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include "UmfSparse.hh"

#include <iostream>
using std::cout;
#include <vector>
using std::vector;
#include <list>
using std::list;
#include <map>
using std::map;
#include <cmath>
#include <cassert>

#include "Visr.hh"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::Point Point;
typedef CDT::Vertex_handle Vertex_handle;

/// insert (closed) polygon of constraint points
vector<Vertex_handle> insert_constraint_poly(const vector<Point>& v, CDT& cdt, bool closed) {
    vector<Vertex_handle> vv;
    for(size_t i = 0; i < v.size(); i++) {
        vv.push_back(cdt.insert(v[i]));
        if(i) cdt.insert_constraint(vv[i-1], vv[i]);
    }
    if(closed && vv.size() >= 2) cdt.insert_constraint(vv[0], vv.back());
    return vv;
}


vector<Vertex_handle> make_circle(Point ccenter, double r, CDT& cdt, int npts = 25) {
    vector<Point> vcircle;
    for(int i = 0; i < npts; i++) {
        double th = double(i)/npts*2*M_PI;
        vcircle.push_back(Point(ccenter.x()+r*cos(th), ccenter.y()+r*sin(th)));
    }
    return insert_constraint_poly(vcircle, cdt, true);
}

/// calculations related to a triangle of the mesh
class MeshTriangle {
public:
    /// Constructor from triangle
    MeshTriangle(const CDT::Face_handle& F) {
        for(int i=0; i<3; i++) vertices[i] = F->vertex(i)->point();
        
        for(int i=0; i<3; i++) {
            int j = (i+1)%3;
            int k = (i+2)%3;
            
            double x_i = vertices[i].x();
            double y_i = vertices[i].y();
            double x_j = vertices[j].x();
            double y_j = vertices[j].y();
            double x_k = vertices[k].x();
            double y_k = vertices[k].y();
            
            if(!i) Delta = x_j*y_k - x_k*y_j - x_i*y_k + x_k*y_i + x_i*y_j - x_j*y_i;
            
            kij[i][i] = ( (y_j-y_k)*(y_j-y_k) + (x_k-x_j)*(x_k-x_j) )/(2*fabs(Delta));
            kij[i][j] = ( (y_j-y_k)*(y_k-y_i) + (x_k-x_j)*(x_i-x_k) )/(2*fabs(Delta));
            kij[i][k] = ( (y_k-y_j)*(y_j-y_i) + (x_j-x_k)*(x_i-x_j) )/(2*fabs(Delta));
            
            alpha[k] = (x_i*y_j - x_j*y_i)/Delta;
            beta[k] = (y_i-y_j)/Delta;
            gamma[k] = (x_j-x_i)/Delta;
        }
    }
    
    Point vertices[3];
    double Delta;
    double alpha[3];
    double beta[3];
    double gamma[3];
    double kij[3][3];
};

/// Finite element calculations using pre-defined mesh
class FEMesh {
public:
    /// Constructor, from mesh
    FEMesh(CDT& M): cdt(M), nbound(0), nfree(0) {
        for(auto it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it) vertex_enum[it] = 0;
       
        // calculate face geometry factors
        for(auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            if(!fit->is_in_domain())continue;
            meshtri.push_back(MeshTriangle(fit));
        }
    }
    
    /// Destructor
    ~FEMesh() { }
    
    /// Enumerate vertices given bound vertex list; generate K, ~K matrices
    void set_boundary_points(const vector<Vertex_handle>& bpts) {
        cout << "Setting " << bpts.size() << " fixed boundary points...\n";
        
        // clear prior enumeration
        for(auto it = vertex_enum.begin(); it != vertex_enum.end(); it++) it->second = 0;
        // enumerate boundary points
        for(nbound = 0; nbound<bpts.size(); nbound++) {
            auto it = vertex_enum.find(bpts[nbound]);
            assert(it != vertex_enum.end());
            it->second = -int(nbound+1);
        }
        // enumerate remaining free points
        nfree = 0;
        for(auto it = vertex_enum.begin(); it != vertex_enum.end(); it++) if(it->second == 0) it->second = nfree++;
        cout << nbound << " bound + " << nfree << " free vertices = " << vertex_enum.size() << " total.\n";
        
        // generate interaction matrices
        cout << "Filling interaction matrices over " << meshtri.size() << " faces...\n";
        int f = 0; // meshtri face number
        for(auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            if(!fit->is_in_domain())continue;
            
            int vnum[3];
            for(int i=0; i<3; i++) vnum[i] = vertex_enum[fit->vertex(i)];
            
            for(int i=0; i<3; i++) {
                int j = (i+1)%3;
                int k = (i+2)%3;
                if(vnum[i] >= 0) {
                    assert(vnum[i] < nfree && vnum[j] < (int)nfree && vnum[k] < (int)nfree);
                    K(vnum[i], vnum[i]) += meshtri[f].kij[i][i];
                    if(vnum[j] >= 0) K(vnum[i], vnum[j]) += meshtri[f].kij[i][j];
                    if(vnum[k] >= 0) K(vnum[i], vnum[k]) += meshtri[f].kij[i][k];
                } else {
                    if(vnum[j] >= 0) tK(vnum[j], -vnum[i]-1) += meshtri[f].kij[j][i];
                    if(vnum[k] >= 0) tK(vnum[k], -vnum[i]-1) += meshtri[f].kij[k][i];
                }
            }
            
            f++;
        }
        K.display();
        tK.display();
        
        cout << "Preparing matrices...\n";
        K.setupSolver();
    }
    
    /// Set boundary conditions
    void set_boundary_values(const vector<double>& v) {
        cout << "Setting boundary values...\n";
        assert(v.size() == nbound);
        bvals = v;
        tK.mul(tKh, bvals);
        for(auto it = tKh.begin(); it != tKh.end(); it++) *it *= -1;
    }
    
    /// solve
    void solve() {
        cout << "Solving system...\n";
        assert(tKh.size() == nfree);
        K.solve(cvals, tKh);
        assert(cvals.size() == nfree);
        cout << "Done.\n";
    }
    
    /// Draw mesh
    void draw() {
        cout << "Drawing mesh...\n";
        int tn = 0;
        for(auto it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); it++) {
            if(!it->is_in_domain()) continue;
            vsr::startLines();
            for(int i=0; i<4; i++) {
                assert(tn < meshtri.size());
                Point P = meshtri[tn].vertices[i%3];
                double z = 0;
                int vn = vertex_enum[it->vertex(i%3)];
                if(vn >= 0 && cvals.size() == nfree) { assert(vn < (int)cvals.size()); z = cvals[vn]; }
                if(vn < 0 && bvals.size() == nbound) { assert(-vn-1 < (int)bvals.size()); z = bvals[-vn-1]; }
                vsr::vertex(vsr::vec3(P.x(), P.y(), z));
            }
            vsr::endLines();
            tn++;
        }
    }

protected:
    CDT& cdt;                                   ///< the mesh
    size_t nbound;                              ///< number of bound degrees of freedom
    size_t nfree;                               ///< number of free degrees of freedom
    map<Vertex_handle, int> vertex_enum;        ///< vertices, enumerated for matrix index
    
    vector<MeshTriangle> meshtri;               ///< mesh triangle calculations, in finite_faces in domain iterator order
    UmfSparse K;
    double Knorm, Kcond;
    UmfSparse tK;
    vector<double> tKh;                         ///< RHS forcing terms of system
    vector<double> bvals;                       ///< boundary vertex values
    vector<double> cvals;                       ///< free vertex values
};


void* mainThread(void* args) {
    CDT cdt;
    
    printf("Hello, world.\n");
    
    list<Point> list_of_seeds;
    vector<Vertex_handle> bpts;
    vector<double> bval;
    
    // outer square boundary
    vector<Point> vsquare;
    double hw = 0.95;
    vsquare.push_back(Point(-hw,-hw));
    vsquare.push_back(Point(hw,-hw));
    vsquare.push_back(Point(hw,hw));
    vsquare.push_back(Point(-hw,hw));
    vector<Vertex_handle> vvsquare = insert_constraint_poly(vsquare, cdt, true);
    bpts.insert(bpts.end(), vvsquare.begin(), vvsquare.end());
    for(size_t i=0; i<vvsquare.size(); i++) bval.push_back(0.1);
    
    // two internal wires
    Point ccenter(0.5, -0.3);
    vector<Vertex_handle> vvcirc = make_circle(ccenter, 0.2, cdt);
    list_of_seeds.push_back(ccenter);
    bpts.insert(bpts.end(), vvcirc.begin(), vvcirc.end());
    for(size_t i=0; i<vvcirc.size(); i++) bval.push_back(0.5);
    
    Point ccenter2(-0.7, 0.1);
    vector<Vertex_handle> vvcirc2 = make_circle(ccenter2, 0.002, cdt);
    list_of_seeds.push_back(ccenter2);
    bpts.insert(bpts.end(), vvcirc2.begin(), vvcirc2.end());
    for(size_t i=0; i<vvcirc2.size(); i++) bval.push_back(-0.3);

    // generate mesh
    cout << "Number of vertices before meshing: " << cdt.number_of_vertices() << std::endl;
    cout << "Meshing the triangulation..." << std::endl;
    CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(), Criteria(0.125, 0.2));
    cout << "Number of vertices after meshing: " << cdt.number_of_vertices() << std::endl;
    
    FEMesh M(cdt);
    M.set_boundary_points(bpts);
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
