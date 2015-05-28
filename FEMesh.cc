#include "FEMesh.hh"

MeshTriangle::MeshTriangle(const CDT::Face_handle& F) {
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

///////////////////////////////////////////////////////

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
    if(vertices.size() < 2) return;
        
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
                newvert.push_back(nextvert);
                break;
            }
            vtxcirc++;
        } while(ntries++ < 100);
        if(ntries == 100) { printf("Error! Unable to find next vertex!\n"); break; }
    } while( nextvert != vertices[0] );
    
    vertices = newvert;
}
    
////////////////////////////////////////////////////////
    
    
FEMesh::FEMesh(CDT& M): cdt(M), nbound(0), nfree(0) {
    for(auto it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it) vertex_enum[it] = 0;
    
    // calculate face geometry factors
    for(auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        if(!fit->is_in_domain())continue;
        meshtri.push_back(MeshTriangle(fit));
    }
}


void FEMesh::set_boundary_points(const vector<Vertex_handle>& bpts) {
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

void FEMesh::set_boundary_values(const vector<double>& v) {
    cout << "Setting boundary values...\n";
    assert(v.size() == nbound);
    bvals = v;
    tK.mul(tKh, bvals);
    for(auto it = tKh.begin(); it != tKh.end(); it++) *it *= -1;
}

void FEMesh::solve() {
    cout << "Solving system...\n";
    assert(tKh.size() == nfree);
    K.solve(cvals, tKh);
    assert(cvals.size() == nfree);
    cout << "Done.\n";
}

void FEMesh::draw() {
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