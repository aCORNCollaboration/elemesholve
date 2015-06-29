#include "MeshSlice.hh"
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h>
#include <stdio.h>
#include "MeshVis.hh"
#include "Visr.hh"

MeshSlice::MeshSlice(const C3t3& MM, const K::Plane_3 PP, const CoordinateTransform* CT): M(MM), P(PP), myIntersectionHandler(PP, slice, CT) {
    printf("Mesh plane slice calculation!\n");
    
    // plane coordinates basis
    K::Vector_3 b0 = P.base1();
    double n0 = 1./sqrt(b0.squared_length());
    pcoords[0] = K::Vector_3(b0.x()*n0, b0.y()*n0, b0.z()*n0);
    K::Vector_3 b1 = P.base2();
    double n1 = 1./sqrt(b1.squared_length());
    pcoords[1] = K::Vector_3(b1.x()*n1, b1.y()*n1, b1.z()*n1);
    K::Vector_3 b2 = P.orthogonal_vector();
    double n2 = 1./sqrt(b2.squared_length());
    pcoords[2] = K::Vector_3(b2.x()*n2, b2.y()*n2, b2.z()*n2); 
    
    // search for one cell intersecting plane
    toProcess.insert(find_first_intersection());
    if(*toProcess.begin() == NULL) {
        printf("Non-intersecting slice plane!\n");
        return;
    }
    
    // "seed fill" on connected cells
    int ntricells = 0, nquadcells = 0;
    while(toProcess.size()) {
        // move from "to process" to "found" list
        auto it = toProcess.begin();
        Tr::Cell_handle C = *it;
        toProcess.erase(it);
        foundCells.insert(C);
        
        // collect all intersections with cell edges
        set< MS_HDS::Vertex_handle > ixn_vertices;
        find_intersections(C, ixn_vertices);
        ntricells += (ixn_vertices.size()==3);
        nquadcells += (ixn_vertices.size()==4);
        
        if(ixn_vertices.size() < 3) continue;   // non-intersecting or degenerate tip, edge intersection
        vector<Tr::Vertex_handle> traversed = make_face(C, ixn_vertices);
        
        // check which of neighboring cells share a traversal vertex
        for(int i=0; i<4; i++) {
            Tr::Cell_handle C2 = C->neighbor(i);
            if(M.is_in_complex(C2) && !foundCells.count(C2)) {
                for(auto it = traversed.begin(); it != traversed.end(); it++)
                    if(C2->has_vertex(*it)) { toProcess.insert(C2); break; }
            }
        }
    }
    
    printf("\nFound %zu cells (%i+%i) intersecting plane.\n", foundCells.size(), ntricells, nquadcells);
}

void MeshSlice::find_intersections(const Tr::Cell_handle& C, set< MS_HDS::Vertex_handle >& ixn_vertices) {
    for(int i = 0; i < 3; i++) {
        for(int j = i+1; j < 4; j++) {
            MS_HDS::Vertex_handle v = myIntersectionHandler.edgeIntersectsPlane(C,i,j);
            if(v != NULL) {
                ixn_vertices.insert(v);
                if(myIntersectionHandler.degeneracy == 2) ixn_vertices.insert(myIntersectionHandler.other_vtx);
            }
        }
    }
}

MS_HDS::Halfedge_handle MeshSlice::new_edge(MS_HDS::Vertex_handle vfrom, MS_HDS::Vertex_handle vto) {
    MSEdge e(vfrom, vto);
    assert(!foundEdges.count(e));
    MS_HDS::Halfedge h1, h2;
    h1.set_vertex(vto);
    h2.set_vertex(vfrom);
    MS_HDS::Halfedge_handle h = slice.edges_push_back(h1,h2);
    foundEdges[e] = h;
    if(verbose >= 3) std::cout << "\tNew edge " << &*vfrom << " to " << &*vto << "\n";
    if(verbose >= 4) vsr::line(MSVtoV(h->vertex()), MSVtoV(h->opposite()->vertex()));
    return h;
}

MS_HDS::Halfedge_handle MeshSlice::get_opposite(MSEdge e) const {
    auto it = foundEdges.find(e);
    if(it != foundEdges.end()) return it->second->opposite();
    return NULL;
}

MS_HDS::Halfedge_handle MeshSlice::order_face_points(set<MS_HDS::Vertex_handle>& vxs, vector<MS_HDS::Vertex_handle>& v) const {
    if(vxs.size() < 3 || vxs.size() > 4) printf("Unexpected vertex list size: %zu\n", vxs.size());
    
    // search for starting dangling handle between any pair of vertices
    MS_HDS::Halfedge_handle h = NULL;
    for(auto it1 = vxs.begin(); it1 != vxs.end(); it1++) {
        auto it2 = it1;
        it2++;
        for(; it2 != vxs.end(); it2++) {
            h = get_opposite(MSEdge(*it1, *it2));
            if(h != NULL) {
                if(foundEdges.size() == 1) h = h->opposite(); // special case for first edge
                vxs.erase(it1);
                vxs.erase(it2);
                break;
            }
        }
        if(h != NULL) break;
    }
    if(h == NULL) return h;      // TODO not this! fix degenerate hits problem.
    assert(h != NULL); // make sure we've found a starting point
   
    v.push_back(h->vertex()); // starting vertex, pointed to by handle
    assert(vxs.size() == 1 || vxs.size() == 2);
    if(vxs.size() == 2) { // decide between order if two remaining
        auto it = vxs.begin();
        if(!(*it)->mySeg.intersects(v[0]->mySeg)) {
            it++; // try the other one?
            assert((*it)->mySeg.intersects(v[0]->mySeg));
        } else assert(!(*vxs.rbegin())->mySeg.intersects(v[0]->mySeg));
        v.push_back(*it);
        vxs.erase(it);
    }
    v.push_back(*vxs.begin()); // second to last
    v.push_back(h->opposite()->vertex()); // final vertex, pointing back to start
    return h;
}

Tr::Vertex_handle traversal_vertex(const MS_HDS::Halfedge_handle& h) {
    Tr3Edge e1 = h->vertex()->mySeg;
    Tr3Edge e2 = h->opposite()->vertex()->mySeg;
    if(e1.first == e2.first) return e1.first;
    if(e1.second == e2.first) return e1.second;
    if(e1.first == e2.second) return e1.first;
    if(e1.second == e2.second) return e1.second;
    //assert(false); // TODO fix why this happens on degenerate point hits.
    return NULL;
}

vector<Tr::Vertex_handle> MeshSlice::make_face(const Tr::Cell_handle& C, set<MS_HDS::Vertex_handle>& vxs) {
    vector<Tr::Vertex_handle> traversed;
    // vertex loop around cell
    vector<MS_HDS::Vertex_handle> v;
    MS_HDS::Halfedge_handle hprev = order_face_points(vxs, v);
    if(hprev == NULL) return traversed; // shouldn't happen... but just in case...
    traversed.push_back(traversal_vertex(hprev));
    // new face for this cell
    MS_HDS::Face_handle f =  slice.faces_push_back(MS_HDS::Face());
    f->myCell = C;
    f->set_halfedge(hprev);
    MS_HDS::Halfedge_handle hinit = hprev;
    if(verbose >= 3) {
        std::cout << v.size() << " vertices in " << &*C << " ...\n";
        std::cout << "\tStarting edge " << &*hinit->opposite()->vertex() << " to " << &*hinit->vertex() << "\n";
    }
    if(verbose >= 4) vsr::startRecording();
    for(size_t i = 0; i < v.size()-1; i++) {
        MS_HDS::Halfedge_handle h = get_opposite(MSEdge(v[i], v[i+1])); // check for previously-created dangling edge
        if(h == NULL) {
            h = new_edge(v[i], v[i+1]);       // else create new edge
        } else {
            if(verbose >= 3) std::cout << "\tFound edge " << &*v[i] << " to " << &*v[i+1] << " = " << &*h->vertex() << "\n";
            if(verbose >= 4 && h->vertex() != v[i+1]) { vsr::stopRecording(); vsr::pause(); }
            assert(h->vertex() == v[i+1]);             // verify we are pointing to correct vertex
        }
        h->set_face(f);         // assign to this face
        hprev->set_next(h);     // link to previous
        hprev = h;              // be previous for next
        
        Tr::Vertex_handle tv = traversal_vertex(h);
        if(tv != NULL) traversed.push_back(tv);
    }
    if(verbose >= 4) {
        vsr::stopRecording();
        vsr::pause();
    }
    hprev->set_next(hinit); // close the loop
    return traversed;
}

Tr::Cell_handle MeshSlice::find_first_intersection() {
    for(auto it = M.cells_in_complex_begin(); it != M.cells_in_complex_end(); it++) { 
        set< MS_HDS::Vertex_handle > ixn_vertices;
        find_intersections(it, ixn_vertices);
        if(ixn_vertices.size() != 3) continue;
        new_edge(*ixn_vertices.begin(), *ixn_vertices.rbegin());
        if(verbose >= 3) std::cout << "First intersecting cell: " << &*it << "\n";
        return it;
    }
    return NULL;
}

MS_HDS::Vertex_handle MeshSlice::Intersection_handler::operator()(const K::Point_3& p) {
    MS_HDS::Vertex_handle v = M.vertices_push_back(MS_HDS::Vertex(e));
    v->point() = p;
    degeneracy = 1;
    if(p == e.first->point()) {
        v->c = 0;
        degenerate_points[e.first] = v;
    } else if(p == e.second->point()) {
        v->c = 1;
        degenerate_points[e.second] = v;
    } else {
        degeneracy = 0;
        double c0 = sqrt(CGAL::squared_distance(e.first->point(), P));
        double c1 = sqrt(CGAL::squared_distance(e.second->point(), P));
        v->c = c0/(c0+c1);
        intersections[v->mySeg] = v;
    }
    if(degeneracy) printf("Degenerate point hit!\n");
    return v;
}

MS_HDS::Vertex_handle MeshSlice::Intersection_handler::operator()(const K::Segment_3& s) {
    MS_HDS::Vertex_handle v1 = M.vertices_push_back(MS_HDS::Vertex(Tr3Edge(e.first, e.first)));
    v1->point() = e.first->point();
    MS_HDS::Vertex_handle v2 =  M.vertices_push_back(MS_HDS::Vertex(Tr3Edge(e.second, e.second)));
    v2->point() = e.second->point();
    intersections[v1->mySeg] = v1;
    degenerate_points[e.first] = v1;
    intersections[v2->mySeg] = v2;
    degenerate_points[e.second] = v2;
    degeneracy = 2;
    printf("Degenerate segment hit!\n");
    other_vtx = v2;
    return v1;
}

MS_HDS::Vertex_handle MeshSlice::Intersection_handler::edgeIntersectsPlane(Tr::Cell_handle C, int k, int l) {
    // check for known degenerate points
    degeneracy = 0;
    e = Tr3Edge(C->vertex(k), C->vertex(l));
    auto it1 = degenerate_points.find(e.first);
    if(it1 != degenerate_points.end()) {
        degeneracy++;
        it1->second->mySeg = e;
    }
    auto it2 = degenerate_points.find(e.second);
    if(it2 != degenerate_points.end()) {
        if(degeneracy) other_vtx = it1->second;
        degeneracy++;
        it2->second->mySeg = e;
        return it2->second;
    }
    if(degeneracy) return it1->second;
    
    auto it = intersections.find(e);
    if(it != intersections.end()) return it->second;
    
    if(T) {
        auto ixn = intersection(P, K::Segment_3(T->mesh2phys(e.first->point()),T->mesh2phys(e.second->point())));
        if(ixn) return boost::apply_visitor(*this, *ixn);
    } else {
        auto ixn = intersection(P, K::Segment_3(e.first->point(), e.second->point()));
        if(ixn) return boost::apply_visitor(*this, *ixn);
    }
    return NULL;
}

void MeshSlice::draw() const {
    for(auto it = foundEdges.begin(); it != foundEdges.end(); it++)
        vsr::line(MSVtoV(it->second->vertex()), MSVtoV(it->second->opposite()->vertex()));
}
