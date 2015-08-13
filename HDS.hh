/// \file HDS.hh Halfedge data structure with file I/O
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef HDS_HH
#define HDS_HH

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <iostream>
using std::ostream;
using std::istream;
using std::cout;
#include <cassert>

/// HDS vertex
template<size_t D, typename val_tp>
class HDS_Vertex {
public:
    /// Constructor
    HDS_Vertex() { }
    val_tp x[D];        ///< vertex data

    /// load from file
    virtual void read(istream& is) { is.read((char*)x, D*sizeof(x[0])); }
    /// write to file
    virtual void write(ostream& o) const { o.write((char*)x, D*sizeof(x[0])); }
};

/// index type
typedef size_t index_tp;

/// HDS edge base type
class HDS_Edge {
public:
    /// Constructor
    HDS_Edge() { }
    
    index_tp next = 0;          ///< index of next edge
    index_tp opposite = 0;      ///< index of opposite edge
    index_tp vtx = 0;           ///< index of incident vertex
    
    /// load from file
    virtual void read(istream& is) {
        is.read((char*)&next, sizeof(next));
        is.read((char*)&opposite, sizeof(opposite));
        is.read((char*)&vtx, sizeof(vtx));
    }
    /// write to file
    virtual void write(ostream& o) const {
        o.write((char*)&next, sizeof(next));
        o.write((char*)&opposite, sizeof(opposite));
        o.write((char*)&vtx, sizeof(vtx));
    }
};

/// HDS edge with "richer" information
class HDS_FatEdge: public HDS_Edge {
public:
    /// Constructor
    HDS_FatEdge() { }
    
    index_tp prev = 0;          ///< index of previous edge
    index_tp face = 0;          ///< index of associated face
};

/// HDS face type
template<size_t D, typename val_tp>
class HDS_Face {
public:
    /// Constructor
    HDS_Face() { }
    
    index_tp edge = 0;  ///< index of one incident halfedge
    val_tp x[D];        ///< face data
    
    /// load from file
    virtual void read(istream& is) {
        is.read((char*)&edge, sizeof(edge));
        is.read((char*)x, D*sizeof(x[0]));
    }
    /// write to file
    virtual void write(ostream& o) const {
        o.write((char*)&edge, sizeof(edge));
        o.write((char*)x, D*sizeof(x[0]));
    }
};

/// Halfedge data structure (not intended for modification)
template<class vtx_tp, class face_tp, class edge_tp = HDS_Edge>
class HalfedgeDS {
public:
    /// Constructor
    HalfedgeDS() { }
    
    vector<vtx_tp> vertices;    ///< enumerated vertices
    vector<edge_tp> edges;      ///< enumerated halfedges; edge 0 reserved for "unpaired"
    vector<face_tp> faces;      ///< faces
    int verbose = 1;            ///< verbosity level
    
    /// load from file
    virtual void read(istream& is);
    /// write to file
    virtual void write(ostream& o) const;
};

/// Halfedge data structure with extra helper information
template<class vtx_tp, class face_tp, class edge_tp = HDS_FatEdge>
class HalfedgeFatDS: public HalfedgeDS<vtx_tp, face_tp, edge_tp> {
public:
    /// Constructor
    HalfedgeFatDS() { }
    /// load from file
    virtual void read(istream& is) { HalfedgeDS<vtx_tp, face_tp, edge_tp>::read(is); fill_FatEdges(); }
protected:
    /// compute face and reverse info in HDS_FatEdge
    virtual void fill_FatEdges();
};

/// Helper class for building HDS from alternate version
template<class vtx_tp, class face_tp, typename vtx_index, typename edge_index>
class HalfedgeDS_Builder {
public:
    /// Constructor
    HalfedgeDS_Builder() { HDS.edges.resize(1); }
    
    /// add new vertex; assign enumeration value
    void add_vertex(const vtx_tp& v, vtx_index i) { HDS.vertices.push_back(v); vtx_enum[i] = HDS.vertices.size()-1; }
    /// define null edge index
    void enumerate_null_edge(edge_index i) { edge_enum[i] = 0; }
    /// assign enumeration for edge
    void enumerate_edge(edge_index i) { HDS.edges.resize(HDS.edges.size()+1); edge_enum[i] = HDS.edges.size()-1; }
    /// get previously-enumerated edge index
    index_tp get_edge_enum(edge_index i) const {
        auto it = edge_enum.find(i);
        assert(it != edge_enum.end());
        return it->second;
    }
    /// get previously-enumerate vertex index
    index_tp get_vtx_enum(vtx_index vtx) const {
        auto it = vtx_enum.find(vtx);
        assert(it != vtx_enum.end());
        return it->second;
    }
    /// set info for an edge
    void setup_edge(edge_index i, edge_index next, edge_index opposite, vtx_index vtx) {
        HDS_Edge& E = HDS.edges[get_edge_enum(i)];
        E.next = get_edge_enum(next);
        E.opposite = get_edge_enum(opposite);
        E.vtx = get_vtx_enum(vtx);
    }
    /// add face (specified by incident edge)
    void add_face(const face_tp& f, edge_index i) { HDS.faces.push_back(f); HDS.faces.back().edge = get_edge_enum(i); }
    
    HalfedgeDS<vtx_tp, face_tp, HDS_Edge> HDS;  ///< data structure being built
    map<vtx_index, index_tp> vtx_enum;          ///< re-enumeration of vertices
    map<edge_index, index_tp> edge_enum;        ///< re-enumeration of edges
    
};

//////////////////////////////////
//////////////////////////////////

template<class vtx_tp, class face_tp, class edge_tp>
void HalfedgeDS<vtx_tp, face_tp, edge_tp>::read(istream& is) {
    // load vertices
    index_tp nvertices;
    is.read((char*)&nvertices, sizeof(nvertices));
    if(verbose) cout << "Loading " << nvertices << " vertices..." << std::endl;
    vertices.resize(nvertices);
    for(auto it = vertices.begin(); it != vertices.end(); it++) it->read(is);
    
    // load edges
    index_tp nedges;
    is.read((char*)&nedges, sizeof(nedges));
    if(verbose) cout << "Loading " << nedges << " halfedges..." << std::endl;
    edges.resize(nedges);
    for(auto it = edges.begin(); it != edges.end(); it++) it->read(is);
    
    // load face list
    index_tp nfaces;
    is.read((char*)&nfaces, sizeof(nfaces));
    if(verbose) cout << "Loading " << nfaces << " faces..." << std::endl;
    faces.resize(nfaces);
    for(auto it = faces.begin(); it != faces.end(); it++) it->read(is);
    
    if(verbose) cout << "Done!" << std::endl;
}

template<class vtx_tp, class face_tp, class edge_tp>
void HalfedgeDS<vtx_tp, face_tp, edge_tp>::write(ostream& o) const {
    if(verbose) cout << "Writing HDS to file..." << std::endl;
    
    index_tp nvertices = vertices.size();
    if(verbose) cout << "\twriting " << nvertices << " vertices..." << std::endl;
    o.write((char*)&nvertices, sizeof(nvertices));
    for(auto it = vertices.begin(); it != vertices.end(); it++) it->write(o);
    
    index_tp nedges = edges.size();
    if(verbose) cout << "\twriting " << nedges << " halfedges..." << std::endl;
    o.write((char*)&nedges, sizeof(nedges));
    for(auto it = edges.begin(); it != edges.end(); it++) it->write(o);
    
    index_tp nfaces = faces.size();
    if(verbose) cout << "\twriting " << nfaces << " faces..." << std::endl;
    o.write((char*)&nfaces, sizeof(nfaces));
    for(auto it = faces.begin(); it != faces.end(); it++) it->write(o);
    
    if(verbose) cout << "\tDone!" << std::endl;
}

template<class vtx_tp, class face_tp, class edge_tp>
void HalfedgeFatDS<vtx_tp, face_tp, edge_tp>::fill_FatEdges() {
    if(this->verbose) cout << "Calculating additional connections..." << std::endl;
    index_tp facenum = 0;
    for(auto it = this->faces.begin(); it != this->faces.end(); it++) {
        index_tp init_edge = it->edge;
        index_tp current_edge = init_edge;
        do {
            index_tp prev_edge = current_edge;
            current_edge = this->edges[current_edge].next;
            this->edges[current_edge].prev = prev_edge;
            this->edges[current_edge].face = facenum;
        } while(current_edge != init_edge);
        facenum++;
    }
}

#endif
