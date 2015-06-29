// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef FEMESH3_HH
#define FEMESH3_HH

#include "FEMesh.hh"
#include "CGAL_Types.hh"
#include "GeomDef.hh"

/// 3-dimensional solver class
class FEMesh3: public FEMeshSolver<3> {
public:
    /// Do-nothing constructor
    FEMesh3() { }
    /// Constructor, from 3D triangulation mesh
    FEMesh3(C3t3& M, const CoordinateTransform* T = NULL);
    /// Destructor
    ~FEMesh3();
    
protected:
    /// vertex position dump subroutine
    virtual void dump_vertex_position(const void* v, ostream& o) const;
};

#endif
