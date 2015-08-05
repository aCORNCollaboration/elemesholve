/// \file CalculationProcess.hh
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef CALCULATIONPROCESS_HH
#define CALCULATIONPROCESS_HH

#include "GeomSetup.hh"
#include "MeshVis.hh"
#include "FEMesh3Slice.hh"

/// Container for classes/steps involved in specific mesh calculation
class CalculationProcess {
public:
    CalculationProcess();

    YSqueezeTransform YSQ;
    EMirrorGeom G;
    GeomDomainFunctionWrapper GW;
    GeomDomainMeshsizeWrapper RadiusMesh;
    Mesh_domain domain;
    
    Edge_criterea edge_criterea;
    Facet_criteria facet_criteria;
    Cell_criteria cell_criteria;
    Mesh_criteria criteria;
    
    C3t3 c3t3;
    
    void gen_mesh();
    void refine_mesh();
    
    FEMesh3* M = NULL;
    void setup_solver();
    
};

#endif
