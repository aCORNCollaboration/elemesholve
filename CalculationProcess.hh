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
    /// Constructor
    CalculationProcess();

    YSqueezeTransform YSQ;                      ///< coordinate transform for coarser meshing along wire length
    CoordinateTransform* myCT;                  ///< coordinate transform in use
    
    SphereTestGeom STG;                         ///< spherical test geometry
    EMirrorGeom EMG;                            ///< mirror geometry to mesh
    GeomSetup* myGeom;                          ///< geometry being calculated
   
    
    GeomDomainFunctionWrapper GW;               ///< wrapper for mesh domain labeling
    GeomDomainMeshsizeWrapper RadiusMesh;       ///< wrapper for mesh size info
    Mesh_domain domain;                         ///< domain to mesh
    
    Edge_criterea edge_criterea;                ///< mesh edge constraints
    Facet_criteria facet_criteria;              ///< mesh facet constraints
    Cell_criteria cell_criteria;                ///< mesh cell constraints
    Mesh_criteria criteria;                     ///< mesh generation constraints
    
    C3t3 c3t3;                                  ///< generated mesh
    FEMesh3* M = NULL;                          ///< solver for systems on mesh
    
    /// mesh generation step
    void gen_mesh();
    /// mesh refinement step
    void refine_mesh();
    /// field solver step
    void setup_solver();
};

#endif
