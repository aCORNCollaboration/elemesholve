/// \file CalculationProcess.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#include "CalculationProcess.hh"
#include <stdio.h>

CalculationProcess::CalculationProcess():
YSQ(0, 10., 0.7),
myCT(NULL),
EMG(myCT),
STG(myCT),
myGeom(&EMG),
GW(myGeom->theWorld),
RadiusMesh(myGeom->theWorld, 0.5),
domain(GW, myGeom->theWorld->myBounds, 1e-6),
edge_criterea(RadiusMesh),
facet_criteria(29,           // angle bound
               RadiusMesh,   // radius bound field
               2e-3,         // distance bound
               CGAL::FACET_VERTICES_ON_SURFACE ),     // facet topology requirement
               //CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK ),
cell_criteria(3,             // radius-edge ratio
              RadiusMesh),   // sizing field
criteria(edge_criterea, facet_criteria, cell_criteria)
{
    YSQ.z0 = EMG.myWorld.WC.gridz;
    myGeom->add_features(domain);
}

void CalculationProcess::gen_mesh() {
    printf("Generating mesh...\n");
    c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());
}

void CalculationProcess::refine_mesh() {
    printf("Perturbing mesh...\n");
    CGAL::perturb_mesh_3(c3t3, domain, time_limit = 10);
    printf("Exuding mesh...\n");
    CGAL::exude_mesh_3(c3t3, time_limit = 10);
}

void CalculationProcess::setup_solver() {
    printf("Setting up solver...\n");
    myGeom->calc_bvals(c3t3);
    M = new FEMesh3(c3t3, myCT);
    M->set_boundary_points(*myGeom);
    M->set_boundary_values(*myGeom);
}

