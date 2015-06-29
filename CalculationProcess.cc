#include "CalculationProcess.hh"
#include <stdio.h>

CalculationProcess::CalculationProcess():
YSQ(0, 10., 0.7),
G(&YSQ),
GW(G.theWorld),
RadiusMesh(G.theWorld, 0.5),
domain(GW, G.theWorld->myBounds, 1e-6),
edge_criterea(RadiusMesh),
facet_criteria(15,           // angle bound
               RadiusMesh,   // radius bound field
               1e-2,         // distance bound
               CGAL::FACET_VERTICES_ON_SURFACE ),     // facet topology requirement
               //CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK ),
cell_criteria(5,              // radius-edge ratio
              RadiusMesh),      // sizing field
criteria(edge_criterea, facet_criteria, cell_criteria)
{
    YSQ.z0 = G.myWorld.WC.gridz;
    G.add_features(domain);
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
    G.calc_bvals(c3t3);
    M = new FEMesh3(c3t3, &YSQ);
    M->set_boundary_points(G);
    M->set_boundary_values(G);
}

