/// \file BrianFields.hh Loader for gridded electric field data

#ifndef BRIANFIELDS_HH
#define BRIANFIELDS_HH

#include "NCubicGrid.hh"

typedef NCubicGrid<3,float> TriCubic;

void load_Brian_mirror(TriCubic G[3], bool reload = false);

#endif
