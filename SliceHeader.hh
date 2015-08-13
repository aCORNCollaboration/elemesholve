/// \file SliceHeader.hh Header info for mesh slice data dumps
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
// 
// -- Michael P. Mendenhall, 2015

#ifndef SLICEHEADER_HH
#define SLICEHEADER_HH

/// Header information for a (hyper)plane
template<size_t D, typename val_tp>
class SliceHeader {
public:
    /// Constructor
    SliceHeader() { }
    
    val_tp x[D];                ///< point in plane
    val_tp basis[D][D];         ///< projection coordinate basis vectors (last is plane normal) 
    
    /// load from file
    void read(istream& is) {
        is.read((char*)x, D*sizeof(x[0]));
        for(size_t i=0; i<D; i++) is.read((char*)basis[i], D*sizeof(basis[0][0]));
    }
    /// write to file
    void write(ostream& o) const {
        o.write((char*)x, D*sizeof(x[0]));
        for(size_t i=0; i<D; i++) o.write((char*)basis[i], D*sizeof(basis[0][0]));
    }
};

#endif
