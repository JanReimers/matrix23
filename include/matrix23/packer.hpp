// File: packer.hpp  Define the packer concept and related classes.
#pragma once

#include "matrix23/shaper.hpp"
#include <cassert>

//
// Packers define how matrix elements are arranged in memory.  For most packings (except diagonal)
// column major and row major arrangements are possible. In addition triangular and banded
// matrices can also be efficiently packed so that no memory is not wasted on definitely zero elements.
// If indices are (row,col) then row major and column major layouts look like this: :
// row major means the linear data is stored in this order:
//    [  (0,0),(0,1),(0,2)...(0,nc-1),(1,0), (1,1)...(1,nc-1) .....(nr-1,nc-1) ]
// col major means the linear data is stored in this order:
//    [  (0,0),(1,0),(2,0)...(nr-1,0), (0,1),(1,1)...(nr-1,1) .....(nr-1,nc-1) ]
//    col major is how Fortran, Lapack and Blas store a matrix.
//
namespace matrix23
{
// c++20 concept definition for a shaper.
template <class P> concept isPacker = requires (P const p,P pnc,size_t i, size_t j, bool b)
{
    pnc.resize(i,j); //nr nc
    b=p.is_stored(i,j); 
    i=p.offset(i,j);  // linear 1D offset for a 2D index pair.
    i=p.stored_size();
    p.shaper();
    p.transpose();
};

// Common base class for all packers.
class PackerCommon
{
public:
    PackerCommon(const size_t& _nrows, const size_t& _ncols) : nrows(_nrows), ncols(_ncols) {};
    void resize(size_t nr, size_t nc) {nrows=nr;ncols=nc;}
    size_t nr() const {return nrows;}
    size_t nc() const {return ncols;}
    void range_check(size_t i, size_t j) const
    {
        assert(i<nrows && "   Row index ot of bounds");
        assert(j<ncols && "Column index ot of bounds");
    }
protected:
    size_t nrows,ncols;
};

class FullPacker           : public PackerCommon
{
public:
    using PackerCommon::PackerCommon; // Inherit constructors
    bool is_stored(size_t i, size_t j) const {range_check(i,j);return true;} // Full matrix, all elements are stored
    size_t stored_size() const {return nrows * ncols;} // Total number of elements
    FullShaper shaper() const {return FullShaper(nr(),nc());}
};
class FullPackerCM         : public FullPacker
{
public:
    using FullPacker::FullPacker; // Inherit constructors
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return i + j*nrows;
    }
    auto transpose() const {return FullPackerCM(nc(),nr());}
};
class FullPackerRM         : public FullPacker
{
public:
    using FullPacker::FullPacker; // Inherit constructors
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return j + i*ncols;
    }
    auto transpose() const {return FullPackerRM(nc(),nr());}
};

class UpperTriangularPacker : public PackerCommon
{
public:
    using PackerCommon::PackerCommon; // Inherit constructors
    bool is_stored(size_t i, size_t j) const {return i <= j;}
    size_t stored_size() const
    {
        size_t n=std::min(nrows,ncols);
        return n*(n+1)/2 + (ncols>nrows ? (ncols-nrows)*nrows : 0); // Total number of elements
    }
    UpperTriangularShaper shaper() const {return UpperTriangularShaper(nr(),nc());}
};
class UpperTriangularPackerRM : public UpperTriangularPacker
{
public:
    using UpperTriangularPacker::UpperTriangularPacker; // Inherit constructors
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return j + i*(2*ncols-i-1)/2;
    }
    auto transpose() const;
};
class UpperTriangularPackerCM : public UpperTriangularPacker
{
public:
    using UpperTriangularPacker::UpperTriangularPacker; // Inherit constructors
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return i + j*(        j+1)/2;
    }
    auto transpose() const;
};

class LowerTriangularPacker : public PackerCommon
{
public:
    using PackerCommon::PackerCommon; // Inherit constructors
    bool is_stored(size_t i, size_t j) const {return j <= i;}
    size_t stored_size() const
    {
        size_t n=std::min(nrows,ncols);
        return n*(n+1)/2  + (nrows>ncols ?  (nrows-ncols)*ncols :  0); // Total number of elements
    }
    LowerTriangularShaper shaper() const {return LowerTriangularShaper(nr(),nc());}
};
class LowerTriangularPackerRM : public LowerTriangularPacker
{
public:
    using LowerTriangularPacker::LowerTriangularPacker; // Inherit constructors
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return j + i*(        i+1)/2;
    }
    auto transpose() const;
};
class LowerTriangularPackerCM : public LowerTriangularPacker
{
public:
    using LowerTriangularPacker::LowerTriangularPacker; // Inherit constructors
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return i + j*(2*nrows-j-1)/2;
    }
    auto transpose() const;
};

inline auto UpperTriangularPackerRM::transpose() const {return LowerTriangularPackerRM(nc(),nr());}
inline auto UpperTriangularPackerCM::transpose() const {return LowerTriangularPackerCM(nc(),nr());}
inline auto LowerTriangularPackerRM::transpose() const {return UpperTriangularPackerRM(nc(),nr());}
inline auto LowerTriangularPackerCM::transpose() const {return UpperTriangularPackerCM(nc(),nr());}

class DiagonalPacker        : public PackerCommon
{
public:
    using PackerCommon::PackerCommon; // Inherit constructors
    bool is_stored(size_t i, size_t j) const {return i == j;}
    size_t stored_size() const {return std::min(nrows,ncols);}
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return i;
    }
    DiagonalShaper shaper() const {return DiagonalShaper(nr(),nc());}
    auto transpose() const {return DiagonalPacker(nc(),nr());}
};
class SBandPacker           : public PackerCommon
{
// Packing guide:
//      *    *   a02  a13  a24  a35  
//      *   a01  a12  a23  a34  a45 
//     a00  a11  a22  a33  a44  a55 
//     a10  a21  a32  a43  a54  *   
//     a20  a31  a42  a53  *    *
public:
    //  Handle square only.
    SBandPacker(const size_t& n, size_t _k) : PackerCommon(n,n) , k(_k){};
    bool is_stored(size_t i, size_t j) const {range_check(i,j);return j<=i+k && i<=j+k;}
    size_t stored_size() const {return nrows * (2*k+1);}
    SBandShaper shaper() const {return SBandShaper(nr(),nc(),k);}
    auto transpose() const {return SBandShaper(nc(),nr(),k);}
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return k+i-j+j*(2*k+1);
    }
    size_t bandwidth() const {return k;}
    private:
    friend class SBandShaper;
    const size_t k;
};


} // namespace
