// File: Subscriptors.hpp  Define dense packing and indexing for full,tri,tridiagona,diagonal matrices.
#pragma once

#include <cstddef> //To get size_t
#include <ranges> //To get iota_view.
#include <cassert>

namespace matrix23
{
typedef std::ranges::iota_view<size_t,size_t> iota_view;

// If indices are (row,col) then:
// row major means the linear data is stored in this order:
//    [  (0,0),(0,1),(0,2)...(0,nc-1),(1,0), (1,1)...(1,nc-1) .....(nr-1,nc-1) ]
// col major means the linear data is stored in this order:
//    [  (0,0),(1,0),(2,0)...(nr-1,0), (0,1),(1,1)...(nr-1,1) .....(nr-1,nc-1) ]
//    col major is how Fortran, Lapack and Blas store a matrix.

enum class Indexing {row_major, col_major};

template <class P> concept isPacker = requires (P const p,size_t i, size_t j, bool b)
{
    b=p.is_stored(i,j);
    i=p.offset(i,j);
    i=p.stored_size();
    i=p.stored_row_size(j);
    i=p.stored_col_size(j);
};
template <class S> concept isShaper = requires (S const s,size_t i, size_t j, bool b)
{
    s.nonzero_row_indexes(j);
    s.nonzero_col_indexes(i);  
};

#include <ranges>
class PackerCommon
{
public:
    PackerCommon(const size_t& _nrows, const size_t& _ncols) : nrows(_nrows), ncols(_ncols) {};
    size_t nr() const {return nrows;}
    size_t nc() const {return ncols;}
    void range_check(size_t i, size_t j) const
    {
        assert(i<nrows && "   Row index ot of bounds");
        assert(j<ncols && "Column index ot of bounds");
    }
protected:
    const size_t nrows,ncols;
};

class FullPacker           : public PackerCommon
{
public:
    using PackerCommon::PackerCommon; // Inherit constructors
    bool is_stored(size_t i, size_t j) const {range_check(i,j);return true;} // Full matrix, all elements are stored
    size_t stored_size() const {return nrows * ncols;} // Total number of elements
    size_t stored_row_size(size_t row) const {assert(row<nrows);return ncols;}// Each row has nc elements
    size_t stored_col_size(size_t col) const {assert(col<ncols);return nrows;}// Each col has nr elements
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
    size_t stored_row_size(size_t row_index) const
     {
        assert(row_index < nrows);
        if (row_index >= ncols) return 0; // No elements stored in this row
        return ncols-row_index; // 
    }
    size_t stored_col_size(size_t col_index) const
    {
        assert(col_index < ncols);
        if (col_index >= nrows) return nrows; // No elements stored in this row
        return col_index+1; // 
    }
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
    size_t stored_row_size(size_t row_index) const
    {
        assert(row_index < nrows);
        if (row_index >= ncols) return ncols; // full row below the triangle
        return row_index+1; // in the triangle
    }
    size_t stored_col_size(size_t col_index) const
    {
        assert(col_index < ncols);
        if (col_index >= nrows) return 0; // full row below the triangle
        return nrows-col_index; // in the triangle
    }
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
};

class DiagonalPacker        : public PackerCommon
{
public:
    using PackerCommon::PackerCommon; // Inherit constructors
    bool is_stored(size_t i, size_t j) const {return i == j;}
    size_t stored_size() const {return std::min(nrows,ncols);}
    size_t stored_row_size(size_t row) const {return row<ncols ? 1 : 0;}
    size_t stored_col_size(size_t col) const {return col<nrows ? 1 : 0;}
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return i;
    }
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
    size_t stored_row_size(size_t row) const 
    {
        if (row<k) return row+1+k;
        if (row>nrows-k-1) return nrows-row+k;
        return 2*k+1;  
    }
    size_t stored_col_size(size_t col) const //Don't need the col index.
    {
        if (col<k) return col+1+k;
        if (col>ncols-k-1) return ncols-col+k;
        return 2*k+1; 
    }
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

static_assert(isPacker<           FullPackerCM>);
static_assert(isPacker<           FullPackerRM>);
static_assert(isPacker<UpperTriangularPackerCM>);
static_assert(isPacker<UpperTriangularPackerRM>);
static_assert(isPacker<LowerTriangularPackerCM>);
static_assert(isPacker<LowerTriangularPackerRM>);
static_assert(isPacker<       DiagonalPacker  >);
static_assert(isPacker<          SBandPacker  >);

//--------------------------------------------------------------------------------------
//
//  Shapers. These decide which part/range of a row or column is non-zero; *regardless* of whether or not
//  thoses zeros are actually stored.  For example an upper triangular matrix can be will full packing/storage.
//
class ShaperCommon
{
public:
    ShaperCommon(const PackerCommon& p) : nrows(p.nrows), ncols(p.ncols){};
protected:
    const size_t nrows,ncols;
};
class FullShaper            : public ShaperCommon
{
public:
    using ShaperCommon::ShaperCommon; // Inherit constructors
    iota_view nonzero_row_indexes(size_t col) const {return std::views::iota(size_t(0),nrows);}
    iota_view nonzero_col_indexes(size_t row) const {return std::views::iota(size_t(0),ncols);}   
};
class UpperTriangularShaper : public ShaperCommon
{
public:
    using ShaperCommon::ShaperCommon; // Inherit constructors
    iota_view nonzero_row_indexes(size_t col) const {return std::views::iota(size_t(0),std::min(col+1,nrows));}
    iota_view nonzero_col_indexes(size_t row) const {return std::views::iota(row      ,ncols);} 
};
class LowerTriangularShaper : public ShaperCommon
{
public:
    using ShaperCommon::ShaperCommon; // Inherit constructors
    iota_view nonzero_row_indexes(size_t col) const {return std::views::iota(col      ,nrows);}
    iota_view nonzero_col_indexes(size_t row) const {return std::views::iota(size_t(0),std::min(row+1,ncols));}  
};
class DiagonalShaper        : public ShaperCommon
{
public:
    using ShaperCommon::ShaperCommon; // Inherit constructors
    iota_view nonzero_row_indexes(size_t col) const {return std::views::iota(col,std::min(col+1,nrows));}
    iota_view nonzero_col_indexes(size_t row) const {return std::views::iota(row,std::min(row+1,ncols));}  
};
class SBandShaper           : public ShaperCommon
{
public:
    SBandShaper(const SBandPacker& p) : ShaperCommon(p), k(p.k) {};
    iota_view nonzero_row_indexes(size_t col) const 
    {
        size_t i0=col<k ? 0 : col-k;
        size_t i1=col+k>=nrows ? nrows : col+k+1;
        return std::views::iota(i0 ,i1);
    }
    iota_view nonzero_col_indexes(size_t row) const 
    {
        size_t i0=row<k ? 0 : row-k;
        size_t i1=row+k>=ncols ? ncols : row+k+1;
        return std::views::iota(i0 ,i1);
    }    
    size_t k;
};

static_assert(isShaper<           FullShaper>);
static_assert(isShaper<UpperTriangularShaper>);
static_assert(isShaper<LowerTriangularShaper>);
static_assert(isShaper<       DiagonalShaper>);
static_assert(isShaper<          SBandShaper>);

//
//  Define what matrix shapes result from multiply two matricies.
//  These rules are identical the Packers and Shapers not sure if there is a clever
//  way to define the rules once for both Packer and Shaper?

template <isPacker P1, isPacker P2> struct MatrixProductPackerType;

template <isPacker P> struct MatrixProductPackerType<FullPackerCM,P> {typedef FullPackerCM packer_t;};
template <isPacker P> struct MatrixProductPackerType<P,FullPackerCM> {typedef FullPackerCM packer_t;};
template <isPacker P> struct MatrixProductPackerType<FullPackerRM,P> {typedef FullPackerRM packer_t;};
template <isPacker P> struct MatrixProductPackerType<P,FullPackerRM> {typedef FullPackerRM packer_t;};
template <isPacker P> struct MatrixProductPackerType<P,DiagonalPacker> {typedef P packer_t;};
template <isPacker P> struct MatrixProductPackerType<DiagonalPacker,P> {typedef P packer_t;};
template <> struct MatrixProductPackerType<FullPackerCM,FullPackerCM> {typedef FullPackerCM packer_t;};
template <> struct MatrixProductPackerType<FullPackerRM,FullPackerRM> {typedef FullPackerRM packer_t;};
template <> struct MatrixProductPackerType<DiagonalPacker,FullPackerCM> {typedef FullPackerCM packer_t;};
template <> struct MatrixProductPackerType<FullPackerCM,DiagonalPacker> {typedef FullPackerCM packer_t;};
template <> struct MatrixProductPackerType<DiagonalPacker,FullPackerRM> {typedef FullPackerRM packer_t;};
template <> struct MatrixProductPackerType<FullPackerRM,DiagonalPacker> {typedef FullPackerRM packer_t;};
template <> struct MatrixProductPackerType<DiagonalPacker,DiagonalPacker> {typedef DiagonalPacker packer_t;};
template <> struct MatrixProductPackerType<UpperTriangularPackerCM,UpperTriangularPackerCM> {typedef UpperTriangularPackerCM packer_t;};
template <> struct MatrixProductPackerType<LowerTriangularPackerCM,LowerTriangularPackerCM> {typedef LowerTriangularPackerCM packer_t;};
template <> struct MatrixProductPackerType<UpperTriangularPackerCM,LowerTriangularPackerCM> {typedef FullPackerCM packer_t;};
template <> struct MatrixProductPackerType<LowerTriangularPackerCM,UpperTriangularPackerCM> {typedef FullPackerCM packer_t;};
template <> struct MatrixProductPackerType<SBandPacker,SBandPacker> {typedef SBandPacker packer_t;}; //Need to add the ks somehow.



template <isShaper P1, isShaper P2> struct MatrixProductShaperType;
template <isShaper S> struct MatrixProductShaperType<FullShaper,S> {typedef FullShaper shaper_t;};
template <isShaper S> struct MatrixProductShaperType<S,FullShaper> {typedef FullShaper shaper_t;};
template <isShaper S> struct MatrixProductShaperType<S,DiagonalShaper> {typedef S shaper_t;};
template <isShaper S> struct MatrixProductShaperType<DiagonalShaper,S> {typedef S shaper_t;};
template <> struct MatrixProductShaperType<FullShaper,FullShaper> {typedef FullShaper shaper_t;};
template <> struct MatrixProductShaperType<DiagonalShaper,FullShaper> {typedef FullShaper shaper_t;};
template <> struct MatrixProductShaperType<FullShaper,DiagonalShaper> {typedef FullShaper shaper_t;};
template <> struct MatrixProductShaperType<DiagonalShaper,DiagonalShaper> {typedef DiagonalShaper shaper_t;};
template <> struct MatrixProductShaperType<UpperTriangularShaper,UpperTriangularShaper> {typedef UpperTriangularShaper shaper_t;};
template <> struct MatrixProductShaperType<LowerTriangularShaper,LowerTriangularShaper> {typedef LowerTriangularShaper shaper_t;};
template <> struct MatrixProductShaperType<UpperTriangularShaper,LowerTriangularShaper> {typedef FullShaper shaper_t;};
template <> struct MatrixProductShaperType<LowerTriangularShaper,UpperTriangularShaper> {typedef FullShaper shaper_t;};
template <> struct MatrixProductShaperType<SBandShaper,SBandShaper> {typedef SBandShaper shaper_t;}; //Need to add the ks somehow.


} // namespace matrix23