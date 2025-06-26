// File: Subscriptors.hpp  Define dense packing and indexing for full,tri,tridiagona,diagonal matrices.
#pragma once

#include <cstddef> //To get size_t
#include <ranges> //To get iota_view.
#include <cassert>

namespace matrix23
{
typedef std::ranges::iota_view<size_t,size_t> iota_view;

template <class S> 
concept isSubscriptor = requires (S const s,size_t i, size_t j, bool b)
{
    b=s.is_stored(i, j);
    i=s.offset(i,j);
    i=s.size();
    i=s.stored_row_size(j);
    s.nonzero_row_indexes(j);
    s.nonzero_col_indexes(i);
    i=s.nr();
    i=s.nc();

};

class SubsciptorCommon
{
public:
    SubsciptorCommon(size_t _nrows, size_t _ncols) : nrows(_nrows), ncols(_ncols) {};
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
// If indices are (row,col) the row major means the linear data is stored in this order:
// [  (0,0),(0,1),(0,2)...(0,nc-1),(1,0), (1,1)...(1,nc-1) .....(nr-1,nc-1) ]
class FullRowMajorSubsciptor : public SubsciptorCommon
{
public:
    using SubsciptorCommon::SubsciptorCommon; // Inherit constructors.
    bool is_stored(size_t i, size_t j) const
    {
        range_check(i,j);
        return true; // Full matrix, all elements are stored
    }
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return i * ncols + j;
    }
    size_t size() const
    {
        return nrows * ncols; // Total number of elements
    }
    size_t stored_row_size(size_t) const //Don't need the row index.
    {
        return ncols; // Each row has nc elements
    }
    size_t stored_col_size(size_t) const //Don't need the col index.
    {
        return nrows; // Each col has nr elements
    }
    iota_view nonzero_row_indexes(size_t col) const
    {
        return std::views::iota(size_t(0),nrows);
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        return std::views::iota(size_t(0),ncols);
    }
    
};

// If indices are (row,col) the col major means the linear data is stored in this order:
// [  (0,0),(1,0),(2,0)...(nr-1,0), (0,1),(1,1)...(nr-1,1) .....(nr-1,nc-1) ]
// This is how Fortran, Lapack and Blas store a matrix.
class FullColMajorSubsciptor  : public SubsciptorCommon
{
public:
    using SubsciptorCommon::SubsciptorCommon; // Inherit constructors.
    bool is_stored(size_t i, size_t j) const
    {
        range_check(i,j);
        return true; // Full matrix, all elements are stored
    }
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return i + j*nrows;
    }
    size_t size() const
    {
        return nrows * ncols; // Total number of elements
    }
    size_t stored_row_size(size_t) const //Don't need the row index.
    {
        return ncols; // Each row has ncols elements
    }
    size_t stored_col_size(size_t) const //Don't need the row index.
    {
        return nrows; // Each col has nr elements
    }
    iota_view nonzero_row_indexes(size_t col) const
    {
        return std::views::iota(size_t(0),nrows);
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        return std::views::iota(size_t(0),ncols);
    }
};

class UpperTriangularRowMajorSubsciptor  : public SubsciptorCommon
{
public:
    using SubsciptorCommon::SubsciptorCommon; // Inherit constructors.
    bool is_stored(size_t i, size_t j) const
    {
        return i <= j; 
    }
    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        return i * (2*ncols-i-1) / 2 + j; // Upper triangular matrix        
    }
    size_t size() const
    {
        return ncols<=nrows ? ncols*(ncols+1)/2 : nrows*(nrows+1)/2 +(ncols-nrows)*nrows; // Total number of elements
    }
    size_t stored_row_size(size_t row_index) const
    {
        assert(row_index < nrows);
        if (row_index >= ncols) return 0; // No elements stored in this row
        return ncols-row_index + ((ncols>nrows) ? ncols-nrows : 0); // Each row has ncols elements
    }
    template <std::ranges::range R> auto view(R& r) const
    {
        auto row = [r,this](int i) mutable 
        {
            return r | std::views::drop(offset(i,i)) | std::views::take(stored_row_size(i));
        };
        return std::views::iota(0,(int)nrows) | std::views::transform(row);
    }
    iota_view nonzero_row_indexes(size_t col) const
    {
        return std::views::iota(size_t(0),col+1);
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        return std::views::iota(row,ncols);
    }
    
};

class UpperTriangularColMajorSubsciptor  : public SubsciptorCommon
{
public:
    using SubsciptorCommon::SubsciptorCommon; // Inherit constructors.
    bool is_stored(size_t i, size_t j) const
    {
        return i <= j; 
    }
    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        return i + (j+1)*j/2; // Upper triangular matrix        
    }
    size_t size() const
    {
        return ncols<=nrows ? ncols*(ncols+1)/2 : nrows*(nrows+1)/2 +(ncols-nrows)*nrows; // Total number of elements
    }
    size_t stored_row_size(size_t row_index) const
    {
        assert(row_index < nrows);
        if (row_index >= ncols) return 0; // No elements stored in this row
        return ncols-row_index + ((ncols>nrows) ? ncols-nrows : 0); // Each row has ncols elements
    }
    template <std::ranges::range R> auto view(R& r) const
    {
        auto row = [r,this](int i) mutable 
        {
            return r | std::views::drop(offset(i,i)) | std::views::take(stored_row_size(i));
        };
        return std::views::iota(0,(int)nrows) | std::views::transform(row);
    }
    iota_view nonzero_row_indexes(size_t col) const
    {
        return std::views::iota(size_t(0),col+1);
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        return std::views::iota(row,ncols);
    }
    
};
class LowerTriangularColMajorSubsciptor  : public SubsciptorCommon
{
public:
    using SubsciptorCommon::SubsciptorCommon; // Inherit constructors.
    bool is_stored(size_t i, size_t j) const
    {
        return j <= i; 
    }
    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        return i + (j)*(2*nrows-j-1)/2; // Upper triangular matrix        
    }
    size_t size() const
    {
        return ncols*(ncols+1)/2  + (ncols<nrows ?  (nrows-ncols)*ncols :  0); // Total number of elements
    }
    size_t stored_row_size(size_t row_index) const
    {
        assert(row_index < nrows);
        if (row_index >= ncols) return ncols; // full row below the triangle
        return row_index+1; // in the triangle
    }
    template <std::ranges::range R> auto view(R& r) const
    {
        auto row = [r,this](int i) mutable 
        {
            return r | std::views::drop(offset(i,0)) | std::views::take(stored_row_size(i));
        };
        return std::views::iota(0,(int)nrows) | std::views::transform(row);
    }
    iota_view nonzero_row_indexes(size_t col) const
    {
        return std::views::iota(col,nrows); // All rows from col+1 to nrows
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        return std::views::iota(size_t(0),row+1);
    }
    
};
class LowerTriangularRowMajorSubsciptor  : public SubsciptorCommon
{
public:
    using SubsciptorCommon::SubsciptorCommon; // Inherit constructors.
    bool is_stored(size_t i, size_t j) const
    {
        return j <= i; 
    }
    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        return j + (i+1)*i/2;    
    }
    size_t size() const
    {
        return ncols*(ncols+1)/2  + (ncols<nrows ?  (nrows-ncols)*ncols :  0); // Total number of elements
    }
    size_t stored_row_size(size_t row_index) const
    {
        assert(row_index < nrows);
        if (row_index >= ncols) return ncols; // full row below the triangle
        return row_index+1; // in the triangle
    }
    template <std::ranges::range R> auto view(R& r) const
    {
        auto row = [r,this](int i) mutable 
        {
            return r | std::views::drop(offset(i,0)) | std::views::take(stored_row_size(i));
        };
        return std::views::iota(0,(int)nrows) | std::views::transform(row);
    }
    iota_view nonzero_row_indexes(size_t col) const
    {
        return std::views::iota(col,nrows); // All rows from col+1 to nrows
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        return std::views::iota(size_t(0),row+1);
    }
    
};

class DiagonalSubsciptor  : public SubsciptorCommon
{
    public:
    using SubsciptorCommon::SubsciptorCommon; // Inherit constructors.
    bool is_stored(size_t i, size_t j) const
    {
        return i == j; 
    }
    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        return i; // Diagonal matrix
    }
    size_t size() const
    {
        return std::min(nrows,ncols); // Total number of elements
    }
    size_t stored_row_size(size_t) const
    {
        return 1;
    }
    iota_view nonzero_row_indexes(size_t col) const
    {
        return std::views::iota(col,col+1);
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        return std::views::iota(row,row+1);
    }
   
};

class TriDiagonalSubsciptor  : public SubsciptorCommon
{
    public:
    using SubsciptorCommon::SubsciptorCommon; // Inherit constructors.
    bool is_stored(size_t i, size_t j) const
    {
        return i<=j+1 && j<=i+1; // Main diagonal and two adjacent diagonals
    }
    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        return 2*i+j; // Tri-diagonal matrix
    }
    size_t size() const
    {
        assert(nrows==ncols);
        
        return nrows>0 ? 3*nrows-2 : 0; // Total number of elements
    }
    size_t stored_row_size(size_t row_index) const
    {
        return (row_index==0 || row_index==nrows-1) ? 2 : 3;
    }
    
    iota_view nonzero_row_indexes(size_t col) const
    {
        size_t c0=std::max(size_t(1),col)-1, c1=std::min(nrows,col+2);
        return std::views::iota(c0,c1);
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        size_t r0=std::max(size_t(1),row)-1,r1=std::min(ncols,row+2);
        return std::views::iota(r0,r1);
    }
   
};

class BandedSubsciptor
{
    public:
    BandedSubsciptor(size_t nrows, size_t ncols, size_t _k)
        :  nr(nrows), nc(ncols), k(_k) {}

    bool is_stored(size_t i, size_t j) const
    {
        return (i<=j+k && j<=i+k);
    }

    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        size_t dj=j;
        if (i>k) dj+=i;
        return (i * ( k + 1)) + dj; // Banded matrix
    }
    size_t size() const
    {
        return (nr * (2*k + 1)) - k * (k + 1) / 2; // Total number of elements
    
    }
    size_t stored_row_size(size_t row_index) const
    {
        if (row_index < k) return row_index + 1; // First k rows have increasing size
        if (row_index >= nr - k) return nc - row_index + k; // Last k rows have decreasing size
        return 2 * k + 1; // Middle rows have full band width
    }
    
    size_t nr,nc,k;

};

static_assert(isSubscriptor<FullRowMajorSubsciptor>);
static_assert(isSubscriptor<FullColMajorSubsciptor>);
static_assert(isSubscriptor<UpperTriangularRowMajorSubsciptor>);
static_assert(isSubscriptor<DiagonalSubsciptor>);
static_assert(isSubscriptor<TriDiagonalSubsciptor>);
// static_assert(isSubscriptor<BandedSubsciptor>);

} // namespace matrix23