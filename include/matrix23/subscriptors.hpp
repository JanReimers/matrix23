// File: Subscriptors.hpp  Define dense packing and indexing for full,tri,tridiagona,diagonal matrices.
#pragma once

#include <cstddef> //To get size_t
#include <ranges> //To get iota_view.
#include <cassert>

namespace matrix23
{
typedef std::ranges::iota_view<size_t,size_t> iota_view;

enum class Packing  {full, utri, ltri, tridiag, diag, sband, uband, lband};
enum class Indexing {row_major, col_major};
enum class Symmetry {none,symmetric, anit_symmetric, hermitian, anit_hermitian};

template <class P> concept isPacker  = requires (P const p,size_t i, size_t j, bool b)
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
template <class S> concept isSubscriptor = requires (S const s,size_t i, size_t j, bool b)
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
typedef std::function<size_t(size_t, size_t)> indexer_t;

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

class PackerCommon
{
public:
    PackerCommon(const size_t& _nrows, const size_t& _ncols, const indexer_t& ind) : nrows(_nrows), ncols(_ncols), indexer(ind) {};
    size_t nr() const {return nrows;}
    size_t nc() const {return ncols;}
    void range_check(size_t i, size_t j) const
    {
        assert(i<nrows && "   Row index ot of bounds");
        assert(j<ncols && "Column index ot of bounds");
    }
    size_t offset(size_t i, size_t j) const
    {
        range_check(i,j);
        return indexer(i,j);
    }

protected:
    const size_t& nrows,ncols;
    indexer_t indexer; // Function to calculate the index based on row and column

};
// If indices are (row,col) the row major means the linear data is stored in this order:
// [  (0,0),(0,1),(0,2)...(0,nc-1),(1,0), (1,1)...(1,nc-1) .....(nr-1,nc-1) ]

class FullPacker : public PackerCommon
{
public:
    FullPacker(const size_t& _nrows, const size_t& _ncols, Indexing ind) 
        : PackerCommon(_nrows,_ncols, make_indexer(ind,nrows,ncols)) {};
    bool is_stored(size_t i, size_t j) const {range_check(i,j);return true;} // Full matrix, all elements are stored
    size_t stored_size() const {return nrows * ncols;} // Total number of elements
    size_t stored_row_size(size_t) const {return ncols;}// Each row has nc elements
    size_t stored_col_size(size_t) const {return nrows;} // Each col has nr elements
private:
    static indexer_t make_indexer(Indexing ind, const size_t& nrows, const size_t& ncols)
    {
        switch (ind)
        {
            case Indexing::row_major:
                return [&ncols](size_t i, size_t j) -> size_t {return j + i*ncols;};
            case Indexing::col_major:
                return [&nrows](size_t i, size_t j) -> size_t {return i + j*nrows;};
        
        }
    }
};
class UpperTriangularPacker : public PackerCommon
{
public:
    UpperTriangularPacker(const size_t& _nrows, const size_t& _ncols, Indexing ind) 
        : PackerCommon(_nrows,_ncols,make_indexer(ind,nrows,ncols)) {};
    bool is_stored(size_t i, size_t j) const {return i <= j;}
    size_t stored_size() const
    {
        return nrows*(nrows+1)/2 + (ncols>nrows ? (ncols-nrows)*nrows : 0); // Total number of elements
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
private:
    static indexer_t make_indexer(Indexing ind, const size_t& nrows, const size_t& ncols)
    {
        switch (ind)
        {
            case Indexing::row_major:
                return [&ncols](size_t i, size_t j) -> size_t {return j + i*(2*ncols-i-1)/2;};
            case Indexing::col_major:
                return [      ](size_t i, size_t j) -> size_t {return i + j*(        j+1)/2;};
        
        }
    }

};
class LowerTriangularPacker : public PackerCommon
{
public:
    LowerTriangularPacker(const size_t& _nrows, const size_t& _ncols, Indexing ind) 
        : PackerCommon(_nrows,_ncols,make_indexer(ind,nrows,ncols)) {};
    bool is_stored(size_t i, size_t j) const {return j <= i;}
    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        return j + (i+1)*i/2;    // row major
        //return i + (j)*(2*nrows-j-1)/2; // col major   
    }
    size_t stored_size() const
    {
        return ncols*(ncols+1)/2  + (nrows>ncols ?  (nrows-ncols)*ncols :  0); // Total number of elements
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
private:
    static indexer_t make_indexer(Indexing ind, const size_t& nrows, const size_t& ncols)
    {
        switch (ind)
        {
            case Indexing::row_major:
                return [      ](size_t i, size_t j) -> size_t {return j + i*(        i+1)/2;};
            case Indexing::col_major:
                return [&nrows](size_t i, size_t j) -> size_t {return i + j*(2*nrows-j-1)/2;};
        
        }
    }
  
 
};
class DiagonalPacker : public PackerCommon
{
public:
    DiagonalPacker(const size_t& _nrows, const size_t& _ncols, Indexing ind) 
        : PackerCommon(_nrows,_ncols, [](size_t i, size_t) -> size_t {return i;}) {};
    bool is_stored(size_t i, size_t j) const {return i == j;}
    size_t stored_size() const {return std::min(nrows,ncols);}
    size_t stored_row_size(size_t) const {return 1;}
    size_t stored_col_size(size_t) const {return 1;}
};
//
//      *    *   a02  a13  a24  a35  
//      *   a01  a12  a23  a34  a45 
//     a00  a11  a22  a33  a44  a55 
//     a10  a21  a32  a43  a54  *   
//     a20  a31  a42  a53  *    *
class SBandPacker : public PackerCommon
{
public:
    //  Handle square only.
    SBandPacker(const size_t& n, size_t _k, Indexing ind) //k=band width,
        : PackerCommon(n,n, [&n,_k](size_t i, size_t j) -> size_t {return _k+i-j+j*n;}), k(_k){};
    bool is_stored(size_t i, size_t j) const {range_check(i,j);return i+k<=j && j+k<=i;}
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
private:
    size_t k;
};


static_assert(isPacker<FullPacker>);
static_assert(isPacker<UpperTriangularPacker>);
static_assert(isPacker<LowerTriangularPacker>);
static_assert(isPacker<DiagonalPacker>);
static_assert(isPacker<SBandPacker>);

} // namespace matrix23