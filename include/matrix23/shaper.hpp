// File shaper.hpp   Define the shaper concept and related classes.
#pragma once

#include <cstddef> //To get size_t
#include <ranges>  //To get iota_view.

//
// A matrix shape in concerned with what elements are non-zero. But a shape is unconcerned
// with if or how these elements are actually stored in memory. For example an upper triangular
// matrix has all zeros below the diagonal.  It can implemented using efficient 
// packed N*(N+1)/2) storage (where the zeros are not stored) or full N*N storage (where all 
// the zeros are stored).  In both cases the *shape* is upper triangular.
//
namespace matrix23
{
typedef std::ranges::iota_view<size_t,size_t> iota_view; //Used for specifying index ranges.

// c++20 concept definition for a shaper.
template <class S> concept isShaper = requires (S const s,size_t i, size_t j, bool b)
{
    s.nonzero_row_indexes(j); //range of row indices in column j that are non zero.
    s.nonzero_col_indexes(i); //range of column indices in row i that are non zero.
    s.nr(); // Total number of rows.
    s.nc(); // Total number of columns.
};

// Common base class for all shapers.
class ShaperCommon
{
public:
    ShaperCommon(size_t nr, size_t nc) : nrows(nr), ncols(nc){};
    size_t nr() const {return nrows;}
    size_t nc() const {return ncols;}
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
    SBandShaper(size_t nr, size_t nc, size_t _k) : ShaperCommon(nr,nc), k(_k) {};
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
    size_t bandwidth() const {return k;}
    size_t k;
};

static_assert(isShaper<           FullShaper>);
static_assert(isShaper<UpperTriangularShaper>);
static_assert(isShaper<LowerTriangularShaper>);
static_assert(isShaper<       DiagonalShaper>);
static_assert(isShaper<          SBandShaper>);

}; //namespace matrix23



