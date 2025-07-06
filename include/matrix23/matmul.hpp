// File: matmul.hpp  Infrastructure required for general matrix pultiplication.
#pragma once

#include "matrix23/matrix1.hpp"

//
//  The requirements we want to meet are:
//      1) Support all combinations shapes, packings and symmetries.
//      2) Support chained operations A*B*C*D without creation of temporaries.  (lazy eval)
//      3) Support mixed element types Matrix<double>*Matrix<std::complex<float>>.
//      4) Don't waste time on definite zeros for triangular, band and diagonal matrices
//      5) propagate packings and shapes corretly.  For example full=upper*lower, upper=upper*upper
//      6) Achieve the same performance as the hand coded reference copy/cache algorithm:
// M mymul_wcopy(const M& A, const M& B)
// {
//     M C(A.nr(),B.nc());
//     for (size_t i=0;i<A.nr();i++)
//     {
//         Vector<double> Ai=A.row(i); //Cache the row which is expesive to index for col major packing.
//         for (size_t j=0;j<B.nc();j++)
//         {
//             double t=0.0;
//             for (size_t k=0;k<A.nc();k++)
//                 t+=Ai(k)*B(k,j);
//             C(i,j)=t;
//         }
//     }
//     return C;
// }
//

namespace matrix23
{
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
//
//  Create product packers and shapers.
//
template <isPacker A, isPacker B> auto MatrixProductPacker(const A& a, const B& b)
{
    return typename MatrixProductPackerType<A,B>::packer_t(a.nr(),b.nc());
}
template <isShaper A, isShaper B> auto MatrixProductShaper(const A& a, const B& b)
{
    return typename MatrixProductShaperType<A,B>::shaper_t(a.nr(),b.nc());
}
template <> inline auto MatrixProductPacker(const SBandPacker& a, const SBandPacker& b)
{
    assert(a.nr()==b.nr());
    return SBandPacker(a.nr(),a.bandwidth()+b.bandwidth());
}
template <> inline auto MatrixProductShaper(const SBandShaper& a, const SBandShaper& b)
{
    assert(a.nr()==b.nr());
    return SBandShaper(a.nr(),b.nc(),a.bandwidth()+b.bandwidth());
}

//
//  Lazy evaluated view of a matrix product.
//
template <std::ranges::viewable_range R, std::ranges::viewable_range C, isPacker P, isShaper S> class MatrixProductView
{
public:
    typedef std::ranges::range_value_t<R> Rv; //rows value type which is a column range.
    typedef std::ranges::range_value_t<Rv> value_type; //column range value type which should be scalar (double etc.)
    MatrixProductView(const R& _rows, const C& _cols,P _packer, S _shaper )
    : a_rows(_rows), b_cols(_cols), itsPacker(_packer), itsShaper(_shaper)
    {
        assert(nr()==itsPacker.nr());
        assert(nc()==itsPacker.nc());
    }
  
    size_t size() const { return  nr()*nc(); }
    size_t nr  () const { return std::ranges::size(a_rows); }
    size_t nc  () const { return std::ranges::size(b_cols); }

    value_type operator()(size_t i, size_t j) const
    {
        // assert(subsciptor.is_stored(i,j) && "Index out of range for MatrixView");
        return a_rows[i]*b_cols[j]; //VectorView*VectorView
    }
   

    auto rows() const
    {
        auto outerp=std::views::cartesian_product(a_rows,b_cols) | std::views::transform([](auto tuple) {return get<0>(tuple) * get<1>(tuple);}); // uses op*(const Vector<T>& v,const Matrix<T,S>& m)
        return outerp | std::views::chunk(nc()) | std::views::transform([](auto chunk) {return VectorView(std::move(chunk));});
    }

    auto cols() const
    {
        auto outerp=std::views::cartesian_product(a_rows,b_cols) | std::views::transform([](auto tuple) {return get<0>(tuple) * get<1>(tuple);}); // uses op*(const Vector<T>& v,const Matrix<T,S>& m)
        return  std::views::iota(size_t(0), nc()) | std::views::transform
            ([outerp,this](size_t j) 
                {
                    auto colj = outerp | std::views::drop(j) | std::views::stride(nc());
                    return VectorView(std::move(colj));
                }
            );
    }
    P packer() const {return itsPacker;}
    S shaper() const {return itsShaper;}

protected:
    R a_rows; //a as a range fo rows.
    C b_cols; //b as a range of cols.
    P itsPacker; //packing for the product.
    S itsShaper; // shape for the product
};

// Special version for full matrix products.  Skips indice interesction analysis the row[i]*col[j] dot products.
template <std::ranges::viewable_range R, std::ranges::viewable_range C> class FullMatrixCMProductView
: MatrixProductView<R,C,FullPackerCM,FullShaper>
{
public:
    using Base=MatrixProductView<R,C,FullPackerCM,FullShaper>;
    // These are all required in order to statisfy the isMatrix concept. With template classes they don't
    // automatically get pulled in from the base class :( 
    using value_type=Base::value_type;
    using Base::nr;
    using Base::nc;
    using Base::rows;
    using Base::cols;
    using Base::packer;
    using Base::shaper;
    //protected    
    using Base::a_rows;
    using Base::b_cols;

   FullMatrixCMProductView(const R& rows, const C& cols,FullPackerCM packer, FullShaper shaper )
    : Base(rows,cols,packer,shaper), i_cache(nr()) , ai_cache(0) {}
    value_type operator()(size_t i, size_t j) const
    {
        if (i!=i_cache)
        {
            size_t anc=a_rows[i].size(); //# of stored values in this row.
            if (ai_cache.size()!=anc) 
                ai_cache=default_data_type<value_type>(anc);
            i_cache=i;
            auto aij_cache=std::ranges::begin(ai_cache);
            for (auto aij:a_rows[i]) *aij_cache++=aij;
        }
        return inner_product(ai_cache,b_cols[j]); //Skip all indices intersections and checking
    }

private:
    mutable size_t i_cache;
    mutable default_data_type<value_type> ai_cache;
};

// general overloaded op* for matricies.
auto operator*(const isMatrix auto& a,const isMatrix auto& b)
{
    assert(a.nc() == b.nr() && "Matrix dimensions do not match for multiplication");
    auto p=MatrixProductPacker(a.packer(),b.packer());
    auto s=MatrixProductShaper(a.shaper(),b.shaper());
    return MatrixProductView(a.rows(),b.cols(),p,s);
}

// Special version for full matrix products.
template <class T> auto operator*(const FullMatrixCM<T>& a,const FullMatrixCM<T>& b)
{
    assert(a.nc() == b.nr() && "Matrix dimensions do not match for multiplication");
    auto p=MatrixProductPacker(a.packer(),b.packer());
    auto s=MatrixProductShaper(a.shaper(),b.shaper());
    return FullMatrixCMProductView(a.rows(),b.cols(),p,s); //Cache friendly version
}

} // namespace