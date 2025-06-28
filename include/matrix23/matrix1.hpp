// File: Matrix.hpp Define a matrix class for linear algebra operations.
#pragma once

#include "matrix23/matrix.hpp"
#include "matrix23/subscriptors.hpp"

namespace matrix23
{

auto make_packer(const size_t& nr,const size_t& nc, Packing p)
{

}
template <typename T, isPacker P, isShaper S> class Matrix1
{

    public:
    typedef T value_t;
    using il_t = std::initializer_list<std::initializer_list<T>>;
    static size_t nr(const il_t& il) {return il.size();}
    static size_t nc(const il_t& il) {return il.begin()->size();}

    Matrix1(P p, S s) : itsPacker(p), itsShaper(s), data(itsPacker.stored_size()) {};
    Matrix1(const il_t& init,P p, S s) : Matrix1(p,s)
    {
        load(init);
    }
    // Matrix1(std::initializer_list<std::initializer_list<T>> init, size_t k)
    // : itsSubscriptor(init.size(), init.begin()->size(),k), data(itsSubscriptor.size())
    // {
    //     load(init);
    // }
    // template <isMatrix M> Matrix(const M& m) 
    // : itsSubscriptor(m.nr(),m.nc()), data(itsSubscriptor.size())
    // {
    //     for (size_t i = 0; i < itsSubscriptor.nr(); ++i)
    //     {
    //         for (size_t j = 0; j < itsSubscriptor.nc(); ++j)
    //         {
    //             if (itsSubscriptor.is_stored(i, j))
    //                 data[itsSubscriptor.offset(i, j)] = m(i, j);
                
    //         }
    //     }
    // }
    
    //
    //  2D data access.
    //
    T  operator()(size_t i, size_t j) const
    {
        // std::cout << "i,j=" << i << " " << j << "  offfset=" << itsPacker.offset(i,j) << std::endl;
        return itsPacker.is_stored(i,j) ? data[itsPacker.offset(i,j)] : T(0);
    }
    T& operator()(size_t i, size_t j)
    {
        assert(itsPacker.is_stored(i,j));
        return data[itsPacker.offset(i,j)];
    }
    size_t nr() const { return itsPacker.nr(); }
    size_t nc() const { return itsPacker.nc(); }

    // 1D range access.
    size_t size() const
    {
        return itsPacker.stored_size();
    }
    auto begin()       { return std::begin(data); }
    auto end  ()       { return std::end  (data); }
    auto begin() const { return std::begin(data); }
    auto end  () const { return std::end  (data); }
    //
    //  2D range access.
    //
    auto row(size_t i) const
    {
        assert(i < itsPacker.nr());
        auto indices=itsShaper.nonzero_col_indexes(i);
        auto v=  indices | std::views::transform([i,this](size_t j){return operator()(i,j);});
        return VectorView<decltype(v)>(std::move(v),indices);
    }
    auto col(size_t j) const
    {
        assert(j < itsPacker.nc());
        auto indices=itsShaper.nonzero_row_indexes(j);
        auto v= indices | std::views::transform([j,this](size_t i){return operator()(i,j);});
        return VectorView<decltype(v)>(std::move(v),indices);
    }
    auto rows() const //Assumes all rows are non-zero.
    {
        return std::views::iota(size_t(0), itsPacker.nr()) | std::views::transform([this](size_t i){return row(i);});
    }
    auto cols() const //Assumes all cols are non-zero.
    {
        return std::views::iota(size_t(0), itsPacker.nc()) | std::views::transform([this](size_t j){return col(j);});
    }
    P subscriptor() const
    {
        return itsPacker;
    }
    void print() const
    {
        for (size_t i = 0; i < nr(); ++i)
        {
            for (size_t j = 0; j < nc(); ++j)
            {
                // if (itsPacker.is_stored(i, j))
                    std::cout << (*this)(i, j) << " ";
                // else
                //     std::cout << "0 "; // Print 0 for non-stored elements
            }
            std::cout << std::endl;
        }   
    }
private:
    void load(std::initializer_list<std::initializer_list<T>> init)
    {
        assert(init.size() == nr() && "Initializer list size does not match subscriptor row count");
        assert(init.begin()->size() == nc() && "Initializer list row size does not match subscriptor column count");
        size_t i = 0;
        for (const auto& row : init)
        {
            size_t j = 0;
            for (const auto& val : row)
            {
                if (itsPacker.is_stored(i, j))
                    (*this)(i, j) = val;
                else
                    assert(val==0.0);
                ++j;
            }
            ++i;
        }
    }

    P itsPacker;
    S itsShaper;
    std::valarray<T> data;
};


template <class T> class FullMatrix1 : public Matrix1<T,FullPacker,FullShaper>
{
    public:
        using Base = Matrix1<T,FullPacker,FullShaper>;
        using il_t=Base::il_t;
        using Base::nr;
        using Base::nc;
        FullMatrix1(size_t nr, size_t nc, Indexing ind=Indexing::col_major) : Base(FullPacker(nr,nc,ind), FullShaper(nr,nc)) {};
        FullMatrix1(const il_t& il, Indexing ind=Indexing::col_major) 
            : Base(il,FullPacker(nr(il),nc(il),ind), FullShaper(nr(il),nc(il))) {};
};
template <class T> class UpperTriangularMatrix1 : public Matrix1<T,UpperTriangularPacker,UpperTriangularShaper>
{
    public:
        using Base = Matrix1<T,UpperTriangularPacker,UpperTriangularShaper>;
        using il_t=Base::il_t;
        using Base::nr;
        using Base::nc;
        UpperTriangularMatrix1(size_t nr, size_t nc, Indexing ind=Indexing::col_major) : Base(UpperTriangularPacker(nr,nc,ind), UpperTriangularShaper(nr,nc)) {};
        UpperTriangularMatrix1(const il_t& il, Indexing ind=Indexing::col_major) : Base(il,UpperTriangularPacker(nr(il),nc(il),ind), UpperTriangularShaper(nr(il),nc(il))) {};
};
template <class T> class LowerTriangularMatrix1 : public Matrix1<T,LowerTriangularPacker,LowerTriangularShaper>
{
    public:
        using Base = Matrix1<T,LowerTriangularPacker,LowerTriangularShaper>;
        using il_t=Base::il_t;
        using Base::nr;
        using Base::nc;
        LowerTriangularMatrix1(size_t nr, size_t nc, Indexing ind=Indexing::col_major) : Base(LowerTriangularPacker(nr,nc,ind), LowerTriangularShaper(nr,nc)) {};
        LowerTriangularMatrix1(const il_t& il, Indexing ind=Indexing::col_major) : Base(il,LowerTriangularPacker(nr(il),nc(il),ind), LowerTriangularShaper(nr(il),nc(il))) {};
};
template <class T> class DiagonalMatrix1 : public Matrix1<T,DiagonalPacker,DiagonalShaper>
{
    public:
        using Base = Matrix1<T,DiagonalPacker,DiagonalShaper>;
        using il_t=Base::il_t;
        using Base::nr;
        using Base::nc;
        DiagonalMatrix1(size_t nr, size_t nc) : Base(DiagonalPacker(nr,nc), DiagonalShaper(nr,nc)) {};
        DiagonalMatrix1(const il_t& il) : Base(il,DiagonalPacker(nr(il),nc(il)), DiagonalShaper(nr(il),nc(il))) {};
};
template <class T> class SBandMatrix1 : public Matrix1<T,SBandPacker,SBandShaper>
{
    public:
        using Base = Matrix1<T,SBandPacker,SBandShaper>;
        using il_t=Base::il_t;
        using Base::nr;
        using Base::nc;
        SBandMatrix1(size_t n, size_t k) : Base(SBandPacker(n,k), SBandShaper(n,k)) 
        {
            
        };
        SBandMatrix1(const il_t& il,size_t k) : Base(il,SBandPacker(nr(il),k), SBandShaper(nr(il),k))
        {
            assert(nr(il)==nc(il)); //SBand only supports square matricies.
        };
};

// typedef Matrix1<double,    FullRowMajorSubsciptor>            FullMatrix;
// typedef Matrix<double, UpperTriangularRowMajorSubsciptor> UpperTriangularMatrix;
// typedef Matrix<double,        DiagonalSubsciptor>        DiagonalMatrix;
// typedef Matrix<double,     TriDiagonalSubsciptor>     TriDiagonalMatrix;
// // typedef Matrix<double, BandedSubsciptor> BandedMatrix; // 

// auto operator*(const isMatrix auto& m, const isVector auto& v)
// {
//     auto rows=m.rows();
//     auto mv=v.indices() | std::views::transform([rows,v](size_t i) {return rows[i] * v;});
//     static_assert(isVector<VectorView<decltype(mv)>>,"Matrix-vector multiplication should satisfy isVector concept requirements");
//     return VectorView(std::move(mv), v.indices());
// }

// auto operator*(const isVector auto& v,const isMatrix auto& m)
// {
//     auto cols=m.cols();
//     auto vm=v.indices() | std::views::transform([cols,v](size_t j) {return v * cols[j];});
//     static_assert(isVector<VectorView<decltype(vm)>>,"Matrix-vector multiplication should satisfy isVector concept requirements");
//     return VectorView(std::move(vm), v.indices());
// }

// template <std::ranges::viewable_range R, std::ranges::viewable_range C, class S> class MatrixProductView
// {
// public:
//     typedef std::ranges::range_value_t<R> Rv; //rows value type which is a column range.
//     typedef std::ranges::range_value_t<Rv> value_t; //column range value type which should be scalar (double etc.)
//     MatrixProductView(const R& _rows, const C& _cols,S _subsciptor )
//     : a_rows(_rows), b_cols(_cols), itsSubscriptor(_subsciptor)
//     {
//         assert(nr()==itsSubscriptor.nr());
//         assert(nc()==itsSubscriptor.nc());
//     }
  
//     size_t size() const { return  nr()*nc(); }
//     size_t nr  () const { return std::ranges::size(a_rows); }
//     size_t nc  () const { return std::ranges::size(b_cols); }

//     value_t operator()(size_t i, size_t j) const
//     {
//         // assert(subsciptor.is_stored(i,j) && "Index out of range for MatrixView");
//         return a_rows[i]*b_cols[j]; //VectorView*VectorView
//     }

//     auto rows() const
//     {
//         auto outerp=std::views::cartesian_product(a_rows,b_cols) | std::views::transform([](auto tuple) {return get<0>(tuple) * get<1>(tuple);}); // uses op*(const Vector<T>& v,const Matrix<T,S>& m)
//         return outerp | std::views::chunk(nc()) | std::views::transform([](auto chunk) {return VectorView(std::move(chunk));});
//     }

//     auto cols() const
//     {
//         auto outerp=std::views::cartesian_product(a_rows,b_cols) | std::views::transform([](auto tuple) {return get<0>(tuple) * get<1>(tuple);}); // uses op*(const Vector<T>& v,const Matrix<T,S>& m)
//         return  std::views::iota(size_t(0), nc()) | std::views::transform
//             ([outerp,this](size_t j) 
//                 {
//                     auto colj = outerp | std::views::drop(j) | std::views::stride(nc());
//                     return VectorView(std::move(colj));
//                 }
//             );
//     }
//     S subscriptor() const
//     {
//         return itsSubscriptor;
//     }
// private:
//     R a_rows; //a as a range fo rows.
//     C b_cols; //b as a range of cols.
//     S itsSubscriptor; //packing for the product.
// };

// auto operator*(const isMatrix auto& a,const isMatrix auto& b)
// {
//     assert(a.nc() == b.nr() && "Matrix dimensions do not match for multiplication");
//     static_assert(isMatrix<MatrixProductView<decltype(a.rows()),decltype(b.rows()),decltype(a.subscriptor())>>,"Matrix-Matrix multiplication should satisfy isMatrix concept requirements");
//     return MatrixProductView(a.rows(),b.cols(),a.subscriptor());
// }

// //
// //  operator+ introduces some interesting issues.
// //    1) For Matrix we want to do linear/flat traverse of the underlying 1D data.  i.e. no expensive op(i,j) indexing.
// //    2) For MatrixProductView we can add rows which should skip zero-elements as dictated by the subscriptor.
// //
// // template <class T, class S> auto operator+(const Matrix<T,S>& a, const Matrix<T,S>& b)
// // {
// //     auto ab=std::views::zip_transform([](const auto& ia, const auto& ib) { return ia + ib; },a,b);
// //     std::op_plus<T>;
// //     return ???;
// // }

// template <isMatrix Ma, isMatrix Mb, class Op> class MatrixBinOpView
// {
// public:
//     typedef std::remove_cvref_t<Ma>::value_t value_t;
//     MatrixBinOpView(const Ma& _a, const Mb& _b, const Op& _op)
//     : a(_a), b(_b), op(_op)
//     {
//         assert(a.nr()==b.nr());
//         assert(a.nc()==b.nc());
//     }
  
//     size_t size() const { return  nr()*nc(); }
//     size_t nr  () const { return a.nr(); }
//     size_t nc  () const { return a.nc(); }

//     value_t operator()(size_t i, size_t j) const
//     {
//         return op(a(i,j),b(i,j));
//     }

//     auto rows() const
//     {
//        return std::views::zip_transform(op,a.rows(),b.rows());
//     }

//     auto cols() const
//     {
//        return std::views::zip_transform(op,a.cols(),b.cols());
//     }
//     auto subscriptor() const
//     {
//         return a.subscriptor();;
//     }
// private:
//     Ma a; 
//     Mb b; 
//     Op op; 
// };

// template <isMatrix Ma, class Op> class MatrixOpView
// {
// public:
//     typedef std::remove_cvref_t<Ma>::value_t value_t;
//     MatrixOpView(const Ma& _a, const Op& _op)
//     : a(_a), op(_op)
//     {
       
//     }
  
//     size_t size() const { return  nr()*nc(); }
//     size_t nr  () const { return a.nr(); }
//     size_t nc  () const { return a.nc(); }

//     value_t operator()(size_t i, size_t j) const
//     {
//         return op(a(i,j));
//     }

//     auto rows() const
//     {
//        return std::views::transform(a.rows(),op);
//     }

//     auto cols() const
//     {
//        return std::views::transform(a.cols(),op);
//     }
//     auto subscriptor() const
//     {
//         return a.subscriptor();;
//     }
// private:
//     Ma a; 
//     Op op; 
// };


// auto operator+(const isMatrix auto& a,const isMatrix auto& b)
// {
//     static_assert(isMatrix<MatrixBinOpView<decltype(a),decltype(b),decltype([](const auto& ia, const auto& ib){return ia+ib;})>>,"Matrix-Matrix addition should satisfy isMatrix concept requirements");
//     return MatrixBinOpView(a,b,[](const auto& ia, const auto& ib){return ia+ib;});                    
// }

// auto operator-(const isMatrix auto& a,const isMatrix auto& b)
// {
//     return MatrixBinOpView(a,b,[](const auto& ia, const auto& ib){return ia-ib;});                    
// }

// auto operator*(const isMatrix auto& a,const arithmetic auto& b)
// {
//     return MatrixOpView(a,[b](const auto& ia){return ia*b;}); 
// }
// auto operator*(const arithmetic auto& b,const isMatrix auto& a)
// {
//     return MatrixOpView(a,[b](const auto& ia){return b*ia;}); 
// }
// auto operator/(const isMatrix auto& a,const arithmetic auto& b)
// {
//     return MatrixOpView(a,[b](const auto& ia){return ia/b;}); 
// }


} //namespace matrix23