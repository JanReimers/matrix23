// File: Matrix.hpp Define a matrix class for linear algebra operations.
#pragma once

#include "matrix23/vector.hpp"
#include "matrix23/subscriptors.hpp"

namespace matrix23
{

template <class M> 
concept isMatrix = requires (M m,size_t i, size_t j, std::remove_cvref_t<M>::value_t t)
{
    t=m.operator()(i, j);
    m.rows();
    m.cols();
    m.subscriptor();
    m.nr();
    m.nc();
};

template <typename T, isSubscriptor S> class Matrix
{
    public:
    typedef T value_t;
    Matrix(size_t nr, size_t nc) : itsSubscriptor(nc,nr), data(itsSubscriptor.size()) {};
    Matrix(size_t nr, size_t nc, size_t k) : itsSubscriptor(nc,nr,k), data(itsSubscriptor.size()) {};
    Matrix(std::initializer_list<std::initializer_list<T>> init)
        : itsSubscriptor(init.size(), init.begin()->size()), data(itsSubscriptor.size())
    {
        load(init);
    }
    Matrix(std::initializer_list<std::initializer_list<T>> init, size_t k)
    : itsSubscriptor(init.size(), init.begin()->size(),k), data(itsSubscriptor.size())
    {
        load(init);
    }
    template <isMatrix M> Matrix(const M& m) 
    : itsSubscriptor(m.nr(),m.nc()), data(itsSubscriptor.size())
    {
        for (size_t i = 0; i < itsSubscriptor.nr; ++i)
        {
            for (size_t j = 0; j < itsSubscriptor.nc; ++j)
            {
                if (itsSubscriptor.is_stored(i, j))
                    data[itsSubscriptor.offset(i, j)] = m(i, j);
                
            }
        }
    }
    
    T  operator()(size_t i, size_t j) const
    {
        return itsSubscriptor.is_stored(i,j) ? data[itsSubscriptor.offset(i,j)] : T(0);
    }
    T& operator()(size_t i, size_t j)
    {
        assert(itsSubscriptor.is_stored(i,j));
        return data[itsSubscriptor.offset(i,j)];
    }
    size_t size() const
    {
        return itsSubscriptor.size();
    }
    size_t nr() const { return itsSubscriptor.nr; }
    size_t nc() const { return itsSubscriptor.nc; }
    // auto begin()       { return std::begin(data); }
    // auto end  ()       { return std::end  (data); }
    // auto begin() const { return std::begin(data); }
    // auto end  () const { return std::end  (data); }

    auto row(size_t i) const
    {
        assert(i < itsSubscriptor.nr);
        auto indices=itsSubscriptor.nonzero_col_indexes(i);
        auto v=  indices | std::views::transform([i,this](size_t j){return operator()(i,j);});
        return VectorView<decltype(v)>(std::move(v),indices);
    }
    auto col(size_t j) const
    {
        assert(j < itsSubscriptor.nc);
        auto indices=itsSubscriptor.nonzero_row_indexes(j);
        auto v= indices | std::views::transform([j,this](size_t i){return operator()(i,j);});
        return VectorView<decltype(v)>(std::move(v),indices);
    }
    auto rows() const //Assumes all rows are non-zero.
    {
        return std::views::iota(size_t(0), itsSubscriptor.nr) | std::views::transform([this](size_t i){return row(i);});
    }
    auto cols() const //Assumes all cols are non-zero.
    {
        return std::views::iota(size_t(0), itsSubscriptor.nc) | std::views::transform([this](size_t j){return col(j);});
    }
    S subscriptor() const
    {
        return itsSubscriptor;
    }
    void print() const
    {
        for (size_t i = 0; i < itsSubscriptor.nr; ++i)
        {
            for (size_t j = 0; j < itsSubscriptor.nc; ++j)
            {
                if (itsSubscriptor.is_stored(i, j))
                    std::cout << (*this)(i, j) << " ";
                else
                    std::cout << "0 "; // Print 0 for non-stored elements
            }
            std::cout << std::endl;
        }   
    }
private:
    void load(std::initializer_list<std::initializer_list<T>> init)
    {
        assert(init.size() == itsSubscriptor.nr && "Initializer list size does not match subscriptor row count");
        assert(init.begin()->size() == itsSubscriptor.nc && "Initializer list row size does not match subscriptor column count");
        size_t i = 0;
        for (const auto& row : init)
        {
            size_t j = 0;
            for (const auto& val : row)
            {
                if (itsSubscriptor.is_stored(i, j))
                    data[itsSubscriptor.offset(i, j)] = val;
                else
                    assert(val==0.0);
                ++j;
            }
            ++i;
        }
    }

    S itsSubscriptor;
    std::valarray<T> data;
};

typedef Matrix<double,            FullSubsciptor>            FullMatrix;
typedef Matrix<double, UpperTriangularSubsciptor> UpperTriangularMatrix;
typedef Matrix<double,        DiagonalSubsciptor>        DiagonalMatrix;
typedef Matrix<double,     TriDiagonalSubsciptor>     TriDiagonalMatrix;
// typedef Matrix<double, BandedSubsciptor> BandedMatrix; // 

auto operator*(const isMatrix auto& m, const isVector auto& v)
{
    auto rows=m.rows();
    auto mv=v.indices() | std::views::transform([rows,v](size_t i) {return rows[i] * v;});
    static_assert(isVector<VectorView<decltype(mv)>>,"Matrix-vector multiplication should satisfy isVector concept requirements");
    return VectorView(std::move(mv), v.indices());
}

auto operator*(const isVector auto& v,const isMatrix auto& m)
{
    auto cols=m.cols();
    auto vm=v.indices() | std::views::transform([cols,v](size_t j) {return v * cols[j];});
    static_assert(isVector<VectorView<decltype(vm)>>,"Matrix-vector multiplication should satisfy isVector concept requirements");
    return VectorView(std::move(vm), v.indices());
}

template <std::ranges::viewable_range R, std::ranges::viewable_range C, class S> class MatrixProductView
{
public:
    typedef std::ranges::range_value_t<R> Rv; //rows value type which is a column range.
    typedef std::ranges::range_value_t<Rv> value_t; //column range value type which should be scalar (double etc.)
    MatrixProductView(const R& _rows, const C& _cols,S _subsciptor )
    : a_rows(_rows), b_cols(_cols), itsSubscriptor(_subsciptor)
    {
        assert(nr()==itsSubscriptor.nr);
        assert(nc()==itsSubscriptor.nc);
    }
  
    size_t size() const { return  nr()*nc(); }
    size_t nr  () const { return std::ranges::size(a_rows); }
    size_t nc  () const { return std::ranges::size(b_cols); }

    value_t operator()(size_t i, size_t j) const
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
    S subscriptor() const
    {
        return itsSubscriptor;
    }
private:
    R a_rows; //a as a range fo rows.
    C b_cols; //b as a range of cols.
    S itsSubscriptor; //packing for the product.
};

auto operator*(const isMatrix auto& a,const isMatrix auto& b)
{
    assert(a.nc() == b.nr() && "Matrix dimensions do not match for multiplication");
    static_assert(isMatrix<MatrixProductView<decltype(a.rows()),decltype(b.rows()),decltype(a.subscriptor())>>,"Matrix-Matrix multiplication should satisfy isMatrix concept requirements");
    return MatrixProductView(a.rows(),b.cols(),a.subscriptor());
}

//
//  operator+ introduces some interesting issues.
//    1) For Matrix we want to do linear/flat traverse of the underlying 1D data.  i.e. no expensive op(i,j) indexing.
//    2) For MatrixProductView we can add rows which should skip zero-elements as dictated by the subscriptor.
//
// template <class T, class S> auto operator+(const Matrix<T,S>& a, const Matrix<T,S>& b)
// {
//     auto ab=std::views::zip_transform([](const auto& ia, const auto& ib) { return ia + ib; },a,b);
//     std::op_plus<T>;
//     return ???;
// }

template <isMatrix Ma, isMatrix Mb, class Op> class MatrixBinOpView
{
public:
    typedef std::remove_cvref_t<Ma>::value_t value_t;
    MatrixBinOpView(const Ma& _a, const Mb& _b, const Op& _op)
    : a(_a), b(_b), op(_op)
    {
        assert(a.nr()==b.nr());
        assert(a.nc()==b.nc());
    }
  
    size_t size() const { return  nr()*nc(); }
    size_t nr  () const { return a.nr(); }
    size_t nc  () const { return a.nc(); }

    value_t operator()(size_t i, size_t j) const
    {
        return op(a(i,j),b(i,j));
    }

    auto rows() const
    {
       return std::views::zip_transform(op,a.rows(),b.rows());
    }

    auto cols() const
    {
       return std::views::zip_transform(op,a.cols(),b.cols());
    }
    auto subscriptor() const
    {
        return a.subscriptor();;
    }
private:
    Ma a; 
    Mb b; 
    Op op; 
};

template <isMatrix Ma, class Op> class MatrixOpView
{
public:
    typedef std::remove_cvref_t<Ma>::value_t value_t;
    MatrixOpView(const Ma& _a, const Op& _op)
    : a(_a), op(_op)
    {
       
    }
  
    size_t size() const { return  nr()*nc(); }
    size_t nr  () const { return a.nr(); }
    size_t nc  () const { return a.nc(); }

    value_t operator()(size_t i, size_t j) const
    {
        return op(a(i,j));
    }

    auto rows() const
    {
       return std::views::transform(a.rows(),op);
    }

    auto cols() const
    {
       return std::views::transform(a.cols(),op);
    }
    auto subscriptor() const
    {
        return a.subscriptor();;
    }
private:
    Ma a; 
    Op op; 
};


auto operator+(const isMatrix auto& a,const isMatrix auto& b)
{
    static_assert(isMatrix<MatrixBinOpView<decltype(a),decltype(b),decltype([](const auto& ia, const auto& ib){return ia+ib;})>>,"Matrix-Matrix addition should satisfy isMatrix concept requirements");
    return MatrixBinOpView(a,b,[](const auto& ia, const auto& ib){return ia+ib;});                    
}

auto operator-(const isMatrix auto& a,const isMatrix auto& b)
{
    return MatrixBinOpView(a,b,[](const auto& ia, const auto& ib){return ia-ib;});                    
}

auto operator*(const isMatrix auto& a,const arithmetic auto& b)
{
    return MatrixOpView(a,[b](const auto& ia){return ia*b;}); 
}
auto operator*(const arithmetic auto& b,const isMatrix auto& a)
{
    return MatrixOpView(a,[b](const auto& ia){return b*ia;}); 
}
auto operator/(const isMatrix auto& a,const arithmetic auto& b)
{
    return MatrixOpView(a,[b](const auto& ia){return ia/b;}); 
}


} //namespace matrix23