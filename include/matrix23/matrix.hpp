// File: Matrix.hpp Define a matrix class for linear algebra operations.
#pragma once

#include "matrix23/vector.hpp"
#include "matrix23/subscriptors.hpp"

namespace matrix23
{

template <class M> 
concept isMatrix = requires (M m,size_t i, size_t j)
{
    m.operator()(i, j);
    m.rows();
    m.cols();
    m.subscriptor();
    m.nr();
    m.nc();
};

template <typename T, class S> class Matrix
{
    public:
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
    T  operator()(size_t i, size_t j) const
    {
        assert(itsSubscriptor.is_stored(i,j));
        return data[itsSubscriptor.offset(i,j)];
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

template <isMatrix M, isVector V> auto operator*(const M& m, const V& v)
{
    auto rows=m.rows();
    auto mv=v.indices() | std::views::transform([rows,v](size_t i) {return rows[i] * v;});
    return VectorView(std::move(mv), v.indices());
}

template <isVector V, isMatrix M> auto operator*(const V& v,const M& m)
{
    auto cols=m.cols();
    auto vm=v.indices() | std::views::transform([cols,v](size_t j) {return v * cols[j];});
    return VectorView(std::move(vm), v.indices());
}

template <std::ranges::viewable_range R, std::ranges::viewable_range C, class S> class MatrixProductView
{
public:
    typedef std::ranges::range_value_t<R> Rv; //rows value type which is a column range.
    typedef std::ranges::range_value_t<Rv> T; //column range value type which should be scalar (double etc.)
    MatrixProductView(const R& _rows, const C& _cols,S _subsciptor )
    : a_rows(_rows), b_cols(_cols), itsSubscriptor(_subsciptor)
    {
        assert(nr()==itsSubscriptor.nr);
        assert(nc()==itsSubscriptor.nc);
    }
  
    size_t size() const { return  nr()*nc(); }
    size_t nr  () const { return std::ranges::size(a_rows); }
    size_t nc  () const { return std::ranges::size(b_cols); }

    T operator()(size_t i, size_t j) const
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

template <isMatrix A, isMatrix B> auto operator*(const A& a,const B& b)
{
    assert(a.nc() == b.nr() && "Matrix dimensions do not match for multiplication");
    return MatrixProductView(a.rows(),b.cols(),a.subscriptor());
}





} //namespace matrix23