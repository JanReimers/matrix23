// File: Matrix.hpp Define a matrix class for linear algebra operations.
#pragma once

#include "matrix23/subscriptors.hpp"
#include "matrix23/vector.hpp"

namespace matrix23
{

template <class M> 
concept isMatrix = requires (M m,size_t i, size_t j, std::remove_cvref_t<M>::value_t t)
{
    t=m.operator()(i, j);
    m.rows();
    m.cols();
    m.packer();
    m.shaper();
    m.nr();
    m.nc();
};

template <typename T, isPacker P, isShaper S> class Matrix
{

    public:
    typedef T value_t;
    using il_t = std::initializer_list<std::initializer_list<T>>;
    static size_t nr(const il_t& il) {return il.size();}
    static size_t nc(const il_t& il) {return il.begin()->size();}

    Matrix(P p) : itsPacker(p), itsShaper(itsPacker), data(itsPacker.stored_size()) {};
    Matrix(P p, S s) : itsPacker(p), itsShaper(s), data(itsPacker.stored_size()) {};
    Matrix(const il_t& init,P p) : Matrix(p)
    {
        load(init);
    }
    Matrix(const il_t& init,P p, S s) : Matrix(p,s)
    {
        load(init);
    }
    template <isMatrix M> Matrix(const M& m,P p, S s) 
        : itsPacker(p), itsShaper(s), data(itsPacker.stored_size())
    {
        for (size_t i = 0; i < nr(); ++i)
            for (size_t j = 0; j < nc(); ++j)
                if (itsPacker.is_stored(i, j)) 
                    (*this)(i,j) = m(i, j);
                else
                    assert(m(i,j)==0.0); //Make sure we are not throwing away data.
    }
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
    P packer() const {return itsPacker;}
    S shaper() const {return itsShaper;}

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
protected:
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



template <class T> class FullMatrixRM : public Matrix<T,FullPackerRM,FullShaper>
{
public:
    using Base = Matrix<T,FullPackerRM,FullShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    FullMatrixRM(size_t nr, size_t nc) : Base(FullPackerRM(nr,nc)) {};
    FullMatrixRM(const il_t& il) : Base(il,FullPackerRM(nr(il),nc(il))) {};
    template <isMatrix M> FullMatrixRM(const M& m) : Base(m,m.packer(), m.shaper()) {};
    // template <std::ranges::range Range> static void print(Range v)
    // {
    //     std::cout << "[";
    //     for (auto element : v) 
    //         std::cout << element << " ";

    //     std::cout << "]\n";
    // }
    auto rows() const //Assumes all rows are non-zero.
    {
        iota_view indices=std::views::iota(size_t(0), nc());
        return Base::data 
            | std::views::chunk(nc()) 
            | std::views::transform(
                [indices](auto row){return VectorView<decltype(row)>(std::move(row),indices);});
    }

    auto row(size_t i) const 
    {
        iota_view indices=std::views::iota(size_t(0), nc());
        auto rowi=Base::data | std::views::drop(i*nc()) | std::views::take(nc());
        return VectorView<decltype(rowi)>(std::move(rowi),indices);
    }

    auto col(size_t j) const 
    {
        iota_view indices=std::views::iota(size_t(0), nr());
        auto colj=Base::data | std::views::drop(j) | std::views::stride(nc());
        return VectorView<decltype(colj)>(std::move(colj),indices);
    }

};
template <class T> class FullMatrixCM : public Matrix<T,FullPackerCM,FullShaper>
{
public:
    using Base = Matrix<T,FullPackerCM,FullShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    FullMatrixCM(size_t nr, size_t nc) : Base(FullPackerCM(nr,nc)) {};
    FullMatrixCM(const il_t& il) : Base(il,FullPackerCM(nr(il),nc(il))) {};
    template <isMatrix M> FullMatrixCM(const M& m) : Base(m,m.packer(), m.shaper()) {};

    auto cols() const //Assumes all rows are non-zero.
    {
        iota_view indices=std::views::iota(size_t(0), nr());
        return Base::data 
            | std::views::chunk(nr()) 
            | std::views::transform(
                [indices](auto col){return VectorView<decltype(col)>(std::move(col),indices);});
    }

    auto row(size_t i) const 
    {
        iota_view indices=std::views::iota(size_t(0), nc());
        auto rowi=Base::data | std::views::drop(i) | std::views::stride(nr());
        return VectorView<decltype(rowi)>(std::move(rowi),indices);
    }

    auto col(size_t j) const 
    {
        iota_view indices=std::views::iota(size_t(0), nr());
        auto colj=Base::data | std::views::drop(j*nr()) | std::views::take(nr());
        return VectorView<decltype(colj)>(std::move(colj),indices);
    }
};

template <class T> class UpperTriangularMatrixCM : public Matrix<T,UpperTriangularPackerCM,UpperTriangularShaper>
{
public:
    using Base = Matrix<T,UpperTriangularPackerCM,UpperTriangularShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    UpperTriangularMatrixCM(size_t nr, size_t nc) : Base(UpperTriangularPackerCM(nr,nc)) {};
    UpperTriangularMatrixCM(const il_t& il) : Base(il,UpperTriangularPackerCM(nr(il),nc(il))) {};
    template <isMatrix M> UpperTriangularMatrixCM(const M& m) : Base(m,m.packer(), m.shaper()) {};
};
template <class T> class UpperTriangularMatrixRM : public Matrix<T,UpperTriangularPackerRM,UpperTriangularShaper>
{
public:
    using Base = Matrix<T,UpperTriangularPackerRM,UpperTriangularShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    UpperTriangularMatrixRM(size_t nr, size_t nc) : Base(UpperTriangularPackerRM(nr,nc)) {};
    UpperTriangularMatrixRM(const il_t& il) : Base(il,UpperTriangularPackerRM(nr(il),nc(il))) {};
    template <isMatrix M> UpperTriangularMatrixRM(const M& m) : Base(m,m.packer(), m.shaper()) {};
};

template <class T> class LowerTriangularMatrixCM : public Matrix<T,LowerTriangularPackerCM,LowerTriangularShaper>
{
public:
    using Base = Matrix<T,LowerTriangularPackerCM,LowerTriangularShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    LowerTriangularMatrixCM(size_t nr, size_t nc) : Base(LowerTriangularPackerCM(nr,nc)) {};
    LowerTriangularMatrixCM(const il_t& il) : Base(il,LowerTriangularPackerCM(nr(il),nc(il))) {};
    template <isMatrix M> LowerTriangularMatrixCM(const M& m) : Base(m,m.packer(), m.shaper()) {};
};
template <class T> class LowerTriangularMatrixRM : public Matrix<T,LowerTriangularPackerRM,LowerTriangularShaper>
{
public:
    using Base = Matrix<T,LowerTriangularPackerRM,LowerTriangularShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    LowerTriangularMatrixRM(size_t nr, size_t nc) : Base(LowerTriangularPackerRM(nr,nc)) {};
    LowerTriangularMatrixRM(const il_t& il) : Base(il,LowerTriangularPackerRM(nr(il),nc(il))) {};
    template <isMatrix M> LowerTriangularMatrixRM(const M& m) : Base(m,m.packer(), m.shaper()) {};
};
template <class T> class DiagonalMatrix : public Matrix<T,DiagonalPacker,DiagonalShaper>
{
public:
    using Base = Matrix<T,DiagonalPacker,DiagonalShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    DiagonalMatrix(size_t nr, size_t nc) : Base(DiagonalPacker(nr,nc)) {};
    DiagonalMatrix(const il_t& il) : Base(il,DiagonalPacker(nr(il),nc(il))) {};
    template <isMatrix M> DiagonalMatrix(const M& m) : Base(m,m.packer(), m.shaper()) {};
};
template <class T> class SBandMatrix : public Matrix<T,SBandPacker,SBandShaper>
{
public:
    using Base = Matrix<T,SBandPacker,SBandShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    SBandMatrix(size_t n, size_t k) : Base(SBandPacker(n,k)) 
    {
        
    };
    SBandMatrix(const il_t& il,size_t k) : Base(il,SBandPacker(nr(il),k))
    {
        assert(nr(il)==nc(il)); //SBand only supports square matricies.
    };
    template <isMatrix M> SBandMatrix(const M& m) : Base(m,m.packer(), m.shaper()) {};
};

auto operator*(const isMatrix auto& m, const isVector auto& v)
{
    assert(m.nc()==v.size());
    auto rows=m.rows();
    auto indices=std::views::iota(size_t(0),rows.size());
    auto mv= indices | std::views::transform([rows,v](size_t i) {return rows[i] * v;});
    assert(mv.size()==rows.size());
    static_assert(isVector<VectorView<decltype(mv)>>,"Matrix-vector multiplication should satisfy isVector concept requirements");
    return VectorView(std::move(mv), indices);
}

auto operator*(const isVector auto& v,const isMatrix auto& m)
{
    auto cols=m.cols();
    auto indices=std::views::iota(size_t(0),cols.size());
    auto vm=indices | std::views::transform([cols,v](size_t j) {return v * cols[j];});
    assert(vm.size()==cols.size());
    static_assert(isVector<VectorView<decltype(vm)>>,"Matrix-vector multiplication should satisfy isVector concept requirements");
    return VectorView(std::move(vm), indices);
}

template <std::ranges::viewable_range R, std::ranges::viewable_range C, isPacker P, isShaper S> class MatrixProductView
{
public:
    typedef std::ranges::range_value_t<R> Rv; //rows value type which is a column range.
    typedef std::ranges::range_value_t<Rv> value_t; //column range value type which should be scalar (double etc.)
    MatrixProductView(const R& _rows, const C& _cols,P _packer, S _shaper )
    : a_rows(_rows), b_cols(_cols), itsPacker(_packer), itsShaper(_shaper)
    {
        assert(nr()==itsPacker.nr());
        assert(nc()==itsPacker.nc());
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
    P packer() const {return itsPacker;}
    S shaper() const {return itsShaper;}

private:
    R a_rows; //a as a range fo rows.
    C b_cols; //b as a range of cols.
    P itsPacker; //packing for the product.
    S itsShaper; // shape for the product
};

template <isPacker A, isPacker B> auto MatrixProductPacker(const A& a, const B& b)
{
    using packer_t=MatrixProductPackerType<A,B>::packer_t;
    return packer_t(a.nr(),b.nc());
}
template <> inline auto MatrixProductPacker(const SBandPacker& a, const SBandPacker& b)
{
    assert(a.nr()==b.nr());
    return SBandPacker(a.nr(),a.bandwidth()+b.bandwidth());
}

auto operator*(const isMatrix auto& a,const isMatrix auto& b)
{
    assert(a.nc() == b.nr() && "Matrix dimensions do not match for multiplication");
    // using packer_t=MatrixProductPackerType<decltype(a.packer()),decltype(b.packer())>::packer_t;
    using shaper_t=MatrixProductShaperType<decltype(a.shaper()),decltype(b.shaper())>::shaper_t;
    auto p=MatrixProductPacker(a.packer(),b.packer());
    shaper_t s=shaper_t(p);
    return MatrixProductView(a.rows(),b.cols(),p,s);
}



} //namespace matrix23