// File: Matrix.hpp Define a matrix class for linear algebra operations.
#pragma once

#include "matrix23/vector.hpp"
#include "matrix23/shaper.hpp"
#include "matrix23/packer.hpp"
#include "matrix23/symmetry.hpp"
#include <iostream>

namespace matrix23
{

template <class M> 
concept isMatrix = requires (M m,size_t i, size_t j, std::remove_cvref_t<M>::value_type t)
{
    t=m.operator()(i, j);
    m.rows();
    m.cols();
    m.packer();
    m.shaper();
    m.nr();
    m.nc();
};


//default_data_type is defined in vector.hpp.

template <typename T, isPacker P, isShaper S, typename D=default_data_type<T>, isSymmetry Sym=NoSymmetry<D,P> > class Matrix
{
public:
    typedef T value_type;
    using il_t = std::initializer_list<std::initializer_list<T>>;
    static size_t nr(const il_t& il) {return il.size();}
    static size_t nc(const il_t& il) {return il.begin()->size();}

    Matrix(                    ) : Matrix( 0) {}; //Empty matrix.
    Matrix(size_t n            ) : Matrix(n,n) {}; //Square matrix.
    Matrix(size_t nr, size_t nc) : Matrix(nr,nc,none) {}; //Default to no fill.
    Matrix(size_t n            , fill_t f, T v=T(1)) : Matrix(n,n,f,v) {}; //Square matrix with fill.
    Matrix(size_t nr, size_t nc, fill_t f, T v=T(1)) : Matrix(P(nr,nc),S(nr,nc),f,v) {} //Rect matrix with fill
    Matrix(const il_t& init) : Matrix(P(nr(init),nc(init)),S(nr(init),nc(init)),init) {} //Fill from std::initial_list
    template <isMatrix M>  Matrix(const M& m) : Matrix(P(m.nr(),m.nc()),m.shaper(),m) {} //assign from an expression.

    // Banded matricies have extra paramater(s) k,ku,kl etc. As such they need to pass down already created
    // packers in order to handle these extra parameters.
    Matrix(P p, fill_t f, T v=T(1)) : Matrix(p,p.shaper(),f,v) {}
    Matrix(P p,const il_t& init) : Matrix(p,p.shaper(),init) {} //SBand needs this, only p knows k.
    template <isMatrix M> Matrix(P p,const M& m) : Matrix(p,m.shaper(),m) {} //SBand needs this, only p knows k.

private:
    // All of the public constructors should lead to these generic versions
    Matrix(P p, S s,const il_t& init) : Matrix(p,s) {load(init);}
    Matrix(P p, S s, fill_t f, T v=T(1)) : Matrix(p,s)
    {
        switch (f)
        {
            case none:
                break;
            case zero:
                fillvalue(T{0});
                break;
            case one:
                fillvalue(T{1});
                break;
            case value:
                fillvalue(v);
                break;
            case random:
                fillrandom(v); //v is max abs
                break;
            case unit:
                fillvalue(T{0});
                filldiagonal(T{1});
                break;
        }
    }
    template <isMatrix M> Matrix(P p, S s, const M& m) : Matrix(p,s) {load(m);} //assign from an expression.
    // All of the private generic versions should lead to this root constructor.
    Matrix(P p, S s) : itsPacker(p), itsShaper(s), data(itsPacker.stored_size()), itsSymmetry(data,itsPacker) {};
public:
    template <isMatrix M> auto& operator=(M&& m)
    {
        if (nr()!=m.nr() || nc()!=m.nc())
        {
            packer().resize(m.nr(),m.nc());
            shaper().resize(m.nr(),m.nc());
            data=D(packer().stored_size());
        }
        load(m);
        return *this;
    }


    //
    //  2D data access.
    //
    T  operator()(size_t i, size_t j) const
    {
        return itsSymmetry.apply(i,j);
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
            std::cout << "{";
            for (size_t j = 0; j < nc(); ++j)
            {
                std::cout << (*this)(i, j);
                if (j<nc()-1) std::cout << ",";
            }
            std::cout << "}";
            if (i<nr()-1) std::cout << "," << std::endl;;
        }
        std::cout << std::endl << std::endl;
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
                    assert(val==itsSymmetry.apply(i,j));
                ++j;
            }
            ++i;
        }
    }
    template <isMatrix M> void load(const M& m)
    {
        for (size_t i = 0; i < nr(); ++i)
            for (size_t j = 0; j < nc(); ++j)
                if (itsPacker.is_stored(i, j)) 
                    (*this)(i,j) = m(i, j);
                else
                    assert(m(i,j)==itsSymmetry.apply(i,j)); //Make sure data honours the symmetry.
    }
    void fillvalue(T v) {for (auto& i:data) i=v;}
    void fillrandom(T v) 
    {
        if (v==1)
            for (auto& i:data) i=OMLRandPos<T>();
        else
            for (auto& i:data) i=OMLRandPos<T>()*v;
    }
    void filldiagonal(T v) 
    {
        for (size_t i=0;i<nr()&&i<nc();i++)
            (*this)(i,i)=v;
    }
    P itsPacker;
    S itsShaper;
    D data;
    Sym itsSymmetry;
};



template <class T> struct FullMatrixRM : public Matrix<T,FullPackerRM,FullShaper>
{
    using Matrix<T,FullPackerRM,FullShaper>::Matrix; //Inherit base constructors.
};
template <class T> struct FullMatrixCM : public Matrix<T,FullPackerCM,FullShaper>
{
    using Matrix<T,FullPackerCM,FullShaper>::Matrix; //Inherit base constructors.
};
template <class T> struct UpperTriangularMatrixCM : public Matrix<T,UpperTriangularPackerCM,UpperTriangularShaper>
{
    using Matrix<T,UpperTriangularPackerCM,UpperTriangularShaper>::Matrix; //Inherit base constructors.
};
template <class T> struct UpperTriangularMatrixRM : public Matrix<T,UpperTriangularPackerRM,UpperTriangularShaper>
{
    using Matrix<T,UpperTriangularPackerRM,UpperTriangularShaper>::Matrix; //Inherit base constructors.
};
template <class T> struct LowerTriangularMatrixCM : public Matrix<T,LowerTriangularPackerCM,LowerTriangularShaper>
{
    using Matrix<T,LowerTriangularPackerCM,LowerTriangularShaper>::Matrix; //Inherit base constructors.
};
template <class T> struct LowerTriangularMatrixRM : public Matrix<T,LowerTriangularPackerRM,LowerTriangularShaper>
{
    using Matrix<T,LowerTriangularPackerRM,LowerTriangularShaper>::Matrix; //Inherit base constructors.
};
template <class T> struct DiagonalMatrix : public Matrix<T,DiagonalPacker,DiagonalShaper>
{
    using Matrix<T,DiagonalPacker,DiagonalShaper>::Matrix;  //Inherit base constructors.
};
template <class T> struct SBandMatrix : public Matrix<T,SBandPacker,SBandShaper>
{
public:
    using Base = Matrix<T,SBandPacker,SBandShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    SBandMatrix(                  ) : SBandMatrix(0,0) {}; //nr=nc=n=0, k=0
    SBandMatrix(size_t n, size_t k) : SBandMatrix(n,k,none) {};
    SBandMatrix(size_t n, size_t k, fill_t f, T v=T(1)) : Base(SBandPacker(n,k),f,v) {};
    SBandMatrix(const il_t& il,size_t k) : Base(SBandPacker(nr(il),k),il)
    {
        assert(nr(il)==nc(il)); //SBand only supports square matricies.
    };
    template <isMatrix M> SBandMatrix(const M& m) : Base(m.packer(),m) {}; //Only m.packer() knows the bandwidth k for the expression m

    size_t bandwidth() const {return this->packer().bandwidth();}
};

//
//  With the shaper/packer framework the definition of a Symmetric matrix is straight forward: Upper traiangular packing with full shape.
//  The complexity only arises because for non-defualt symmetry we are now forced to specify the data type.
//
template <class T> struct SymmetricMatrixCM 
: public Matrix<T,
                UpperTriangularPackerCM,
                FullShaper,
                default_data_type<T>,
                Symmetric<default_data_type<T>,UpperTriangularPackerCM>
                >
{
    using Matrix<T,
                UpperTriangularPackerCM,
                FullShaper,
                default_data_type<T>,
                Symmetric<default_data_type<T>,UpperTriangularPackerCM>
                >::Matrix;  //Inherit base constructors.
};

// Triangular with full packing. Wast of space, but that is what blas level 3 supports.
template <class T> struct UpperTriangularMatrixFCM : public Matrix<T,FullPackerCM,UpperTriangularShaper>
{
    using Matrix<T,FullPackerCM,UpperTriangularShaper>::Matrix; //Inherit base constructors.
};
template <class T> struct UpperTriangularMatrixFRM : public Matrix<T,FullPackerRM,UpperTriangularShaper>
{
    using Matrix<T,FullPackerRM,UpperTriangularShaper>::Matrix; //Inherit base constructors.
};
template <class T> struct LowerTriangularMatrixFCM : public Matrix<T,FullPackerCM,LowerTriangularShaper>
{
    using Matrix<T,FullPackerCM,LowerTriangularShaper>::Matrix; //Inherit base constructors.
};
template <class T> struct LowerTriangularMatrixFRM : public Matrix<T,FullPackerRM,LowerTriangularShaper>
{
    using Matrix<T,FullPackerRM,LowerTriangularShaper>::Matrix; //Inherit base constructors.
};






} //namespace matrix23

#include "matrix23/matops.hpp"
#include "matrix23/matmul.hpp"