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

    Matrix(P p, S s) : itsPacker(p), itsShaper(s), data(itsPacker.stored_size()), itsSymmetry(data,itsPacker) {};
    Matrix(P p) : Matrix(p,p.shaper()) {};

    Matrix(const il_t& init,P p) : Matrix(p) {load(init);}
    Matrix(const il_t& init,P p, S s) : Matrix(p,s) {load(init);}
    template <isMatrix M> Matrix(const M& m,P p, S s) : Matrix(p,s) {load(m);}
    template <isMatrix M, isPacker otherP> Matrix(const M& m,otherP p, S s)  : Matrix(P(p.nr(),p.nc()),s) {load(m);}
    Matrix(P p, fill f, T v=T(1)) : Matrix(p)
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
    Matrix(P p, S s, fill f, T v=T(1)) : Matrix(p,s)
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
            if (i<nr()-1) std::cout << ",";
        }
        std::cout << std::endl;
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



template <class T> class FullMatrixRM : public Matrix<T,FullPackerRM,FullShaper>
{
public:
    using Base = Matrix<T,FullPackerRM,FullShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    FullMatrixRM(size_t nr, size_t nc) : Base(FullPackerRM(nr,nc)) {};
    FullMatrixRM(size_t nr, size_t nc, fill f, T v=T(1)) : Base(FullPackerRM(nr,nc),f,v) {};
    FullMatrixRM(const il_t& il) : Base(il,FullPackerRM(nr(il),nc(il))) {};
    template <isMatrix M> FullMatrixRM(const M& m) : Base(m,m.packer(), m.shaper()) {};
    
};
template <class T> class FullMatrixCM : public Matrix<T,FullPackerCM,FullShaper>
{
public:
    using Base = Matrix<T,FullPackerCM,FullShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    FullMatrixCM(size_t nr, size_t nc) : Base(FullPackerCM(nr,nc))  {};
    FullMatrixCM(size_t nr, size_t nc, fill f, T v=T(1)) : Base(FullPackerCM(nr,nc),f,v)  {};
    FullMatrixCM(const il_t& il) : Base(il,FullPackerCM(nr(il),nc(il)))  {};
    template <isMatrix M> FullMatrixCM(const M& m) : Base(m,m.packer(), m.shaper()) {};
};

template <class T> class UpperTriangularMatrixCM : public Matrix<T,UpperTriangularPackerCM,UpperTriangularShaper>
{
public:
    using Base = Matrix<T,UpperTriangularPackerCM,UpperTriangularShaper>;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    UpperTriangularMatrixCM(size_t nr, size_t nc) : Base(UpperTriangularPackerCM(nr,nc)) {};
    UpperTriangularMatrixCM(size_t nr, size_t nc, fill f, T v=T(1)) : Base(UpperTriangularPackerCM(nr,nc),f,v) {};
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
    UpperTriangularMatrixRM(size_t nr, size_t nc, fill f, T v=T(1)) : Base(UpperTriangularPackerRM(nr,nc),f,v) {};
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
    LowerTriangularMatrixCM(size_t nr, size_t nc, fill f, T v=T(1)) : Base(LowerTriangularPackerCM(nr,nc),f,v) {};
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
    LowerTriangularMatrixRM(size_t nr, size_t nc, fill f, T v=T(1)) : Base(LowerTriangularPackerRM(nr,nc),f,v) {};
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
    SBandMatrix(size_t n, size_t k) : Base(SBandPacker(n,k)) {};
    SBandMatrix(size_t n, size_t k, fill f, T v=T(1)) : Base(SBandPacker(n,k),f,v) {};
    SBandMatrix(const il_t& il,size_t k) : Base(il,SBandPacker(nr(il),k))
    {
        assert(nr(il)==nc(il)); //SBand only supports square matricies.
    };
    template <isMatrix M> SBandMatrix(const M& m) : Base(m,m.packer(), m.shaper()) {};

    size_t bandwidth() const {return this->packer().bandwidth();}
};

//
//  With the shaper/packer framework the definition of a Symmetric matrix is straight forward: Upper traiangular packing with full shape.
//  The complexity only arises because for non-defualt symmetry we are now forced to specify the data type.
//
template <class T> class SymmetricMatrixCM 
: public Matrix<T,
                UpperTriangularPackerCM,
                FullShaper,
                default_data_type<T>,
                Symmetric<default_data_type<T>,UpperTriangularPackerCM>
                >
{
public:
    using Base = Matrix<T,
                UpperTriangularPackerCM,
                FullShaper,
                default_data_type<T>,
                Symmetric<default_data_type<T>,UpperTriangularPackerCM>
                >;
    using il_t=Base::il_t;
    using Base::nr;
    using Base::nc;
    SymmetricMatrixCM(size_t n) : Base(UpperTriangularPackerCM(n,n),FullShaper(n,n)) {};
    SymmetricMatrixCM(size_t n, fill f, T v=T{1}) : Base(UpperTriangularPackerCM(n,n),FullShaper(n,n),f,v) {};
    SymmetricMatrixCM(const il_t& il) : Base(il,UpperTriangularPackerCM(nr(il),nc(il)),FullShaper(nr(il),nc(il))) {assert(nr()==nc());};
    template <isMatrix M> SymmetricMatrixCM(const M& m) : Base(m,m.packer(), m.shaper()) {};
};







} //namespace matrix23

#include "matrix23/matops.hpp"
#include "matrix23/matmul.hpp"