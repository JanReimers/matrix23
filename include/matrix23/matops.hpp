// File: matops.hpp  Define operators other than matrix*matrix multiplication for matrices.
#pragma once

#include "matrix23/matrix.hpp"

namespace matrix23
{

auto operator*(const isMatrix auto& m, const isVector auto& v)
{
    assert(m.nc()==v.size());
    auto rows=m.rows();
    auto indices=std::views::iota(size_t(0),rows.size());
    auto mv= indices | std::views::transform([rows,v](size_t i) {return rows[i] * v;});
    assert(mv.size()==rows.size());
    return VectorView(std::move(mv), indices);
}
auto operator*(const isVector auto& v,const isMatrix auto& m)
{
    auto cols=m.cols();
    auto indices=std::views::iota(size_t(0),cols.size());
    auto vm=indices | std::views::transform([cols,v](size_t j) {return v * cols[j];});
    assert(vm.size()==cols.size());
    return VectorView(std::move(vm), indices);
}

template <isMatrix M> bool operator==(const M& a,const std::initializer_list<std::initializer_list<double>>& b)
{
    for (const auto& [ia,ib] : std::views::zip(a.rows(),b)) 
        if (ia != ib) return false;
    return true;
}
template <isMatrix M, isMatrix Mb> bool operator==(const M& a,const Mb& b)
{
    for (const auto& [ia,ib] : std::views::zip(a.rows(),b.rows())) 
        if (ia != ib) return false;
    return true;
}
 

template <isMatrix M, isMatrix Mb, class Op> class MatrixBinOpView
{
public:
    typedef std::remove_cvref_t<M>::value_type value_type;
    MatrixBinOpView(const M& _a, const Mb& _b, const Op& _op) : a(_a), b(_b), op(_op)
    {
        assert(a.nr()==b.nr());
        assert(a.nc()==b.nc());
    }
  
    value_type operator()(size_t i, size_t j) const {return op(a(i,j),b(i,j));}
    size_t size() const {return nr()*nc(); }
    size_t nr  () const {return a.nr(); }
    size_t nc  () const {return a.nc(); }
    auto rows  () const {return std::views::zip_transform(op,a.rows(),b.rows());}
    auto cols  () const {return std::views::zip_transform(op,a.cols(),b.cols());}
    auto packer() const {return MatrixProductPacker(a.packer(),b.packer());}
    auto shaper() const {return MatrixProductShaper(a.shaper(),b.shaper());}
private:
    M a; 
    Mb b; 
    Op op; 
};

template <isMatrix M, class Op> class MatrixOpView
{
public:
    typedef std::remove_cvref_t<M>::value_type value_type;
    MatrixOpView(const M& _a, const Op& _op) : a(_a), op(_op) {}

    value_type operator()(size_t i, size_t j) const {return op(a(i,j));}
    size_t size() const {return nr()*nc(); }
    size_t nr  () const {return a.nr(); }
    size_t nc  () const {return a.nc(); }
    auto rows  () const {return std::views::transform(a.rows(),op);}
    auto cols  () const {return std::views::transform(a.cols(),op);}
    auto packer() const {return a.packer();}
    auto shaper() const {return a.shaper();}
private:
    M a; 
    Op op; 
};



auto operator-(const isMatrix auto& a,const isMatrix auto& b)
{
    return MatrixBinOpView(a,b,[](const auto& ia, const auto& ib){return ia-ib;});                    
}
auto operator+(const isMatrix auto& a,const isMatrix auto& b)
{
    return MatrixBinOpView(a,b,[](const auto& ia, const auto& ib){return ia+ib;});                    
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

// isMatrix cannot support these because it does not store data or have 1D linear iterators..
template <typename T, isPacker P, isShaper S, typename D, isSymmetry Sym> auto& operator+=(Matrix<T,P,S,D,Sym>& a, const arithmetic auto& b)
{
    for (auto& ia:a) ia+=b;
    return a;
}
template <typename T, isPacker P, isShaper S, typename D, isSymmetry Sym> auto& operator-=(Matrix<T,P,S,D,Sym>& a, const arithmetic auto& b)
{
    for (auto& ia:a) ia-=b;
    return a;
}
template <typename T, isPacker P, isShaper S, typename D, isSymmetry Sym> auto& operator*=(Matrix<T,P,S,D,Sym>& a, const arithmetic auto& b)
{
    for (auto& ia:a) ia*=b;
    return a;
}
template <typename T, isPacker P, isShaper S, typename D, isSymmetry Sym> auto& operator/=(Matrix<T,P,S,D,Sym>& a, const arithmetic auto& b)
{
    for (auto& ia:a) ia/=b;
    return a;
}


template <typename Ta, typename Tb,isPacker P, isShaper S, typename D, isSymmetry Sym> 
auto& operator+=(Matrix<Ta,P,S,D,Sym>& a, const Matrix<Tb,P,S,D,Sym>& b)
{
    auto ib=b.begin();
    for (auto& ia:a) ia+=*ib++;
    return a;
}
template <typename Ta, typename Tb,isPacker P, isShaper S, typename D, isSymmetry Sym> 
auto& operator-=(Matrix<Ta,P,S,D,Sym>& a, const Matrix<Tb,P,S,D,Sym>& b)
{
    auto ib=b.begin();
    for (auto& ia:a) ia+=*ib++;
    return a;
}

template <typename T,isPacker P, isShaper S, typename D, isSymmetry Sym> 
auto& operator+=(Matrix<T,P,S,D,Sym>& a, const isMatrix auto& b)
{
    assert(a.nr()==b.nr());
    assert(a.nc()==b.nc());
    for (auto i:iota_view(size_t(0),a.nr()))
        for (auto j:a.shaper().nonzero_col_indexes(i))
            if (a.packer().is_stored(i,j)) a(i,j)+=b(i,j);
    return a;
}
template <typename T,isPacker P, isShaper S, typename D, isSymmetry Sym> 
auto& operator-=(Matrix<T,P,S,D,Sym>& a, const isMatrix auto& b)
{
    assert(a.nr()==b.nr());
    assert(a.nc()==b.nc());
    for (auto i:iota_view(size_t(0),a.nr()))
        for (auto j:a.shaper().nonzero_col_indexes(i))
            if (a.packer().is_stored(i,j)) a(i,j)-=b(i,j);
    return a;
}

template <isMatrix M> class MatrixTransposeView
{
public:
    typedef std::remove_cvref_t<M>::value_type value_type;
    MatrixTransposeView(const M& _m) : m(_m) {}

    value_type operator()(size_t i, size_t j) const {return m(j,i);}
    size_t size() const {return nr()*nc(); }
    size_t nr  () const {return m.nc(); }
    size_t nc  () const {return m.nr(); }
    auto rows  () const {return m.cols();}
    auto cols  () const {return m.rows();}
    auto packer() const {return m.packer().transpose();}
    auto shaper() const {return m.shaper().transpose();}
private:
    M m; 
};

auto Transpose(isMatrix auto& m) {return MatrixTransposeView(m);}
auto operator~(isMatrix auto& m) {return MatrixTransposeView(m);}

template <isMatrix M> auto fnorm(const M& m)
{
    typename M::value_type fn{0};
    for (const auto& row:m.rows())
        for (const auto& ic:row)
            fn+=ic*ic;
    return sqrt(fn);
}

} // namespace
