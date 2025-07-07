// File: matops.hpp  Define operators other than matrix*matrix multiplication for matrices.
#pragma once

#include "matrix23/matrix1.hpp"

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
template <isMatrix Ma, isMatrix Mb> bool operator==(const Ma& a,const Mb& b)
{
    for (const auto& [ia,ib] : std::views::zip(a.rows(),b.rows())) 
        if (ia != ib) return false;
    return true;
}
 

template <isMatrix Ma, isMatrix Mb, class Op> class MatrixBinOpView
{
public:
    typedef std::remove_cvref_t<Ma>::value_type value_type;
    MatrixBinOpView(const Ma& _a, const Mb& _b, const Op& _op) : a(_a), b(_b), op(_op)
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
    auto packer() const {return a.packer();}
    auto shaper() const {return a.shaper();}
private:
    Ma a; 
    Mb b; 
    Op op; 
};

template <isMatrix Ma, class Op> class MatrixOpView
{
public:
    typedef std::remove_cvref_t<Ma>::value_type value_type;
    MatrixOpView(const Ma& _a, const Op& _op) : a(_a), op(_op) {}

    value_type operator()(size_t i, size_t j) const {return op(a(i,j));}
    size_t size() const {return nr()*nc(); }
    size_t nr  () const {return a.nr(); }
    size_t nc  () const {return a.nc(); }
    auto rows  () const {return std::views::transform(a.rows(),op);}
    auto cols  () const {return std::views::transform(a.cols(),op);}
    auto packer() const {return a.packer();}
    auto shaper() const {return a.shaper();}
private:
    Ma a; 
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




} // namespace
