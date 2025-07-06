// File: vector.hpp Define a vector class for linear algebra operations.
#pragma once

#include "matrix23/ran250.h"
#include <valarray>
// #include <vector>
#include <ranges>
#include <cassert>

namespace matrix23
{

template <class V> 
concept isVector = requires (V v)
{
    v.begin();
    v.end();
    v.indices();
    v.size();
};




template <std::ranges::viewable_range R> class VectorView
{
public:
    typedef std::ranges::iota_view<size_t,size_t> iota_view; 
    
    VectorView(R&& r, const iota_view& indices)
    : range(std::forward<R>(r)), itsIndices(indices)
    {
        assert(size() == std::ranges::size(range) && "VectorView stop index out of range");
    }
    VectorView(const R& r, const iota_view& indices)
    : range(r), itsIndices(indices)
    {
        assert(size() == std::ranges::size(range) && "VectorView stop index out of range");
    }
     VectorView(R&& r)
    : range(std::forward<R>(r)), itsIndices(size_t(0),range.size())
    {
    }
    auto begin()       { return std::ranges::begin(range); }
    auto end  ()       { return std::ranges::end  (range); }
    auto begin() const { return std::ranges::begin(range); }
    auto end  () const { return std::ranges::end  (range); }
    // many ranges don't support random access with op[].
    size_t size() const { return  itsIndices.size(); }
    iota_view indices() const { return itsIndices; }
private:
    R range; // only includes data for the non-zero portion of the vector.
    iota_view itsIndices;
};

enum fill {none, zero, one, value, random, unit};


// template <typename T> using default_data_type=std::vector<T>;
template <typename T> using default_data_type=std::valarray<T>;


template <class T, typename Data=default_data_type<T>> class Vector
{
public:
    Vector() : data(0) {};
    Vector(size_t n) : data(n) {}
    Vector(size_t n, fill f, T v=T(1)) : data(n) 
    {
        switch (f)
        {
            case none:
                break;
            case zero:
                fillvalue(T(0.0));
                break;
            case one:
            case unit:
                fillvalue(T(1.0));
                break;
            case value:
                fillvalue(v);
                break;
            case random:
                fillrandom(v); //v is max abs
                break;
                
        }
    }
    Vector(const std::initializer_list<T>& init) : data(init.size()) {assign_from(init);}
    template <std::ranges::range R> 
    Vector(const R& range) : data(range.size()) {assign_from(range);}
    template <std::ranges::range R> 
    Vector(const VectorView<R> view) : data(view.size()) 
    {
        assign_from(view);
    }
    template <isVector V> Vector& operator=(const V& v)
    {
        assign_from(v);
        return *this;
    }

    // template <std::ranges::view V> GenVector& operator=(const VectorView<V>&& view) //const won't work for views.
    // {   
    //     if (data.size() == 0)
    //         data.resize(view.size());
    //     else
    //         assert(size()==view.size());

    //     assign_from(view);
    //     return *this;
    // }

    T operator()(size_t i) const
    {
        assert(i < size());
        return data[i];
    }
    T& operator()(size_t i)
    {
        assert(i < size());
        return data[i];
    }
    size_t size() const { return data.size(); }
    auto   begin()       { return std::begin(data); }
    auto   end  ()       { return std::end  (data); }
    auto   begin() const { return std::begin(data); }
    auto   end  () const { return std::end  (data); }
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
    iota_view indices() const { return iota_view(size_t(0), size()); } //Full view of indices

    // auto view() const
    // {
    //     return VectorView<decltype(data)>(data,iota_view(size_t(0), size())); //Full view
    // }
    // auto view(const iota_view& indices) const //possibly a partial view
    // {
    //     auto v =indices | std::views::transform([this](size_t i){ return data[i]; });
    //     return VectorView<decltype(v)>(std::move(v), indices);
    // }
    // template <std::ranges::viewable_range R> operator VectorView<R>() const
    // {
    //     return VectorView<R>(data, iota_view(size_t(0), size()));
    // }

    
protected:
    template <std::ranges::range R> void assign_from(const R& range)
    {
        size_t i=0;
        for (auto r:range) data[i++] = r; //This should be where the lazy evaluation of all the chained views happens.
    }
    void fillvalue(T v) {for (auto& i:data) i=v;}
    void fillrandom(T v) 
    {
        if (v==1)
            for (auto& i:data) i=OMLRandPos<T>();
        else
            for (auto& i:data) i=OMLRandPos<T>()*v;
    }
    

    Data data;
};


template <std::ranges::range Range1,std::ranges::range Range2> auto inner_product(const Range1& a, const Range2& b)
{
    assert(a.size() == b.size() && "Ranges must be of the same size for dot product");
    std::ranges::range_value_t<Range1> dot(0); //
    for (auto [ia,ib]:std::views::zip(a,b)) dot += ia*ib;
    return dot;
} 

struct intersection
{
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
    intersection(const iota_view& a, const iota_view& b) : indices(0,0), drop1(0), drop2(0)
    {
        if (!a.empty() && !b.empty())
        {
            size_t i0=std::max(a.front(), b.front());
            size_t i1=std::min(a.back (), b.back ())+1;
            if (i0>i1) i1=i0;
            indices=std::ranges::iota_view(i0,i1); //new intersection range
            drop1 = a.front() < i0 ? i0 - a.front() : 0; // how many to drop from a
            drop2 = b.front() < i0 ? i0 - b.front() : 0; // how many to drop from b
        }
    }

    iota_view indices;
    size_t drop1,drop2;
};


auto operator*(const isVector auto& a, const isVector auto& b)
{
    intersection inter(a.indices(),b.indices());
    auto va=a | std::views::drop(inter.drop1) | std::views::take(inter.indices.size());
    auto vb=b | std::views::drop(inter.drop2) | std::views::take(inter.indices.size());
    return inner_product(va,vb);
}

auto operator+(const isVector auto& a, const isVector auto& b)
{
    assert(a.size() == b.size() && "Vectors must be of the sam  e size for addition");
    auto ab=std::views::zip_transform([](const auto& ia, const auto& ib) { return ia + ib; },a,b);
    return VectorView(std::move(ab));
}
auto operator-(const isVector auto& a, const isVector auto& b)
{
    assert(a.size() == b.size() && "Vectors must be of the sam  e size for addition");
    return VectorView(std::views::zip_transform([](const auto& ia, const auto& ib) { return ia - ib; },a,b));
}

// This obviously useful concept should be in the std:: but it is not!
template <typename T> concept arithmetic = std::integral<T> || std::floating_point<T>;

auto operator*(const isVector auto& a, const arithmetic auto& b)
{
    return VectorView(std::views::transform(a,[b](const auto& ia) { return ia*b; }));
}
auto operator*(const arithmetic auto& b,const isVector auto& a)
{
    return VectorView(std::views::transform(a,[b](const auto& ia) { return b*ia; }));
}
auto operator/(const isVector auto& a, const arithmetic auto& b)
{
    return VectorView(std::views::transform(a,[b](const auto& ia) { return ia/b; }));
}


template <isVector V> bool operator==(const V& a,const std::initializer_list<double>& b)
{
    if (a.size() != b.size()) return false;
    for (const auto& [ia,ib] : std::views::zip(a,b)) 
        if (ia != ib) return false;

    return true;
}
template <isVector Va, isVector Vb> bool operator==(const Va& a,const Vb& b)
{
    if (a.size() != b.size()) return false;
    for (const auto& [ia,ib] : std::views::zip(a,b)) 
        if (ia != ib) return false;

    return true;
}

} // namespace matrix23