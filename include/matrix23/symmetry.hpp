// File: symmetry.hpp  Define the symmetry concept and related classes.
#pragma once

#include "matrix23/packer.hpp"
//
// Symmetry decides how matrix elements are related on constrained.  For example for a symetric matrix A(i,j)=A(j,i). 
// Symmetry ties in with packing when considereing storage efficiency.  In particular upper or lower triangula packing
// is an efficient means of storing a symmetric matrix.  Symmetry also require a shape that supports the
// symmetry.  For example a symmetric matrix can have full, sband or diagonal shape, but not triangular shape.  triangle
// shape is inherently non symmetric.
//
namespace matrix23
{
// c++20 concept definition for matrix symmetry.    
template <class S> concept isSymmetry = requires (S const s,size_t i)
{
    s.apply(i,i);
};

template <class D, isPacker P> struct NoSymmetry
{
    NoSymmetry(const D& d, const P& p) : data(d), packer(p) {};
    D::value_type  apply(size_t i, size_t j) const {return packer.is_stored(i,j) ? data[packer.offset(i,j)] : typename D::value_type(0);}
private:
    const D& data;
    const P& packer;
};
template <class D, isPacker P> struct Symmetric
{
    Symmetric(D& d, const P& p) : data(d), packer(p) {};
    D::value_type  apply(size_t i, size_t j) const 
    {
        bool storedij=packer.is_stored(i,j);
        assert(storedij || packer.is_stored(j,i));
        return  storedij ? data[packer.offset(i,j)] : data[packer.offset(j,i)];
    }
private:
    const D& data;
    const P& packer;
};
template <class D, isPacker P> struct AntiSymmetric
{
    AntiSymmetric(const D& d, const P& p) : data(d), packer(p) {};
    D::value_type  apply(size_t i, size_t j) const 
    {
        bool storedij=packer.is_stored(i,j);
        assert(storedij || packer.is_stored(j,i));
        return  storedij ? data[packer.offset(i,j)] : -data[packer.offset(j,i)];
    }
private:
    const D& data;
    const P& packer;
};

// static_assert(isSymmetry<   NoSymmetry<default_data_type<double>,FullPackerCM>> );
// static_assert(isSymmetry<    Symmetric<default_data_type<double>,FullPackerCM>> );
// static_assert(isSymmetry<AntiSymmetric<default_data_type<double>,FullPackerCM>> );


} //namespace
