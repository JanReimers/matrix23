// File: Subscriptors.hpp  Define dense packing and indexing for full,tri,tridiagona,diagonal matrices.
#pragma once

#include "matrix23/packer.hpp"

namespace matrix23
{
//
//  Define what matrix shapes result from multiply two matricies.
//  These rules are identical the Packers and Shapers not sure if there is a clever
//  way to define the rules once for both Packer and Shaper?

template <isPacker P1, isPacker P2> struct MatrixProductPackerType;

template <isPacker P> struct MatrixProductPackerType<FullPackerCM,P> {typedef FullPackerCM packer_t;};
template <isPacker P> struct MatrixProductPackerType<P,FullPackerCM> {typedef FullPackerCM packer_t;};
template <isPacker P> struct MatrixProductPackerType<FullPackerRM,P> {typedef FullPackerRM packer_t;};
template <isPacker P> struct MatrixProductPackerType<P,FullPackerRM> {typedef FullPackerRM packer_t;};
template <isPacker P> struct MatrixProductPackerType<P,DiagonalPacker> {typedef P packer_t;};
template <isPacker P> struct MatrixProductPackerType<DiagonalPacker,P> {typedef P packer_t;};
template <> struct MatrixProductPackerType<FullPackerCM,FullPackerCM> {typedef FullPackerCM packer_t;};
template <> struct MatrixProductPackerType<FullPackerRM,FullPackerRM> {typedef FullPackerRM packer_t;};
template <> struct MatrixProductPackerType<DiagonalPacker,FullPackerCM> {typedef FullPackerCM packer_t;};
template <> struct MatrixProductPackerType<FullPackerCM,DiagonalPacker> {typedef FullPackerCM packer_t;};
template <> struct MatrixProductPackerType<DiagonalPacker,FullPackerRM> {typedef FullPackerRM packer_t;};
template <> struct MatrixProductPackerType<FullPackerRM,DiagonalPacker> {typedef FullPackerRM packer_t;};
template <> struct MatrixProductPackerType<DiagonalPacker,DiagonalPacker> {typedef DiagonalPacker packer_t;};
template <> struct MatrixProductPackerType<UpperTriangularPackerCM,UpperTriangularPackerCM> {typedef UpperTriangularPackerCM packer_t;};
template <> struct MatrixProductPackerType<LowerTriangularPackerCM,LowerTriangularPackerCM> {typedef LowerTriangularPackerCM packer_t;};
template <> struct MatrixProductPackerType<UpperTriangularPackerCM,LowerTriangularPackerCM> {typedef FullPackerCM packer_t;};
template <> struct MatrixProductPackerType<LowerTriangularPackerCM,UpperTriangularPackerCM> {typedef FullPackerCM packer_t;};
template <> struct MatrixProductPackerType<SBandPacker,SBandPacker> {typedef SBandPacker packer_t;}; //Need to add the ks somehow.



template <isShaper P1, isShaper P2> struct MatrixProductShaperType;
template <isShaper S> struct MatrixProductShaperType<FullShaper,S> {typedef FullShaper shaper_t;};
template <isShaper S> struct MatrixProductShaperType<S,FullShaper> {typedef FullShaper shaper_t;};
template <isShaper S> struct MatrixProductShaperType<S,DiagonalShaper> {typedef S shaper_t;};
template <isShaper S> struct MatrixProductShaperType<DiagonalShaper,S> {typedef S shaper_t;};
template <> struct MatrixProductShaperType<FullShaper,FullShaper> {typedef FullShaper shaper_t;};
template <> struct MatrixProductShaperType<DiagonalShaper,FullShaper> {typedef FullShaper shaper_t;};
template <> struct MatrixProductShaperType<FullShaper,DiagonalShaper> {typedef FullShaper shaper_t;};
template <> struct MatrixProductShaperType<DiagonalShaper,DiagonalShaper> {typedef DiagonalShaper shaper_t;};
template <> struct MatrixProductShaperType<UpperTriangularShaper,UpperTriangularShaper> {typedef UpperTriangularShaper shaper_t;};
template <> struct MatrixProductShaperType<LowerTriangularShaper,LowerTriangularShaper> {typedef LowerTriangularShaper shaper_t;};
template <> struct MatrixProductShaperType<UpperTriangularShaper,LowerTriangularShaper> {typedef FullShaper shaper_t;};
template <> struct MatrixProductShaperType<LowerTriangularShaper,UpperTriangularShaper> {typedef FullShaper shaper_t;};
template <> struct MatrixProductShaperType<SBandShaper,SBandShaper> {typedef SBandShaper shaper_t;}; //Need to add the ks somehow.




} // namespace matrix23