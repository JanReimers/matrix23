// File: blas.hpp  Interface for calling blass using matrix23 containers
#pragma once

#include "matrix23/matrix.hpp"
namespace matrix23 {

// https://www.netlib.org/lapack/explore-html/d7/dda/group__gemv_ga4ac1b675072d18f902db8a310784d802.html#ga4ac1b675072d18f902db8a310784d802
// DGEMV  performs one of the matrix-vector operations
//
//    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
//
// where alpha and beta are scalars, x and y are vectors and A is an
// m by n matrix.
//     
template <class T> void gemv(T alpha, const FullMatrixCM<T>& A, const Vector<T>& x, T beta, Vector<T>& y );
template <class T> void gevm(T alpha, const FullMatrixCM<T>& A, const Vector<T>& x, T beta, Vector<T>& y );
template <class T> void gemv(T alpha, const FullMatrixRM<T>& A, const Vector<T>& x, T beta, Vector<T>& y );
template <class T> void gevm(T alpha, const FullMatrixRM<T>& A, const Vector<T>& x, T beta, Vector<T>& y );
template <class T> void gbmv(T alpha, const  SBandMatrix<T>& A, const Vector<T>& x, T beta, Vector<T>& y );
template <class T> void gbvm(T alpha, const  SBandMatrix<T>& A, const Vector<T>& x, T beta, Vector<T>& y );

template <class T> void tpmv(const UpperTriangularMatrixCM<T>& A, Vector<T>& x);
template <class T> void tpmv(const UpperTriangularMatrixRM<T>& A, Vector<T>& x);
template <class T> void tpmv(const LowerTriangularMatrixCM<T>& A, Vector<T>& x);
template <class T> void tpmv(const LowerTriangularMatrixRM<T>& A, Vector<T>& x);



template <class T> void gemm(T alpha, const FullMatrixCM<T>& A, const FullMatrixCM<T>& B, T beta, FullMatrixCM<T>& C);
// template <class T> void gemm(T alpha, const FullMatrixRM<T>& A, const FullMatrixRM<T>& B, T beta, FullMatrixRM<T>& C);

//
//  Convenience helper functions so users don't need to worry about alpha.beta and constructing the return container.
//
template <class T, isMatrix Mat> void gemv(const Mat& A, const Vector<T>& x,  Vector<T>& y) {gemv(T(1),A,x,T(0),y);}
template <class T, isMatrix Mat> void gevm(const Mat& A, const Vector<T>& x,  Vector<T>& y) {gevm(T(1),A,x,T(0),y);}
template <class T> void gemv(const SBandMatrix<T>& A, const Vector<T>& x,  Vector<T>& y) {gbmv(T(1),A,x,T(0),y);}
template <class T> void gevm(const SBandMatrix<T>& A, const Vector<T>& x,  Vector<T>& y) {gbvm(T(1),A,x,T(0),y);}

template <class T, isMatrix Mat> Vector<T> blasmv(const Mat& M, const Vector<T>& v)
{
    assert(M.nc()==v.size());
    Vector<T> Mv(M.nr());
    gemv(M,v,Mv);
    return Mv;
}
template <class T, isMatrix Mat> Vector<T> blasvm(const Vector<T>& v,const Mat& M)
{
    assert(M.nr()==v.size());
    Vector<T> vM(M.nc());
    gevm(M,v,vM);
    return vM;
}



template <class T> FullMatrixCM<T> blasmm(const FullMatrixCM<T>& A, const FullMatrixCM<T>& B)
{
    assert(A.nc()==B.nr());
    FullMatrixCM<T> C(A.nr(),B.nc());
    gemm(1.0,A,B,0.0,C);
    return C;
}
// template <class T> FullMatrixRM<T> blasmm(const FullMatrixRM<T>& A, const FullMatrixRM<T>& B)
// {
//     assert(A.nc()==B.nr());
//     FullMatrixRM<T> C(A.nr(),B.nc());
//     gemm(1.0,A,B,0.0,C);
//     return C;
// }


} //namespace matrix23