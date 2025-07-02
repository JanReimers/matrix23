// File: blas.hpp  Interface for calling blass using matrix23 containers
#pragma once

#include "matrix23/matrix1.hpp"
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

} //namespace matrix23