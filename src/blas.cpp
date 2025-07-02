// File: blas.cpp  Interface for calling blass using matrix23 containers

#include "matrix23/blas.hpp"
extern"C" {
void dgemv_(char* trans,int* m,int* n,double* alpha, const double* A,int* lda, const double* x, int* incx, double* beta, double* y, int* incy);
}

namespace matrix23 {

template <> void gemv(double alpha, const FullMatrixCM<double>& A, const Vector<double>& x, double beta, Vector<double>& y )
{
    assert(A.nc()==x.size());
    assert(A.nr()==y.size());
    char trans='N'; //Don't transpose A.
    int m=A.nr(),n=A.nc(),inc=1;
    dgemv_(&trans,&m,&n,&alpha,A.begin(),&m,x.begin(),&inc,&beta,y.begin(),&inc);
}

} //namespace matrix23