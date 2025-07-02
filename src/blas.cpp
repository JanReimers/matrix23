// File: blas.cpp  Interface for calling blass using matrix23 containers

#include "matrix23/blas.hpp"
extern"C" {
void dgemv_(char* trans,int* m,int* n,double* alpha, const double* A,int* lda, const double* x, int* incx, double* beta, double* y, int* incy);
void dtpmv_(char* uplo, char* trans, char* diag, int* n, const double* A, double* x,  int* incx);
void dgbmv_(char* trans,int* m,int* n,int* kl,int* ku,double* alpha,const double*A ,int* lda,const double* x,int* incx,double*	beta,double* y,int* incy);

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
template <> void gemv(double alpha, const FullMatrixRM<double>& A, const Vector<double>& x, double beta, Vector<double>& y )
{
    assert(A.nc()==x.size());
    assert(A.nr()==y.size());
    char trans='T'; //Don't transpose A.
    int m=A.nc(),n=A.nr(),inc=1;
    dgemv_(&trans,&m,&n,&alpha,A.begin(),&m,x.begin(),&inc,&beta,y.begin(),&inc);
}
template <> void tpmv(const UpperTriangularMatrixCM<double>& A, Vector<double>& x)
{
    assert(A.nc()==x.size());
    assert(A.nr()==x.size());
    char uplo='U', trans='N', diag='N'; //A is upper tri, Don't transpose A, A is not a diagonal unit.
    int n=A.nc(),inc=1;
    dtpmv_(&uplo,&trans,&diag,&n,A.begin(),x.begin(),&inc);
}
template <> void tpmv(const UpperTriangularMatrixRM<double>& A, Vector<double>& x)
{
    assert(A.nc()==x.size());
    assert(A.nr()==x.size());
    char uplo='L', trans='T', diag='N'; //Transposed and pretent lower for row major packing.
    int n=A.nc(),inc=1;
    dtpmv_(&uplo,&trans,&diag,&n,A.begin(),x.begin(),&inc);
}
template <> void tpmv(const LowerTriangularMatrixCM<double>& A, Vector<double>& x)
{
    assert(A.nc()==x.size());
    assert(A.nr()==x.size());
    char uplo='L', trans='N', diag='N'; //A is lower tri, Don't transpose A, A is not a diagonal unit.
    int n=A.nc(),inc=1;
    dtpmv_(&uplo,&trans,&diag,&n,A.begin(),x.begin(),&inc);
}
template <> void tpmv(const LowerTriangularMatrixRM<double>& A, Vector<double>& x)
{
    assert(A.nc()==x.size());
    assert(A.nr()==x.size());
    char uplo='U', trans='T', diag='N'; //Pretend A is upper tri, transpose A for row major packing.
    int n=A.nc(),inc=1;
    dtpmv_(&uplo,&trans,&diag,&n,A.begin(),x.begin(),&inc);
}
template <> void gbmv(double alpha, const SBandMatrix<double>& A, const Vector<double>& x, double beta, Vector<double>& y )
{
    assert(A.nc()==x.size());
    assert(A.nr()==y.size());
    char trans='N'; //Don't transpose A.
    int m=A.nr(),n=A.nc(),k=A.bandwidth(),lda=2*k+1,inc=1;
    dgbmv_(&trans,&m,&n,&k,&k,&alpha,A.begin(),&lda,x.begin(),&inc,&beta,y.begin(),&inc);
}

} //namespace matrix23