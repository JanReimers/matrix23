// File: blas.cpp  Interface for calling blass using matrix23 containers

#include "matrix23/blas.hpp"
extern"C" {
void dgemv_(char* trans,int* m,int* n,double* alpha, const double* A,int* lda, const double* x, int* incx, double* beta, double* y, int* incy);
void dtpmv_(char* uplo, char* trans, char* diag, int* n, const double* A, double* x,  int* incx);
void dgbmv_(char* trans,int* m,int* n,int* kl,int* ku,double* alpha,const double*A ,int* lda,const double* x,int* incx,double* beta,double* y,int* incy);

void dgemm_( char* transa,char* transb,int* m,int* n,int* k,double* alpha,const double* A,int* lda,const double* B,int* ldb,double* beta, double* C,int* ldc );
// dtrmm_ would be for full packing and triangular shape.  i.e. lower zeros are stored but not refrenced.
void dtrmm_( char* side,char* uplo,char* transa,char* diag,int* m,int* n,double* alpha,const double*A,int* lda,const double*B,int* ldb );
}

namespace matrix23 {

template <> void gemv(double alpha, const FullMatrixCM<double>& A, const Vector<double>& x, double beta, Vector<double>& y )
{
    assert(A.nc()==x.size());
    assert(A.nr()==y.size());
    char trans='N'; //Don't transpose A.
    int m=A.nr(),n=A.nc(),inc=1;
    dgemv_(&trans,&m,&n,&alpha,&*A.begin(),&m,&*x.begin(),&inc,&beta,&*y.begin(),&inc);
}
template <> void gemv(double alpha, const FullMatrixRM<double>& A, const Vector<double>& x, double beta, Vector<double>& y )
{
    assert(A.nc()==x.size());
    assert(A.nr()==y.size());
    char trans='T'; //Don't transpose A.
    int m=A.nc(),n=A.nr(),inc=1;
    dgemv_(&trans,&m,&n,&alpha,&*A.begin(),&m,&*x.begin(),&inc,&beta,&*y.begin(),&inc);
}
template <> void gevm(double alpha, const FullMatrixCM<double>& A, const Vector<double>& x, double beta, Vector<double>& y )
{
    assert(A.nr()==x.size());
    assert(A.nc()==y.size());
    char trans='T'; //Do transpose A.
    int m=A.nr(),n=A.nc(),inc=1;
    dgemv_(&trans,&m,&n,&alpha,&*A.begin(),&m,&*x.begin(),&inc,&beta,&*y.begin(),&inc);
}
template <> void gevm(double alpha, const FullMatrixRM<double>& A, const Vector<double>& x, double beta, Vector<double>& y )
{
    assert(A.nr()==x.size());
    assert(A.nc()==y.size());
    char trans='N'; //Don't transpose A.
    int m=A.nc(),n=A.nr(),inc=1;
    dgemv_(&trans,&m,&n,&alpha,&*A.begin(),&m,&*x.begin(),&inc,&beta,&*y.begin(),&inc);
}


template <> void tpmv(const UpperTriangularMatrixCM<double>& A, Vector<double>& x)
{
    assert(A.nc()==x.size());
    assert(A.nr()==x.size());
    char uplo='U', trans='N', diag='N'; //A is upper tri, Don't transpose A, A is not a diagonal unit.
    int n=A.nc(),inc=1;
    dtpmv_(&uplo,&trans,&diag,&n,&*A.begin(),&*x.begin(),&inc);
}
template <> void tpmv(const UpperTriangularMatrixRM<double>& A, Vector<double>& x)
{
    assert(A.nc()==x.size());
    assert(A.nr()==x.size());
    char uplo='L', trans='T', diag='N'; //Transposed and pretent lower for row major packing.
    int n=A.nc(),inc=1;
    dtpmv_(&uplo,&trans,&diag,&n,&*A.begin(),&*x.begin(),&inc);
}
template <> void tpmv(const LowerTriangularMatrixCM<double>& A, Vector<double>& x)
{
    assert(A.nc()==x.size());
    assert(A.nr()==x.size());
    char uplo='L', trans='N', diag='N'; //A is lower tri, Don't transpose A, A is not a diagonal unit.
    int n=A.nc(),inc=1;
    dtpmv_(&uplo,&trans,&diag,&n,&*A.begin(),&*x.begin(),&inc);
}
template <> void tpmv(const LowerTriangularMatrixRM<double>& A, Vector<double>& x)
{
    assert(A.nc()==x.size());
    assert(A.nr()==x.size());
    char uplo='U', trans='T', diag='N'; //Pretend A is upper tri, transpose A for row major packing.
    int n=A.nc(),inc=1;
    dtpmv_(&uplo,&trans,&diag,&n,&*A.begin(),&*x.begin(),&inc);
}
template <> void gbmv(double alpha, const SBandMatrix<double>& A, const Vector<double>& x, double beta, Vector<double>& y )
{
    assert(A.nc()==x.size());
    assert(A.nr()==y.size());
    char trans='N'; //Don't transpose A.
    int m=A.nr(),n=A.nc(),k=A.bandwidth(),lda=2*k+1,inc=1;
    dgbmv_(&trans,&m,&n,&k,&k,&alpha,&*A.begin(),&lda,&*x.begin(),&inc,&beta,&*y.begin(),&inc);
}
template <> void gbvm(double alpha, const SBandMatrix<double>& A, const Vector<double>& x, double beta, Vector<double>& y )
{
    assert(A.nr()==x.size());
    assert(A.nc()==y.size());
    char trans='T'; //Don't transpose A.
    int m=A.nr(),n=A.nc(),k=A.bandwidth(),lda=2*k+1,inc=1;
    dgbmv_(&trans,&m,&n,&k,&k,&alpha,&*A.begin(),&lda,&*x.begin(),&inc,&beta,&*y.begin(),&inc);
}

template <> void gemm(double alpha, const FullMatrixCM<double>& A, const FullMatrixCM<double>& B, double beta, FullMatrixCM<double>& C )
{
    assert(A.nc()==B.nr());
    assert(A.nr()==C.nr());
    assert(B.nc()==C.nc());
    char transa='N', transb='N'; //Don't transpose A or B.
    int m=A.nr(),k=A.nc(),n=B.nc();
    dgemm_(&transa,&transb,&m,&n,&k,&alpha,&*A.begin(),&m,&*B.begin(),&k,&beta,&*C.begin(),&m);
}
//Don't transpose A or B.  The swap A<-->B in the dgemm call!
// Not working.
// template <> void gemm(double alpha, const FullMatrixRM<double>& B, const FullMatrixRM<double>& A, double beta, FullMatrixRM<double>& C )
// {
//     assert(B.nc()==A.nr());
//     assert(B.nr()==C.nr());
//     assert(A.nc()==C.nc());
//     char transa='T', transb='T'; 
//     int m=A.nc(),k=A.nr(),n=B.nr();
//     std::cout << "m,k,n=" << m << " " << k << " " << n << std::endl;
//     dgemm_(&transa,&transb,&m,&n,&k,&alpha,A.begin(),&k,B.begin(),&n,&beta,C.begin(),&m);
// }

template <> void trmm(double alpha, const UpperTriangularMatrixFCM<double>& A, FullMatrixCM<double>& B)
{
    assert(A.nc()==B.nr());
    assert(A.nc()==A.nr()); //A has to square.
    
    char side='L',uplo='U',transa='N', diag='N'; //A*B, A upper, don't tranpose A, A is not diagonal.
    int m=B.nr(),lda=A.nr(),n=B.nc();
    dtrmm_(&side, &uplo,&transa,&diag,&m,&n,&alpha,&*A.begin(),&lda,&*B.begin(),&m);
}
template <> void trmm(double alpha, FullMatrixCM<double>& B, const UpperTriangularMatrixFCM<double>& A)
{
    assert(A.nr()==B.nc());
    assert(A.nc()==A.nr()); //A has to square.
    
    char side='R',uplo='U',transa='N', diag='N'; //B*A, A upper, don't tranpose A, A is not diagonal.
    int m=B.nr(),lda=A.nr(),n=B.nc();
    dtrmm_(&side, &uplo,&transa,&diag,&m,&n,&alpha,&*A.begin(),&lda,&*B.begin(),&m);
}
template <> void trmm(double alpha, const LowerTriangularMatrixFCM<double>& A, FullMatrixCM<double>& B)
{
    assert(A.nc()==B.nr());
    assert(A.nc()==A.nr()); //A has to square.
    
    char side='L',uplo='L',transa='N', diag='N'; //A*B, A lower, don't tranpose A, A is not diagonal.
    int m=B.nr(),lda=A.nr(),n=B.nc();
    dtrmm_(&side, &uplo,&transa,&diag,&m,&n,&alpha,&*A.begin(),&lda,&*B.begin(),&m);
}
template <> void trmm(double alpha, FullMatrixCM<double>& B, const LowerTriangularMatrixFCM<double>& A)
{
    assert(A.nr()==B.nc());
    assert(A.nc()==A.nr()); //A has to square.
    
    char side='R',uplo='L',transa='N', diag='N'; //B*A, A upper, don't tranpose A, A is not diagonal.
    int m=B.nr(),lda=A.nr(),n=B.nc();
    dtrmm_(&side, &uplo,&transa,&diag,&m,&n,&alpha,&*A.begin(),&lda,&*B.begin(),&m);
}


} //namespace matrix23