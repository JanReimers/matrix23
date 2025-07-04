// File: benchmarks.cpp run some timing tests against libblas

#include "matrix23/matrix1.hpp"
#include "matrix23/blas.hpp"
#include "gtest/gtest.h"
#include <iostream>
#include <chrono>

using std::cout;
using std::endl;
using matrix23::Vector;

class Benchmarks : public ::testing::Test
{
public:
    Benchmarks() = default;
    ~Benchmarks() override = default;

    template <std::ranges::range Range> static void print(Range v)
    {
        cout << "[";
        for (auto element : v) 
            std::cout << element << " ";

        cout << "]\n";
    }
    static auto constexpr print2D=[](auto rng){for(auto r:rng)print(r);};
    typedef std::initializer_list<double> il;
    typedef std::initializer_list<il> ilil;
};

using M = matrix23::FullMatrixCM<double>;
M mymul(const M& A, const M& B)
{
    M C(A.nr(),B.nc());
    for (size_t j=0;j<B.nc();j++)
        for (size_t i=0;i<A.nr();i++)
        {
            double t=0.0;
            for (size_t k=0;k<A.nc();k++)
                t+=A(i,k)*B(k,j);
            C(i,j)=t;
        }
    return C;
}

TEST_F(Benchmarks, MatrixMultiply)
{

    
    cout << "  n    blas::gemm(ms)    matrix23::op*(ms)" << endl;
    size_t N=100;
    std::valarray<double> blas(N), m23(N);


    for (size_t i:std::ranges::iota_view(size_t(0),N))
    {
    size_t n=300;
    M A(n, n, matrix23::random);
    M B(n, n, matrix23::random);
    M C(n, n);

        auto start1 = std::chrono::high_resolution_clock::now();
        matrix23::gemm(1.0, A, B, 0.0, C);
        auto stop1 = std::chrono::high_resolution_clock::now();
        auto duration1 = duration_cast<std::chrono::microseconds>(stop1 - start1);
        // cout << "n=" << n << "  " << std::setw(10) << duration1.count() << " " ;
        blas[i]=duration1.count();
        auto start2 = std::chrono::high_resolution_clock::now();
        // M C1 = A * B;
        M C1=mymul(A,B);
        auto stop2 = std::chrono::high_resolution_clock::now();
        auto duration2 = duration_cast<std::chrono::microseconds>(stop2 - start2);
        // cout << std::setw(10) << duration2.count() << std::setw(10) << (double)duration2.count()/duration1.count() << endl;
        m23[i]=duration2.count();
    }
    cout << "blas average: " << blas.sum()/N << " m23 average: " << m23.sum()/N << endl;
}