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

TEST_F(Benchmarks, MatrixMultiply)
{

    using M = matrix23::FullMatrixCM<double>;
    cout << "  n    blas::gemm(ms)    matrix23::op*(ms)" << endl;
    for (size_t n:{100,200,300,400,500,600,700,800,900,1000})
    {
    M A(n, n, matrix23::random);
    M B(n, n, matrix23::random);
    M C(n, n);

        auto start1 = std::chrono::high_resolution_clock::now();
        matrix23::gemm(1.0, A, B, 0.0, C);
        auto stop1 = std::chrono::high_resolution_clock::now();
        auto duration1 = duration_cast<std::chrono::microseconds>(stop1 - start1);
        cout << "n=" << n << "  " << std::setw(10) << duration1.count() << " " ;
        auto start2 = std::chrono::high_resolution_clock::now();
        M C1 = A * B;
        auto stop2 = std::chrono::high_resolution_clock::now();
        auto duration2 = duration_cast<std::chrono::microseconds>(stop2 - start2);
        cout << std::setw(10) << duration2.count() << std::setw(10) << (double)duration2.count()/duration1.count() << endl;
    }
}