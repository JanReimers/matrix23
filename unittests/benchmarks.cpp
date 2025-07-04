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
    for (size_t i=0;i<A.nr();i++)
        for (size_t j=0;j<B.nc();j++)
        {
            double t=0.0;
            for (size_t k=0;k<A.nc();k++)
                t+=A(i,k)*B(k,j);
            C(i,j)=t;
        }
    return C;
}
M mymul_wcopy(const M& A, const M& B)
{
    M C(A.nr(),B.nc());
    for (size_t i=0;i<A.nr();i++)
    {
        Vector<double> Ai=A.row(i);
        for (size_t j=0;j<B.nc();j++)
        {
            double t=0.0;
            for (size_t k=0;k<A.nc();k++)
                t+=Ai(k)*B(k,j);
            C(i,j)=t;
        }
    }
    return C;
}

double average(const std::valarray<double>& a)
{
    return a.sum()/a.size();
}

double stdev(const std::valarray<double>& a)
{
    size_t Nm=a.size()-1;
    double avg=average(a);
    std::valarray<double> da=(a-avg)*(a-avg);
    return sqrt(da.sum())/Nm;
    
}

#ifndef DEBUG

TEST_F(Benchmarks, MatrixMultiply)
{
    
    cout << "  n       blas::gemm(ms)        ranges(ms)         std_mmul           mmul_wcopy" << endl;
    size_t N=10;
    const size_t iblas=0,iranges=1,imymul=2,imymul_wcopy=3;
    for ( size_t n:{100,200,300,400,500,600,700})
    {
    std::valarray<std::valarray<double>> timings(4);
    for (auto& i:timings) i=std::valarray<double>(N);


    for (size_t i:std::ranges::iota_view(size_t(0),N))
    {
       
        M A(n, n, matrix23::random);
        M B(n, n, matrix23::random);
        M C(n, n);
        {
            auto start = std::chrono::high_resolution_clock::now();
            matrix23::gemm(1.0, A, B, 0.0, C);
            auto stop = std::chrono::high_resolution_clock::now();
            timings[iblas][i]=duration_cast<std::chrono::milliseconds>(stop - start).count();
        }
        {
            auto start = std::chrono::high_resolution_clock::now();
            M C1 = A * B;
            auto stop = std::chrono::high_resolution_clock::now();
            timings[iranges][i]= duration_cast<std::chrono::milliseconds>(stop - start).count();
        }
        {
            auto start = std::chrono::high_resolution_clock::now();
            M C1=mymul(A,B);
            auto stop = std::chrono::high_resolution_clock::now();
            // cout << std::setw(10) << duration2.count() << std::setw(10) << (double)duration2.count()/duration1.count() << endl;
            timings[imymul][i]= duration_cast<std::chrono::milliseconds>(stop - start).count();
        }
        {
            auto start = std::chrono::high_resolution_clock::now();
            M C1=mymul_wcopy(A,B);
            auto stop = std::chrono::high_resolution_clock::now();
            // cout << std::setw(10) << duration2.count() << std::setw(10) << (double)duration2.count()/duration1.count() << endl;
            timings[imymul_wcopy][i]= duration_cast<std::chrono::milliseconds>(stop - start).count();
        }

    }

    cout << n << "      ";
    for (auto it:{iblas,iranges,imymul,imymul_wcopy})
    {
        double avg=average(timings[it]);
        double dev=stdev(timings[it]);
        cout << std::setprecision(1) << std::fixed << std::setw(7) << avg << "(" << std::setw(4) << dev << ")      ";
    }
    cout << endl;
    } //for n
}
#endif //DEBUG