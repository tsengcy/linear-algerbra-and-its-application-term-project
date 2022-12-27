#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/LU>
#include "tool.h"
#include <chrono>
#include <fstream>

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    int loop = atoi(argv[2]);
    // load data
    Eigen::MatrixXd A(N,N);
    Eigen::VectorXd b(N);
    std::ifstream fa("A.txt");
    for(int i = 0; i<N; i++)
    {
        for(int j = 0; j<N; j++)
        {
            fa >> A(i,j);
        }
    }
    std::ifstream fb("b.txt");
    for(int i = 0; i<N; i++)
    {
        fb >> b(i);
    }
    // gmres
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    Eigen::FullPivLU<Eigen::MatrixXd> luA(A);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() 
                << " ms\t" << std::endl;
    return 0;
}
