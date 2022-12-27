#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include "tool.h"
#include <chrono>
#include <fstream>

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    int loop = atoi(argv[2]);
    // data create
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

    // GS
    int iter = 1;
    double error;
    double epsilon = 0.0001;
    // split into D M = (L+U)
    Eigen::MatrixXd L(N,N);
    L = Eigen::ArrayXXd::Zero(N,N);
    Eigen::MatrixXd U(N,N);
    U = Eigen::ArrayXXd::Zero(N,N);

    for(int i = 0; i<N; i++)
    {
        for(int j = 0; j<N; j++)
        {
            if(i >= j)
            {
                L(i,j) = A(i,j);
            }
            else
            {
                U(i,j) = A(i,j);
            }
        }
    }

    Eigen::VectorXcd eig(N);
    eig = (L.inverse()*U).eigenvalues();
    
    for(int i = 0; i<N; i++)
    {
        if(pow(eig(i).real(),2) + pow(eig(i).imag(),2) > 1)
        {
            // std::cout << "unable to solve 1" << std::endl;
            return 0;
        }
    }

    std::fstream out;
    std::string str(argv[1]);
    std::string strf = "./output/";
    out.open(strf + str + "gs.csv", std::ios::app);


    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    Eigen::VectorXd x(N);
    x = cr::rvector(N);
    error = (b - A*x).norm();
    while((error > epsilon))
    {
        // update x based on the equation
        x = cr::Lsolver(N, L, ( b - U * x ));
        error = (b - A*x).norm();
        iter++;
    }
    if(true)
    {
        
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::cout << "the iteration number " << iter 
                << "\t error= " << error //<< std::endl;
                << "\t consuming time " 
                << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " ms" << std::endl;
        out << loop << "," << iter << "," << error << "," << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << std::endl;
    }
    else
    {
        // std::cout << "unable to solve" << std::endl;
        // out << loop << "," << -1 << "," << -1 << "," << -1 << std::endl;
    }
    out.close();
    
    return 0;
}






