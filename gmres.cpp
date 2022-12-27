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

    Eigen::FullPivLU<Eigen::MatrixXd> luA(A);
    if(luA.rank() != N)
    {
        std::cout << "unable to solve 1" << std::endl;
        return 0;
    }
    // gmres
    int iter = 1;
    double error;
    double epsilon = 0.0001;

    std::fstream out;
    out.open("gmres.txt", std::ios::app);

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    Eigen::VectorXd x = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd xn(N);
    // x = cr::rvector(N);

    error = (b - A*x).norm();
    // std::cout << error << std::endl;
    bool flags;
    
    while(error > epsilon)
    {
        Eigen::MatrixXd H(0,0);
        Eigen::MatrixXd Q(N,1);
        Eigen::MatrixXd R(0,0);
        Eigen::MatrixXd O(0,0);

        Eigen::VectorXd v = b - A * x;
        double beta = v.norm();

        Q.col(0) = v / beta;
        Eigen::VectorXd y;
        int i = 0;
        for(i; i<N; i++)
        {
            flags = cr::HQ_maker(H, Q, A);
            // if(flags == false)
            // {
            //     // std::cout << i << std::endl;
            //     // break;
            // }
            cr::QR_cal(H, O, R); 
            y = cr::Usolver(i+1, R, beta * O.row(0));
            xn = x + Q.block(0, 0, Q.rows(), i+1) * y;
            iter++;
            // x = xn;
            error = (b - A*xn).norm();
            if(error < epsilon) goto there;
            // std::cout << i << "\t" << error << std::endl;  
        }
        x = xn;
        // y = cr::Usolver(i, R, beta * O.row(0));
        // xn = x + Q.block(0, 0, Q.rows(), i) * y;
        // iter++;
        // x = xn;
        // error = (b - A*x).norm();
        // std::cout << error << std::endl;
    }
    
    there:
    x = xn;
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::cout   << "the iteration number " << iter 
                << "\t error= " << error //<< std::endl;
                << "\t consuming time " 
                << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " ms\t" << std::endl;
    out << loop << "\t" << iter << "\t" << error << "\t" 
        << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << std::endl;
    // out.close();
    
    return 0;
}







