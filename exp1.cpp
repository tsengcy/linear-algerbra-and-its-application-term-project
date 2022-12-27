#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include "tool.h"
#include <chrono>
#include <fstream>

int main()
{
    int N = 3;
    int size[] = {2,5,10};
    int count[] = {0, 0, 0, 0, 0};
    int sor_count[] = {0, 0, 0, 0, 0};
        
    for (int l = 0; l<N; l++)
    {
        Eigen::MatrixXd A(size[l],size[l]);
        for(int k = 0; k<10000; k++)
        {
            A = cr::rand_A(size[l]);
            Eigen::MatrixXd L(size[l],size[l]);
            L = Eigen::ArrayXXd::Zero(size[l],size[l]);
            
            Eigen::MatrixXd U(size[l],size[l]);
            U = Eigen::ArrayXXd::Zero(size[l],size[l]);

            for(int i = 0; i<size[l]; i++)
            {
                for(int j = 0; j<size[l]; j++)
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

            

            Eigen::VectorXcd eig(size[l]);
            eig = (L.inverse()*U).eigenvalues();
    
            for(int i = 0; i<size[l]; i++)
            {
                if(pow(eig(i).real(),2) + pow(eig(i).imag(),2) > 1)
                {
                    goto there;
                }
            }
            count[l]++;
            there:
                int a;
            double omega = 1;
            Eigen::MatrixXd rL(size[l],size[l]);
            rL = Eigen::ArrayXXd::Zero(size[l],size[l]);
            Eigen::MatrixXd rU(size[l],size[l]);
            rU = Eigen::ArrayXXd::Zero(size[l],size[l]);
            Eigen::MatrixXd rD(size[l],size[l]);
            rD = Eigen::ArrayXXd::Zero(size[l],size[l]);
    
            for(int i = 0; i<size[l]; i++)
            {
                for(int j = 0; j<size[l]; j++)
                {
                    if(i > j)
                    {
                        rL(i,j) = A(i,j);
                    }
                    else if(i == j)
                    {
                        rD(i,j) = A(i,j);
                    }
                    else
                    {
                        rU(i,j) = A(i,j);
                    }
                }
            }
            Eigen::MatrixXd DL(size[l],size[l]);
            DL = rL + rD / omega;
            Eigen::MatrixXd DU(size[l],size[l]);
            DU = rD / omega - rD - rU;

            Eigen::VectorXcd reig(size[l]);
            reig = (DL.inverse()*DU).eigenvalues();
    
            for(int i = 0; i<size[l]; i++)
            {
                if(pow(reig(i).real(),2) + pow(reig(i).imag(),2) > 1)
                {
                    goto here;
                }
            }
            sor_count[l]++;
            here:
                int b;
                
        }
        std::cout << "i = " << size[l] << " available: " << count[l] << std::endl;
        std::cout << "i = " << size[l] << " available: " << sor_count[l] << std::endl;
    }
    

}