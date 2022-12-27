#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Core>

#include "tool.h"

int main(int argc, char* argv[])
{
    int N = atoi(argv[1]);
    int seed = atoi(argv[2]);
    srand(seed);
    std::ofstream file("A.txt");
    Eigen::MatrixXd m(N,N);
    m = cr::rand_A(N);

    for(int i = 0; i<N; i++)
    {
        for(int j = 0; j<N; j++)
        {
            file << m(i,j) << std::endl;
        }
    }
    file.close();

    std::ofstream f2("b.txt");
    Eigen::VectorXd x(N);
    x = cr::rvector(N);
    Eigen::VectorXd v(N);
    v = m*x;
    for(int i = 0; i<N; i++)
    {
        f2 << v(i) << std::endl;
    }
    f2.close();
}