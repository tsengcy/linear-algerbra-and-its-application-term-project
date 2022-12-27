#include <iostream>
#include <eigen3/Eigen/Core>
#include <stdlib.h>
#include <math.h>
#include <eigen3/Eigen/LU>
#include "tool.h"
#include <eigen3/Eigen/Eigenvalues>

namespace cr
{
    Eigen::MatrixXd create(int N)
    {
        Eigen::MatrixXd m(N, N);
        Eigen::VectorXd v(N);
        v = rvector(N);
        for(int i=0; i<N; i++)
        {
            for(int j = 0; j<N; j++)
            {
                m(i,j) = pow(v(i), j);
            }
        }
        return m;
    }
    Eigen::VectorXd rvector(int N)
    {
        Eigen::VectorXd v(N);
        int i = 0;
        while(i < N)
        {
            double x = (double) rand() / (RAND_MAX + 1.0) - 0.5;
            if(x > 0)
            {
                x = x + 0.5;
            }
            else
            {
                x = x - 0.5;
            }
            bool flag = true;
            for(int j = 0; j<i; j++)
            {
                if(v(j) ==  x)
                {
                    break;
                }
            }
            if(true)
            {
                v(i) = x;
                i++;
            }
        }
        return v;
    }
    Eigen::MatrixXd orth_basis(int N)
    {
        // based on gram schmit method
        Eigen::MatrixXd m(N, N);
        m = create(N);
        for(int i = 0; i<N; i++)
        {
            for(int j = 0; j<i; j++)
            {
                double a = m.col(i).transpose() * m.col(j);
                m.col(i) = m.col(i) - a * m.col(j);
            }
            m.col(i) = m.col(i) / m.col(i).norm();
        }
        
        return m;
    }
    Eigen::MatrixXd A_data(int N)
    {
        Eigen::MatrixXd Q(N,N);
        Q = orth_basis(N);
        Eigen::MatrixXd D(N,N);
        for(int i = 0; i<N; i++)
        {
            D(i,i) = (double) rand() / (RAND_MAX + 1.0) * 1 - 0.5;
        }

        return Q*D*Q.transpose();
    }
    Eigen::VectorXd Lsolver(int N, Eigen::MatrixXd L, Eigen::VectorXd b)
    {
        Eigen::VectorXd x(N);
        for(int i = 0; i<N; i++)
        {
            double a = b(i);
            for(int j = 0; j<i; j++)
            {
                a = a - L(i,j) * x(j);
            }
            x(i) = a / L(i,i);
        }
        return x;
    }
    Eigen::MatrixXd rand_A(int N)
    {
        Eigen::MatrixXd m = Eigen::MatrixXd::Random(N,N);
        return m;
    }
    bool HQ_maker(Eigen::MatrixXd &H, Eigen::MatrixXd &Q, Eigen::MatrixXd A)
    {
        int N = Q.rows();
        int k = Q.cols();
        Eigen::VectorXd nq = Q.col(Q.cols() - 1);

        H.conservativeResize(k+1, k);
        H.row(k) = Eigen::VectorXd::Zero(k);
        Q.conservativeResize(N, k+1);
        
        Eigen::VectorXd u(N);
        u = A * nq;

        double h;
        for(int i = 0; i<k; i++)
        {
            h =  Q.col(i).transpose() * u;
            H(i, k-1) = h;
            u = u - h * Q.col(i);
        }
        if(u.norm() > 0.001)
        {
            Q.col(k) = u / u.norm();
            H(k, k-1) = u.norm();
            return true;
        }
        else
        {
            Q.col(k) = Eigen::VectorXd::Zero(N);
            H(k, k-1) = 0;
            return false;
        }

    }
    void QR_cal(Eigen::MatrixXd H, Eigen::MatrixXd &Q, Eigen::MatrixXd &R)
    {
        int k = H.cols();
        Eigen::VectorXd v = H.col(k-1);
        Q.conservativeResize(k+1, k);
        Q.row(k) = Eigen::VectorXd::Zero(k);
        Q.col(k-1) = Eigen::VectorXd::Zero(k+1);
        R.conservativeResize(k, k);
        R.row(k-1) = Eigen::VectorXd::Zero(k);
        R.col(k-1) = Eigen::VectorXd::Zero(k);

        double h;
        for(int i = 0; i<k-1; i++)
        {
            h = Q.col(i).transpose() * v;
            v = v - h * Q.col(i);
        }
        Q.col(k-1) = v / v.norm();
        for(int i = 0; i<k; i++)
        {
            R(i, k-1) = Q.col(i).transpose() * H.col(k-1);
        }
    }
    Eigen::VectorXd Usolver(int N, Eigen::MatrixXd U, Eigen::VectorXd b)
    {
        Eigen::VectorXd x(N);
        for(int i = N-1; i>=0; i--)
        {
            double a = b(i);
            for(int j = N-1; j>i; j--)
            {
                a = a - U(i,j) * x(j);
            }
            x(i) = a / U(i,i);
        }
        return x;
    }

}