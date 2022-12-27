#ifndef CREATE_H
#define CREATE_H

#include <eigen3/Eigen/Core>

namespace cr
{
    Eigen::MatrixXd create(int N);
    Eigen::VectorXd rvector(int N);
    Eigen::MatrixXd orth_basis(int N);
    Eigen::MatrixXd A_data(int N);
    Eigen::VectorXd Lsolver(int N, Eigen::MatrixXd L, Eigen::VectorXd b);
    Eigen::MatrixXd rand_A(int N);
    bool HQ_maker(Eigen::MatrixXd &H, Eigen::MatrixXd &Q, Eigen::MatrixXd A);
    void QR_cal(Eigen::MatrixXd H, Eigen::MatrixXd &Q, Eigen::MatrixXd &R);
    Eigen::VectorXd Usolver(int N, Eigen::MatrixXd U, Eigen::VectorXd b);
}

#endif