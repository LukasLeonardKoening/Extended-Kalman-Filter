#include "kalman_filter.h"
#define PI = 3.124
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;



/*
 *   Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
      //std::cout << "F=" << F_ << std::endl;
      //std::cout << "Q=" << Q_ << std::endl;
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    // measurement pre-fit residual
    VectorXd y = z - H_ * x_;
    // prefit residual covariance
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    // optimal Kalman gain
    MatrixXd K = P_ * H_.transpose() * S.inverse();
    // identity matrix
    MatrixXd I = MatrixXd::Identity(4, 4);
    
    // update state and state variance
    x_ = x_ + K * y;
    P_ = (I - K*H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    // Calculation of polar version of current state
    VectorXd h_x = VectorXd(3);
    float norm_pos = sqrt(pow(x_(0),2)+ pow(x_(1),2));
    float atan = atan2(x_(1), x_(0));
    
    // Normalization of phi
    while (atan > M_PI) {
        atan -= 2*M_PI;
    }
    while (atan < -M_PI) {
        atan += 2*M_PI;
    }
    // polar version
    h_x <<  norm_pos,
            atan,
            x_(0)*x_(2)+x_(1)*x_(3)/norm_pos;
    
      //double norm = sqrt(pow(h_x(0),2) + pow(h_x(1),2) + pow(h_x(2),2));
    //h_x = h_x.array()/norm;
      std::cout << "atan2:" << atan << std::endl;
  
    // Calculation of Jacobian matrix
    Tools t;
    MatrixXd H_j = t.CalculateJacobian(x_);
  
      
    
    // measurement pre-fit residual
    VectorXd y = z - h_x;
    // prefit residual covariance
    MatrixXd S = H_j * P_ * H_j.transpose() + R_;
    // optimal Kalman gain
    MatrixXd K = P_ * H_j.transpose() * S.inverse();
    // identity matrix
    MatrixXd I = MatrixXd::Identity(4, 4);

    // update state and state variance
    x_ = x_ + K * y;
    P_ = (I - K*H_j) * P_;
}
