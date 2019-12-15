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
    MatrixXd P_Ht = P_ * H_.transpose();
    
    // measurement pre-fit residual
    VectorXd y = z - H_ * x_;
    // prefit residual covariance
    MatrixXd S = H_ * P_Ht + R_;
    // optimal Kalman gain
    MatrixXd K = P_Ht * S.inverse();
    // identity matrix
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    
    // update state and state variance
    x_ = x_ + K * y;
    P_ = (I - K*H_) * P_;
}

VectorXd CartesianToPolar(const VectorXd &x_state) {
    // Calculation of polar version of given state
    float px = x_state[0];
    float py = x_state[1];
    float vx = x_state[2];
    float vy = x_state[3];
    
    float norm_pos = sqrt(pow(px,2) + pow(py,2));
    float atan = atan2(py, px);
    
    // Normalization of phi
    while (atan > M_PI) {
        atan -= 2*M_PI;
    }
    while (atan < -M_PI) {
        atan += 2*M_PI;
    }
    
    // avoid division by zero
    if (norm_pos < 0.00000001) {
        norm_pos = 0.00000001;
    }
    
    VectorXd h_x = VectorXd(3);
    h_x << norm_pos, atan, (px*vx + py*vy)/norm_pos;
    return h_x;
    
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    // Calculation of polar version of current state
    VectorXd h_x = CartesianToPolar(x_);

    // measurement pre-fit residual
    VectorXd y = z - h_x;
    // prefit residual covariance
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    // optimal Kalman gain
    MatrixXd K = P_ * H_.transpose() * S.inverse();
    // identity matrix
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

    // update state and state variance
    x_ = x_ + K * y;
    P_ = (I - K*H_) * P_;
}
