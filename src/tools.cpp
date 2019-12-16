#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth) {
    VectorXd rsme = VectorXd::Zero(4);
    
    if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
        std::cout << "Error calculating RSME: Dimensions of estimation and ground_truth do not fit!" << std::endl;
        return rsme;
    }

    // Calculate squared residual
    for (unsigned int i = 0; i < estimations.size(); i++) {
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array() * residual.array();
        rsme += residual;
    }
    
    // Calculate mean
    rsme = rsme / estimations.size();
    
    // Calculate root
    rsme = rsme.array().sqrt();
    
    return rsme;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd jacobian = MatrixXd::Zero(3,4);
    double p_x = x_state(0);
    double p_y = x_state(1);
    double v_x = x_state(2);
    double v_y = x_state(3);
    
    if (p_x == 0 && p_y == 0) {
        std::cout << "Error calculating Jacobian Matrix: Division by zero!" << std::endl;
        return jacobian;
    }
    
    double pos_norm_sq = p_x*p_x + p_y*p_y;
    double pos_norm = sqrt(pos_norm_sq);
    double pos_norm_32 = pos_norm_sq * pos_norm;
    
    jacobian(0,0) = p_x / pos_norm;
    jacobian(0,1) = p_y / pos_norm;
    jacobian(1,0) = -p_y / pos_norm_sq;
    jacobian(1,1) = p_x / pos_norm_sq;
    jacobian(2,0) = p_y * (v_x*p_y - v_y*p_x) / pos_norm_32;
    jacobian(2,1) = p_x * (v_y*p_x - v_x*p_y) / pos_norm_32;
    jacobian(2,2) = p_x / pos_norm;
    jacobian(2,3) = p_y / pos_norm;
  
    return jacobian;
}
