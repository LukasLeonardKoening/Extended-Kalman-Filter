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
    for (int i = 0; i < estimations.size(); ++i) {
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array() * residual.array();
        rsme += residual;
    }
    
    // Calculate mean
    rsme = rsme.array() / estimations.size();
    
    // Calculate root
    rsme = rsme.array().sqrt();
    
    return rsme;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
}
