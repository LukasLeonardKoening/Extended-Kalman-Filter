#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 1477010443050000; // war mal 0

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
              0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

    // covariance of process noise
    ekf_.Q_ = MatrixXd(4,4);
    
    // state transition model
    ekf_.F_ = MatrixXd(4,4);
    ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

    // observation model
    H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

VectorXd maxRSME = VectorXd::Zero(4);

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /**
    * Initialization
    */
    if (!is_initialized_) {
        /**
         * Initialize the state ekf_.x_ with the first measurement.
         */

        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;
        VectorXd z = measurement_pack.raw_measurements_;
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            ekf_.x_ << z(0) * cos(z(1)), z(0) * sin(z(1)), 0, 0;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
        }
      
        ekf_.P_ = MatrixXd(4,4);
        ekf_.P_ <<     1,0,0,0,
                      0,1,0,0,
                      0,0,100,0,
                      0,0,0,100;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

  /**
   * Prediction
   */
    
    // Calculate elapsed time since last measurement and implement them into the covariance matricies
    double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
      
    ekf_.F_ <<  1, 0, dt, 0,
                0, 1, 0, dt,
                0, 0, 1, 0,
                0, 0, 0, 1;
      double dt_2 = dt * dt;
      double dt_3 = dt_2 * dt;
    double dt_4 = dt_3 * dt;
    int noise_ax = 9;
    int noise_ay = 9;
    
    double q_11 = dt_4/4 * noise_ax;
    double q_13 = dt_3/2 * noise_ax;
    double q_22 = dt_4/4 * noise_ay;
    double q_24 = dt_3/2 * noise_ay;
    double q_33 = dt_2 * noise_ax;
    double q_44 = dt_2 * noise_ay;
    
      
    ekf_.Q_ << q_11, 0, q_13, 0,
              0, q_22, 0, q_24,
              q_13, 0, q_33, 0,
              0, q_24, 0, q_44;

    ekf_.Predict();

  /**
   * Update
   */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        // Laser updates
        ekf_.R_ = R_laser_;
        ekf_.H_ = H_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
