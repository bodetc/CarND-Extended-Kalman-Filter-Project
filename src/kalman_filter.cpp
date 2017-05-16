#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter():
  I_(MatrixXd::Identity(4,4)),
  x_(VectorXd(4)),
  P_(MatrixXd(4, 4)),
  F_(MatrixXd(4, 4)),
  Q_(MatrixXd(4, 4)),
  H_laser_(MatrixXd(2, 4)),
  R_laser_(MatrixXd(2, 2)) {
    
  /*****************************************************************************
   *  State
   ****************************************************************************/
  
  // x_in Initial state
  x_ << 1, 1, 1, 1;
  // Initial state covariance
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  
  // state transistion matrix
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  // Initial process covariance matrix Q, will be overriden later
  Q_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        1, 0, 1, 0,
        0, 1, 0, 1;
    
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  //measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict(float dt) {
  SetPredictionMatrices(dt);
  
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::SetPredictionMatrices(float dt) {
  //Modify the F matrix so that the time is integrated
  F_(0,2) = dt;
  F_(1,3) = dt;
  
  // Some often used variables
  float dt2 = pow(dt,2);
  float dt3 = .5*pow(dt,3);
  float dt4 = .25*pow(dt,4);
  
  //Acceleration noise
  float noise_ax=9.;
  float noise_ay=9.;
  
  //Set the process covariance matrix Q
  Q_ << dt4*noise_ax, 0.,           dt3*noise_ax, 0.,
        0.,           dt4*noise_ay, 0.,           dt3*noise_ay,
        dt3*noise_ax, 0.,           dt2*noise_ax, 0.,
        0.,           dt3*noise_ay, 0.,           dt2*noise_ay;
}

void KalmanFilter::Update(const VectorXd &z) {
  MatrixXd PHt = P_ * H_laser_.transpose();
  VectorXd y = z - H_laser_ * x_;
  MatrixXd S = H_laser_ * PHt + R_laser_;
  MatrixXd K = PHt * S.inverse();
  x_ = x_ + K*y;
  P_ = (I_ - K * H_laser_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}
