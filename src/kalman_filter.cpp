#include "kalman_filter.h"
#include "tools.h"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter():
  I_(MatrixXd::Identity(4,4)),
  x_(VectorXd(4)),
  P_(MatrixXd(4, 4)),
  F_(MatrixXd(4, 4)),
  Q_(MatrixXd(4, 4)),
  H_laser_(getHLaser()),
  R_laser_(getRLaser()),
  R_radar_(getRRadar()){
    
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
}

KalmanFilter::~KalmanFilter() {}

MatrixXd KalmanFilter::getHLaser() {
  MatrixXd H_laser = MatrixXd(2, 4);
  H_laser << 1, 0, 0, 0,
             0, 1, 0, 0;
  return H_laser;
}

MatrixXd KalmanFilter::getRLaser() {
  MatrixXd R_laser = MatrixXd(2, 2);
  R_laser << 0.0225, 0,
             0, 0.0225;
  return R_laser;
}

MatrixXd KalmanFilter::getRRadar() {
  MatrixXd R_radar = MatrixXd(3, 3);
  R_radar << 0.09, 0, 0,
             0, 0.0009, 0,
             0, 0, 0.09;
  return R_radar;
}

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

void KalmanFilter::GenericUpdate(const VectorXd &y, const MatrixXd &H, const MatrixXd &R) {
  MatrixXd PHt = P_ * H.transpose();
  MatrixXd S = R + H * PHt;
  MatrixXd K = PHt * S.inverse();
  x_ = x_ + K * y;
  P_ = (I_ - K * H) * P_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_laser_ * x_;
  GenericUpdate(y, H_laser_, R_laser_);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd y = z - Tools::CalculateMeasurementFunction(x_);
  y(1)=remainderf(y(1), 2*M_PI); // normalization
  MatrixXd Hj = Tools::CalculateJacobian(x_);
  GenericUpdate(y, Hj, R_radar_);
}
