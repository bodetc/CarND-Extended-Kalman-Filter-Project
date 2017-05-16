#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter():
  I_(MatrixXd::Identity(4,4)),
  F_(MatrixXd(4, 4)),
  Q_(MatrixXd(4, 4)) {
  /*****************************************************************************
   *  State
   ****************************************************************************/
  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in,
                        MatrixXd &H_in, MatrixXd &R_in) {
  x_ = x_in;
  P_ = P_in;
  H_ = H_in;
  R_ = R_in;
}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::SetPredictMatrices(float dt) {
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
    MatrixXd PHt = P_ * H_.transpose();
    VectorXd y = z - H_ * x_;
    MatrixXd S = H_ * PHt + R_;
    MatrixXd K = PHt * S.inverse();
    x_ = x_ + K*y;
    P_ = (I_ - K * H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}
