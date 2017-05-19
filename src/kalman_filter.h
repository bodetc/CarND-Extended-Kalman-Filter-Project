#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
private:
  // Identity matrix of the correct size
  const Eigen::MatrixXd I_;
  
  // state transistion matrix
  Eigen::MatrixXd F_;
  
  // process covariance matrix
  Eigen::MatrixXd Q_;
  
  // measurement matrix
  const Eigen::MatrixXd H_laser_;
  
  // measurement covariance matrix
  const Eigen::MatrixXd R_laser_;
  
  //measurement covariance matrix - radar
  const Eigen::MatrixXd R_radar_;
  
  /**
   * Update the transition matrix F and the process covariance matrix Q to use the proper delta-t
   */
  void SetPredictionMatrices(float dt);
  
  /**
   * Performs a generic update step
   */
  void GenericUpdate(const Eigen::VectorXd &y, const Eigen::MatrixXd &H, const Eigen::MatrixXd &R);
  
  // Constant matrices initialization
  static Eigen::MatrixXd getHLaser();
  static Eigen::MatrixXd getRLaser();
  static Eigen::MatrixXd getRRadar();

public:

  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict(float dt);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

};

#endif /* KALMAN_FILTER_H_ */
