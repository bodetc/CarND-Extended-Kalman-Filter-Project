#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
   * Run the whole flow of the Kalman Filter from here.
   */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);
  
  inline const Eigen::VectorXd& getX() const {
    return ekf_.x_;
  }

private:
  /**
   * Kalman Filter update and prediction math lives in here.
   */
  KalmanFilter ekf_;
  
  // check whether the tracking toolbox was initiallized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;
};

#endif /* FusionEKF_H_ */
