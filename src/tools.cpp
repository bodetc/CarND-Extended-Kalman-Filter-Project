#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size()!=ground_truth.size() || estimations.size()==0) {
    std::cerr << "Error! The input is not valid for calculating the RMSE";
    return rmse;
  }
  
  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd diff = estimations[i] - ground_truth[i];
    VectorXd residual = diff.array() * diff.array();
    rmse += residual;
  }
  
  //calculate the mean
  rmse = rmse/estimations.size();
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  
  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  // get the distance from the car
  float r2 = px*px + py*py;
  float r = sqrt(r2);
  
  //check division by zero
  if(r2<0.0000001) {
    std::cerr << "Error! Impossible to divide by zero";
    return Hj;
  }
  
  // position on the unit circle
  float pxr = px/r;
  float pyr = py/r;
  
  // other variables needed in the Jacobian
  float diff = (vx*py-vy*px);
  float a = pyr*diff/r2;
  float b = -pxr*diff/r2;
  
  //compute the Jacobian matrix
  Hj << pxr, pyr, 0., 0.,
  -py/r2, px/r2, 0., 0.,
  a, b, pxr, pyr;
  
  return Hj;
}
