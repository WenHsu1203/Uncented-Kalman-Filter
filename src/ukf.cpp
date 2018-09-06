#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1,0,0,0,0,
        0,1,0,0,0,
        0,0,1,0,0,
        0,0,0,1,0,
        0,0,0,0,1;
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  is_initialized_ = false;

  n_x_ = 5;

  n_aug_ = n_x_+2;

  lambda_ = 3 - n_aug_;

  time_us_ = 0;

  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int i = 1 ; i < 2*n_aug_+1; i++)
  {
    weights_(i) = 0.5/(lambda_ + n_aug_);
  }

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  NIS_laser_ = NIS_radar_ = 0.0;

  R_LIDAR_ = MatrixXd(2,2);
  R_LIDAR_ << std_laspx_*std_laspx_, 0.0,
              0.0, std_laspy_*std_laspy_;

  R_RADAR_ = MatrixXd(3,3);
  R_RADAR_ << std_radr_*std_radr_, 0.0, 0.0,
              0.0, std_radphi_*std_radphi_, 0.0,
              0.0, 0.0, std_radr_*std_radr_;
  /**
  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if(!is_initialized_)
  {
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      // double d_rho = meas_package.raw_measurements_[2];
      x_(0) = rho*cos(phi);
      x_(1) = rho*sin(phi);
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  /********************************************************************
   *** PREDICT
   *********************************************************************/

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  /********************************************************************
   *** UPDATE
   *********************************************************************/
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    UpdateRadar(meas_package);
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    UpdateLidar(meas_package); 
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // Augmentation
  VectorXd x_aug = VectorXd(7);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  MatrixXd P_aug = MatrixXd(7,7);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  // Calculate the sqrt of P_
  MatrixXd L = P_aug.llt().matrixL();
  // Find the sigma points 
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_ ; i++)
  {
    Xsig_aug.col(i+1) = x_aug+sqrt(lambda_+n_aug_)*L.col(i);
    Xsig_aug.col(i+n_aug_+1) = x_aug-sqrt(lambda_+n_aug_)*L.col(i);
  }

  for (int i = 0; i <2*n_aug_+1; i ++)
  {
    const double px = Xsig_aug(0,i);
    const double py = Xsig_aug(1,i);
    const double v  = Xsig_aug(2,i);
    const double yaw  = Xsig_aug(3,i);
    const double yawd = Xsig_aug(4,i);
    const double nu_a = Xsig_aug(5,i);
    const double nu_yawdd = Xsig_aug(6,i);

    double px_p, py_p;
    if(fabs(yawd)>0.0001)
    {
      px_p = px + v/yawd*(sin(yaw + yawd*delta_t)-sin(yaw)) + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
      py_p = py + v/yawd*(-cos(yaw + yawd*delta_t)+cos(yaw)) + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
    }
    else 
    {
      px_p = px+ v*cos(yaw)*delta_t + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
      py_p = py+ v*sin(yaw)*delta_t + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
    }

    double v_p = v+delta_t*nu_a;
    double yaw_p = yaw + yawd*delta_t + 0.5*delta_t*delta_t*nu_yawdd;
    double yawd_p = yawd + delta_t*nu_yawdd;
    Xsig_pred_(0,i)= px_p;
    Xsig_pred_(1,i)= py_p;
    Xsig_pred_(2,i)= v_p;
    Xsig_pred_(3,i)= yaw_p;
    Xsig_pred_(4,i)= yawd_p;
  }
  // Predict x_
  x_.fill(0.0);
  for (int i = 0 ; i < 2*n_aug_+1; i ++)
  {
    x_ = x_ +weights_(i)*Xsig_pred_.col(i);
  }

  // Predict P_
  P_.fill(0.0);
  for (int i = 0 ; i < 2*n_aug_+1; i ++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));
    P_ = P_ +weights_(i)*x_diff * x_diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i = 0 ; i < 2*n_aug_+1; i ++)
  {
      double p_x = Xsig_pred_(0,i);
      double p_y = Xsig_pred_(1,i);
      Zsig(0,i) = p_x;
      Zsig(1,i) = p_y;
  }
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0 ; i < 2*n_aug_+1; i ++)
  {
      z_pred = z_pred + weights_(i)*Zsig.col(i);
  }
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  MatrixXd T = MatrixXd(n_x_,n_z);
  T.fill(0.0);
  for (int i = 0 ; i < 2*n_aug_+1; i ++)
  {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      S = S + weights_(i) * z_diff * z_diff.transpose();
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      T = T + weights_(i) * x_diff * z_diff.transpose();
  }
  
  S = S + R_LIDAR_;
  MatrixXd K = T*S.inverse();
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_laser_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i = 0 ; i < 2*n_aug_+1; i ++)
  {
      double p_x = Xsig_pred_(0,i);
      double p_y = Xsig_pred_(1,i);
      double v   = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);

      Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
      Zsig(1,i) = atan2(p_y,p_x);
      Zsig(2,i) = (p_x*v*cos(yaw)+p_y*v*sin(yaw))/(sqrt(p_x*p_x+p_y*p_y));
  }
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0 ; i < 2*n_aug_+1; i ++)
  {
      z_pred = z_pred + weights_(i)*Zsig.col(i);
  }
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  MatrixXd T = MatrixXd(n_x_,n_z);
  T.fill(0.0);
  for (int i = 0 ; i < 2*n_aug_+1; i ++)
  {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));
      S = S + weights_(i) * z_diff * z_diff.transpose();
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));
      T = T + weights_(i) * x_diff * z_diff.transpose();
  }
  S = S + R_RADAR_;
  MatrixXd K = T*S.inverse();
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));
  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_radar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
}
