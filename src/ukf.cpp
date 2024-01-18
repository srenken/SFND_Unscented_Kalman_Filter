#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  /* Vector Elements: 
    p_x 
    p_y
    v
    yaw (psi)
    yawd (psi_dot)
  */
  x_ = VectorXd(5); 

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2.5;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  // Until first measurement isn't received / processed KF isn't initializied
  is_initialized_ = false;

  // Set dimension
  n_x_ = x_.size();

  // Augmented dimension (add process noise to mean state vector)
  n_aug_ = n_x_ + 2;  

  // Define Lamda
  lambda_ = 3 - n_x_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  // Set Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  // Weight for i = 0
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  for (int i = 1; i < 2 * n_aug_ + 1; ++i)
  {
    weights_(i) = 1.0 / (2.0 * (lambda_ + n_aug_) );
  }

  // init timestamp
  time_us_ = 0.0;


}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  std::cout << "x_" << std::endl << x_ << std::endl;
  std::cout << "P_" << std::endl << P_ << std::endl;
  
  // First measurement?
  if(!is_initialized_) 
  {
    
    // State
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) 
    {
      double radarmeas_radial_dist  = meas_package.raw_measurements_[0];
      double radarmeas_phi          = meas_package.raw_measurements_[1];
      double radarmeas_radial_v     = meas_package.raw_measurements_[2];

      double p_x = radarmeas_radial_dist * cos(radarmeas_phi);
      double p_y = radarmeas_radial_dist * sin(radarmeas_phi);
      double v_x = radarmeas_radial_v * cos(radarmeas_phi);
      double v_y = radarmeas_radial_v * sin(radarmeas_phi);

      double v_vehicle = sqrt(v_x * v_x + v_y * v_y);

      
      x_ << 
          p_x,                // p_x
          p_y,                // p_y
          v_vehicle,          // v
          0,                 // vehicle yaw (psi)
          0;                 // vehicle yaw accelleration (psi_dot)
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) 
    {      
      x_ << 
          meas_package.raw_measurements_[0], 
          meas_package.raw_measurements_[1], 
          0, 
          0, 
          0;
    }

    // Covariance
    // P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
    //       0, std_laspy_ * std_laspy_, 0, 0, 0,
    //       0, 0, std_radphi_ * std_radphi_, 0, 0,
    //       0, 0, 0, std_radr_ * std_radr_, 0,
    //       0, 0, 0, 0, std_radrd_ * std_radrd_;
    
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
    
     
    // Initialize Time stamp
    time_us_ = meas_package.timestamp_;

    // State and Covariance Initialized
    is_initialized_ = true;
  }
  else
  {
    // compute the time elapsed between the current and previous measurements
    // dt - expressed in seconds
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    // Predict (based on motion model)
    this->Prediction(dt);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Update
      this->UpdateRadar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      // Update
      this->UpdateLidar(meas_package);
    }
  }

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */


  /******** Generate (augmented) Sigma Points ********/

  
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);


  // create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;  // normal distributed longitudinal acceleration noise with mean = 0
  x_aug(6) = 0;  // normal distributed yaw acceleration noise with mean = 0

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  // calculate square root of P
  MatrixXd A_aug = P_aug.llt().matrixL();

  // first column = x
  Xsig_aug.col(0) << x_aug;

  // set sigma points as columns of matrix Xsig

  for (int i = 0; i < n_aug_; ++i) {
    /* from 2 to n_x +1 */
    Xsig_aug.col(i + 1)     = x_aug + sqrt(lambda_ + n_aug_) * A_aug.col(i);

    /* from n_x + 2 to 2 * n_x + 1 */
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A_aug.col(i);
  }




  /******** Sigma Points Prediction ********/


  // predict sigma points
  double p_x = 0.0;
  double p_y = 0.0;
  double v = 0.0;
  double yaw = 0.0;
  double yawd = 0.0;
  double noise_a = 0.0;
  double noise_yawdd = 0.0;


  double p_x_pred = 0.0;
  double p_y_pred = 0.0;
  double v_pred = 0.0;
  double yaw_pred = 0.0;
  double yawd_pred = 0.0;

  Xsig_pred_.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    p_x = Xsig_aug(0,i);
    p_y = Xsig_aug(1,i);
    v   = Xsig_aug(2,i);
    yaw = Xsig_aug(3,i);
    yawd = Xsig_aug(4,i);
    noise_a = Xsig_aug(5,i);
    noise_yawdd = Xsig_aug(6,i);


    // avoid division by zero
    if(fabs(yawd) > 0.001)
    {
      p_x_pred = p_x + v/yawd * ( sin(yaw + yawd * delta_t) - sin(yaw) );
      p_y_pred = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
    } 
    else
    {
      p_x_pred = p_x + v * delta_t * cos(yaw);
      p_y_pred = p_y + v * delta_t * sin(yaw);
    }

    v_pred = v;
    yaw_pred = yaw + yawd * delta_t;
    yawd_pred = yawd;


    // add noise
    p_x_pred = p_x_pred + 0.5 * noise_a * delta_t * delta_t * cos(yaw);
    p_y_pred = p_y_pred + 0.5 * noise_a * delta_t * delta_t * sin(yaw);
    v_pred = v_pred + noise_a * delta_t;

    yaw_pred = yaw_pred + 0.5 * noise_yawdd * delta_t * delta_t;
    yawd_pred = yawd_pred + noise_yawdd * delta_t;


    // write predicted sigma points into matrix columns
    Xsig_pred_(0,i) = p_x_pred;
    Xsig_pred_(1,i) = p_y_pred;
    Xsig_pred_(2,i) = v_pred;
    Xsig_pred_(3,i) = yaw_pred;
    Xsig_pred_(4,i) = yawd_pred; 

  }



  /******** Predict Mean and Covariance *********/

  // predict state mean

  x_.fill(0.0); // Initialize all elements to 0.0

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) // loop over all sigma points
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);  
  }


  // predict state covariance matrix

  P_.fill(0.0); // Initialize all elements to 0.0

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) // loop over all sigma points
  {
    // Ok, but state contains angle 
    // P = P + weights(i) * (Xsig_pred.col(i) - x) * (Xsig_pred.col(i) - x).transpose();  
    // Normalization required when calculating the difference between angles

    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3) > M_PI)  x_diff(3)-=2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  // set measurement dimension, Lidar can measure (in this example) position x and y
  int n_z = 2;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ Transform prediction into measurement / sensor space +++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // ** In this case LiDAR can measure in the state space ** 

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) // loop over all sigma points
  {

    Zsig.col(i)(0) = Xsig_pred_.col(i)(0);
    Zsig.col(i)(1) = Xsig_pred_.col(i)(1);
  }

  // ***** calculate mean predicted measurement *****

  z_pred.fill(0.0); // Initialize all elements to 0.0

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) // loop over all sigma points
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);  
  }

  // ***** calculate covariance matrix S predicted measurement *****

  // measurement noise covariance matrix R
  MatrixXd R = MatrixXd(n_z,n_z);

  R.fill(0.0); // Init all elements with 0.0

  R(0,0) = std_laspx_ * std_laspx_; // variance of x position measured by laser
  R(1,1) = std_laspy_ * std_laspy_; // variance of y position measured by laser
  


  S.fill(0.0);  // Init all elements with 0.0

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) // loop over all sigma points
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization -> for LiDAR not needed
    //while (z_diff(1) > M_PI)  z_diff(1)-=2.*M_PI;
    //while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;

    S = S + ( weights_(i) * z_diff * z_diff.transpose() );
  }

  // add measurement noise covariance

  S = S + R;



  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ Update state with measurement                        +++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  // Current Lidar Measurement into vector

  VectorXd z = VectorXd(n_z);
  z <<
     meas_package.raw_measurements_[0], // px
     meas_package.raw_measurements_[1]; // py
     

  // calculate cross correlation matrix
  Tc.fill(0.0);


  for (int i = 0; i < 2 * n_aug_ + 1; ++i) // loop over all sigma points
  {
    // Calc z and x difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization -> for LiDAR not needed
    //while (z_diff(1) > M_PI)  z_diff(1)-=2.*M_PI;
    //while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;

    while (x_diff(3) > M_PI)  x_diff(3)-=2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3)+=2.*M_PI;

    // calc actuall cross correlation matrix
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  // calculate Kalman gain K;

  // create matrix for Kalman gain K
  MatrixXd K = MatrixXd(n_x_, n_z);

  K = Tc * S.inverse();

  // update state mean and covariance matrix

  // Calc z difference
  VectorXd z_diff = z - z_pred;

  // angle normalization -> for LiDAR not needed
  //while (z_diff(1) > M_PI)  z_diff(1)-=2.*M_PI;
  //while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;


  // update state mean
  x_ = x_ + K * z_diff;

  // update covariance matrix
  P_ = P_ - K * S * K.transpose();

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ Transform prediction into measurement / sensor space +++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // ***** transform sigma points into measurement space *****


  for (int i = 0; i < 2 * n_aug_ + 1; ++i) // loop over all sigma points
  {
    double p_x  = Xsig_pred_.col(i)(0);
    double p_y  = Xsig_pred_.col(i)(1);
    double v    = Xsig_pred_.col(i)(2);
    double yaw  = Xsig_pred_.col(i)(3);
    double yawd = Xsig_pred_.col(i)(4);

    double radar_radial_dist = sqrt(p_x * p_x + p_y * p_y);
    // double radar_phi = atan(p_y/p_x);
    double radar_yaw = atan2(p_y, p_x);
    double radar_radial_v = (p_x * cos(yaw) * v + p_y * sin(yaw) * v) / sqrt(p_x * p_x + p_y * p_y);

    Zsig.col(i)(0) = radar_radial_dist;
    Zsig.col(i)(1) = radar_yaw;
    Zsig.col(i)(2) = radar_radial_v;
  }


  // ***** calculate mean predicted measurement *****


  z_pred.fill(0.0); // Initialize all elements to 0.0

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) // loop over all sigma points
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);  
  }


  // ***** calculate covariance matrix S predicted measurement *****

  // measurement noise covariance matrix R
  MatrixXd R = MatrixXd(n_z,n_z);

  R.fill(0.0); // Init all elements with 0.0

  R(0,0) = std_radr_ * std_radr_; // variance of radial dist
  R(1,1) = std_radphi_ * std_radphi_; // variance of angle phi
  R(2,2) = std_radrd_ * std_radrd_; // variance of radial velocity

  

  S.fill(0.0);  // Init all elements with 0.0

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) // loop over all sigma points
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI)  z_diff(1)-=2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;

    S = S + ( weights_(i) * z_diff * z_diff.transpose() );
  }

  // add measurement noise covariance

  S = S + R;


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++ Update state with measurement                        +++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  // Current Radar Measurement into vector

  VectorXd z = VectorXd(n_z);
  z <<
     meas_package.raw_measurements_[0], // rho in m
     meas_package.raw_measurements_[1], // phi in rad
     meas_package.raw_measurements_[2]; // rho_dot in m/s


  // calculate cross correlation matrix
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) // loop over all sigma points
  {
    // Calc z and x difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI)  z_diff(1)-=2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;

    while (x_diff(3) > M_PI)  x_diff(3)-=2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3)+=2.*M_PI;

    // calc actuall cross correlation matrix
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  // calculate Kalman gain K;

  // create matrix for Kalman gain K
  MatrixXd K = MatrixXd(n_x_, n_z);

  K = Tc * S.inverse();

  // update state mean and covariance matrix

  // Calc z difference
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1) > M_PI)  z_diff(1)-=2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;


  // update state mean
  x_ = x_ + K * z_diff;

  // update covariance matrix
  P_ = P_ - K * S * K.transpose();
}