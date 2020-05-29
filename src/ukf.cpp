#include "ukf.h"
#include "Eigen/Dense"

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

  is_initialized_ = false;

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  // check if the first initialization is done
  // if not initialized, initialize with the first measurement
  if(false == is_initialized_)
  {
    if(MeasurementPackage::SensorType::LASER == meas_package.sensor_type_)
    {
      //initialize state
      x_.block(0,0,2,0) = meas_package.raw_measurements_;
      use_laser_ = false;
    }
    else if (MeasurementPackage::SensorType::RADAR == meas_package.sensor_type_)
    {
      x_(0) = meas_package.raw_measurements_(0)*cos(meas_package.raw_measurements_(1));
      x_(1) = meas_package.raw_measurements_(0)*sin(meas_package.raw_measurements_(1));
      x_(2) = meas_package.raw_measurements_(2);
      
      use_radar_ = false;

    }
    else
    {
      /* do nothing, unsupported type */
    }
    
    //initialize covariance 

    P_(0,0) = 1;
    P_(1,1) = 1;
    P_(2,2) = 1;
    P_(3,3) = 1;
    P_(4,4) = 1;

    //set time of last measurement
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
  }
  //if this is not the first measurement
  else
  {
    //check if the current measurement is the one to use
    if((MeasurementPackage::SensorType::LASER == meas_package.sensor_type_ && true == use_laser_)
    || (MeasurementPackage::SensorType::RADAR == meas_package.sensor_type_ && true == use_radar_))
    {
      // set time difference for prediction
      double dt = meas_package.timestamp_ - time_us_;

      //run prediction step
      Prediction(dt);

      //check which measurement to process'
      if(true == use_laser_)
      {
        UpdateLidar(meas_package);
        use_laser_ = false;
        use_radar_ = true;
      }
      else if(true == use_radar_)
      {
        UpdateRadar(meas_package);
        use_radar_ = false;
        use_laser_ = true;
      }
      else
      {
        /* this should not happen */
      }
      
    }
    
    
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  //create augmented state and covariance
  VectorXd x_aug(n_aug_);
  MatrixXd P_aug(n_aug_,n_aug_);

  x_aug.block(0,0,n_x_,0) = x_;
  P_aug.block(0,0,n_x_,n_x_) = P_;

  //generate sigma points
  int n_sigma = (2*n_aug_) + 1;
  MatrixXd X_sigma(n_x_,n_sigma);
  MatrixXd X_sigma_aug(n_aug_,n_sigma);
  MatrixXd P_sqrt = P_aug.llt().matrixL();

  //insert mean as first sigma point
  X_sigma.col(0) = x_aug;

  //generate the rest of sigma points
  for(int i = 1; i < n_aug_; i++)
  {
    X_sigma_aug.col(i) = x_aug + sqrt(lambda_+n_aug_)*P_sqrt.col(i-1);
    X_sigma_aug.col(n_aug_ + i) = x_aug - sqrt(lambda_+n_aug_)*P_sqrt.col(i-1);
  }

  //calculate weights
  double w_0 = lambda_/(lambda_+n_aug_);
  double w_i = 1/(2*(lambda_+n_aug_));

  //calculate prediction for each sigma point
  for(int i = 0; i < n_sigma; i++)
  {
      double px = X_sigma_aug(0,i);
      double py = X_sigma_aug(1,i);
      double v =  X_sigma_aug(2,i);
      double yaw = X_sigma_aug(3,i);
      double yaw_d = X_sigma_aug(4,i);
      double new_a = X_sigma_aug(5,i);
      double new_phi = X_sigma_aug(6,i);
      
      //avoid division by zero
      if(yaw_d == 0)
      {
          X_sigma(0, i) = px + v * cos(yaw)*delta_t + 0.5*delta_t*delta_t*cos(yaw)*new_a;
          X_sigma(1, i) = py + v * sin(yaw)*delta_t + 0.5*delta_t*delta_t*cos(yaw)*new_a;   
      }
      else
      {
          X_sigma(0, i) = px + (v/yaw_d)*(sin(yaw+yaw_d*delta_t)-sin(yaw)) + 0.5*delta_t*delta_t*cos(yaw)*new_a;
          X_sigma(1, i) = py + (v/yaw_d)*(-cos(yaw+yaw_d*delta_t)+cos(yaw)) + 0.5*delta_t*delta_t*sin(yaw)*new_a;
          
      }
      
        X_sigma(2, i) = v + 0 + delta_t*new_a;
        X_sigma(3, i) = yaw + yaw_d * delta_t + 0.5*delta_t*delta_t*new_phi;
        X_sigma(4, i) = yaw_d + 0 + delta_t*new_phi;
  }

  //calculate the predicted mean and covariance
    
  for(int i = 0; i < n_sigma; i++ )
  {
    double w = (i == 0 ? w_0 : w_i);
    x_ = x_ + (w*X_sigma.col(i));   
  }
  
  for(int i = 0; i < n_sigma; i++ )
  {
    double w = (i == 0 ? w_0 : w_i);
    VectorXd x_diff = X_sigma.col(i) - x_; 
    P_ = P_ +  w*x_diff*x_diff.transpose();     
  }


}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}