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
      }
      else if(true == use_radar_)
      {
        UpdateRadar(meas_package);
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