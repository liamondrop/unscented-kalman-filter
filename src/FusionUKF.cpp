#include "./FusionUKF.h"

#include <iostream>

/**
 * Constructor.
 */
FusionUKF::FusionUKF() {
  use_radar_ = true;
  use_laser_ = true;
  n_x_ = 5;

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

  // laser measurement covariance matrix
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  // radar measurement covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_,   0,   0,
              0, std_radphi_ * std_radphi_, 0,
              0,  0,  std_radrd_ * std_radrd_;
}

/**
 * Destructor.
 */
FusionUKF::~FusionUKF() {}

VectorXd FusionUKF::getInitialX_(const MeasurementPackage &meas_package) {
  VectorXd initial_x(n_x_);
  initial_x.fill(0.0);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    const float rho = meas_package.raw_measurements_[0];
    const float phi = meas_package.raw_measurements_[1];
    initial_x(0) = rho * sin(phi);
    initial_x(1) = rho * cos(phi);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    initial_x(0) = meas_package.raw_measurements_[0];
    initial_x(1) = meas_package.raw_measurements_[1];
  }
  return initial_x;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void FusionUKF::processMeasurement(const MeasurementPackage &meas_package) {
  if ((!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) ||
      (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)) {
    return;
  }

  // compute the time elapsed between the current and previous measurements
  const float delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1e+6;
  previous_timestamp_ = meas_package.timestamp_;

  // check if is initialized and the elapsed time is a sane value
  if (!is_initialized_ || delta_t <= 0 || delta_t > 100) {
    ukf_.x_ = getInitialX_(meas_package);

    // done initializing, no need to predict or update
    is_initialized_ = true;
  } else {
    // prediction step
    ukf_.predict(delta_t);

    // update step
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      ukf_.updateLidar(meas_package.raw_measurements_, R_laser_);
    } else {
      ukf_.updateRadar(meas_package.raw_measurements_, R_radar_);
    }
  }

  std::cout << "x =\n" << ukf_.x_ << std::endl;
  std::cout << "P =\n" << ukf_.P_ << std::endl;
  std::cout << "=========================" << std::endl;
}
