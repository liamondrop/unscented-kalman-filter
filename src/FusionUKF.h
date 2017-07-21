#ifndef SRC_FUSIONUKF_H_
#define SRC_FUSIONUKF_H_

#include "./ukf.h"
#include "./measurement_package.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

class FusionUKF {
 public:
  /**
   * Constructor.
   */
  FusionUKF();

  /**
   * Destructor.
   */
  virtual ~FusionUKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void processMeasurement(const MeasurementPackage &meas_package);

  /**
   * Unscented Kalman Filter
   */
  UKF ukf_;

 private:
  bool is_initialized_;
  int64_t previous_timestamp_;

  int n_x_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // Laser measurement noise standard deviation position1 in m
  float std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  float std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  float std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  float std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  float std_radrd_ = 0.3;

  // Measurement noise for Lidar
  MatrixXd R_laser_;

  // Measurement noise for Radar
  MatrixXd R_radar_;

  /**
   * Get the initial state vector from the measurement
   * @param meas_package The measurement at k+1
   */
  VectorXd getInitialX_(const MeasurementPackage &meas_package);
};

#endif  // SRC_FUSIONUKF_H_
