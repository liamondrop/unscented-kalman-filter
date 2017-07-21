#ifndef SRC_UKF_H_
#define SRC_UKF_H_

#include <fstream>
#include <string>
#include <vector>

#include "./measurement_package.h"
#include "./tools.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
 public:
  Tools tools;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  // state covariance matrix
  MatrixXd P_;

  // predicted sigma points matrix
  MatrixXd Xsig_pred_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Weights of sigma points
  VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Number of sigma points
  int n_sig_;

  // Sigma point spreading parameter
  double lambda_;

  // the NIS for radar
  double NIS_radar_;

  // the NIS for laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * Prediction Predicts sigma points, the state,
   * and the state covariance matrix
   * @param delta_t Time between k and k+1 in s
   */
  void predict(const float &delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param z The measurement at k+1
   * @param R The measurement covariance matrix
   */
  void updateLidar(const VectorXd &z, const MatrixXd &R);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param z The measurement at k+1
   * @param R The measurement covariance matrix
   */
  void updateRadar(const VectorXd &z, const MatrixXd &R);
};

#endif  // SRC_UKF_H_
