#include <iostream>

#include "./ukf.h"

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  n_x_ = 5;
  n_aug_ = 7;
  n_sig_ = 2 * n_aug_ + 1;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

  // create vector for weights
  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_sig_; i++) {  // 2n+1 weights
    weights_(i) = 0.5 / (n_aug_ + lambda_);
  }

  // NIS for radar
  NIS_radar_ = 0.0;

  // NIS for laser
  NIS_laser_ = 0.0;
}


UKF::~UKF() {}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {float} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::predict(const float &delta_t) {
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  // create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_,   n_x_)   = std_a_ * std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

  // create square root matrix
  MatrixXd P_root = P_aug.llt().matrixL();
  P_root *= sqrt(lambda_ + n_aug_);

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + P_root.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - P_root.col(i);
  }

  // predict sigma points
  for (int i = 0; i < n_sig_; i++) {
    double p_x      = Xsig_aug(0, i);
    double p_y      = Xsig_aug(1, i);
    double v        = Xsig_aug(2, i);
    double yaw      = Xsig_aug(3, i);
    double yawd     = Xsig_aug(4, i);
    double nu_a     = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p += nu_a * delta_t;

    yaw_p += 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p += nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_.col(i) << px_p,
                         py_p,
                         v_p,
                         yaw_p,
                         yawd_p;
  }

  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = tools.normalizeAngle(x_diff(3));

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}


void UKF::updateLidar(const VectorXd &z, const MatrixXd &R) {
  // set measurement dimension, lidar can measure p_x and p_y
  size_t n_z = z.size();

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    // measurement model
    Zsig.col(i) << p_x,
                   p_y;
    z_pred += weights_(i) * Zsig.col(i);
  }

  // calculate cross correlation matrix Tc and measurement covariance matrix S
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  Tc.fill(0.0);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc += weights_(i) * x_diff * z_diff.transpose();
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise
  S += R;

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // calculate NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}


void UKF::updateRadar(const VectorXd &z, const MatrixXd &R) {
  // set measurement dimension, radar can measure r, phi, and r_dot
  size_t n_z = z.size();

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    double rho = sqrt(p_x*p_x + p_y*p_y);
    double phi = atan2(p_y, p_x);
    double rho_dot = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);

    // measurement model
    Zsig.col(i) << rho, phi, rho_dot;

    z_pred += weights_(i) * Zsig.col(i);
  }

  // calculate cross correlation matrix Tc and measurement covariance matrix S
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  Tc.fill(0.0);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = tools.normalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = tools.normalizeAngle(x_diff(3));

    S += weights_(i) * z_diff * z_diff.transpose();
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // add measurement noise
  S += R;

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;
  z_diff(1) = tools.normalizeAngle(z_diff(1));

  // calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}
