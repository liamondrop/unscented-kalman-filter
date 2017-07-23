import numpy as np
from scipy import linalg

TWO_PI = 2 * np.pi

class UKF(object):
    def __init__(self, n, m, r_var=0.1, p_var=1e+3, q_var=0.1):
        """
        Creates an Unscented Kalman Filter with n states and m observations
        """
        # size of the state vector
        self.n_x = n

        # size of augmented covariance matrix
        self.n_aug = n + m

        # number of sigma points (2n+1)
        self.n_sig = 2 * self.n_aug + 1

        # predicted sigma points matrix
        self.X_sigma = np.zeros((n, self.n_sig))

        # sigma point spreading parameter
        self._lambda_ = 3. - self.n_aug

        # weight vector
        self.weights = np.zeros((1, self.n_sig))
        self.weights[:, :1] = self._lambda_ / (self._lambda_ + self.n_aug)
        self.weights[:, 1:] = .5 / (self._lambda_ + self.n_aug)

        # measurement noise covariance
        self.R = np.eye(m) * r_var

        # process noise covariance
        self.Q = np.eye(m) * q_var

        self.x = np.ones((n, 1))   # initial state
        self.P = np.eye(n) * p_var # initial prediction uncertainty covariance

        # Normalized Innovation Squared (calculated in update step)
        # this value will help us evaluate the consistency of our predictions
        # and determine whether our noise parameters are good
        self.NIS = 0.

    def predict(self, delta_t):
        """
        Predicts sigma points matrix for time step t+1
        """
        # create augmented mean vector
        x_aug = np.zeros((self.n_aug, 1))
        x_aug[:self.n_x] = self.x

        # create augmented state covariance
        P_aug = np.zeros((self.n_aug, self.n_aug))
        P_aug[:self.n_x, :self.n_x] = self.P
        P_aug[self.n_x:, self.n_x:] = self.Q

        # take square root (cholesky) of P_aug and modify by spreading param
        P_aug_root = linalg.sqrtm((self._lambda_ + self.n_aug) * P_aug)

        # create augmented sigma point matrix
        X_sigma_aug = np.zeros((self.n_aug, self.n_sig))
        X_sigma_aug[:, :1] = x_aug
        X_sigma_aug[:, 1:self.n_aug+1] = x_aug + P_aug_root
        X_sigma_aug[:,  self.n_aug+1:] = x_aug - P_aug_root

        # predict sigma points
        for col in range(self.n_sig):
            px, py, v, yaw, yawd, nu_a, nu_yawdd = X_sigma_aug[:, col]

            # predicted state values for CTRV motion model
            # if yawd is a very small value, it implies a nearly straight line,
            # so degenerate to simpler CV motion model in that case
            if abs(yawd) > 0.001:
                px_p = px + v / yawd * (np.sin(yaw + yawd * delta_t) - np.sin(yaw))
                py_p = py + v / yawd * (np.cos(yaw) - np.cos(yaw + yawd * delta_t))
            else:
                px_p = px + v * delta_t * np.cos(yaw)
                py_p = py + v * delta_t * np.sin(yaw)

            v_p = v
            yaw_p = yaw + yawd * delta_t
            yawd_p = yawd

            # add acceleration noise
            px_p   += .5 * nu_a * delta_t**2 * np.cos(yaw)
            py_p   += .5 * nu_a * delta_t**2 * np.sin(yaw)
            v_p    += nu_a * delta_t
            yaw_p  += .5 * nu_yawdd * delta_t**2
            yawd_p += nu_yawdd * delta_t

            self.X_sigma[:, col] = px_p, py_p, v_p, yaw_p, yawd_p

        # predicted x state vector
        self.x = np.sum(self.X_sigma * self.weights, axis=1).reshape(-1, 1)

        # state residual
        x_diff = self.X_sigma - self.x

        # normalize the predicted yaw component to between -pi & pi
        x_diff[3] -= TWO_PI * np.floor((x_diff[3] + np.pi) / TWO_PI)

        # predicted P uncertainty covariance
        self.P = (self.weights * x_diff).dot(x_diff.T)

    def update(self, z):
        """
        Update predictions with observed measurement z
        """
        nz = len(z)         # measurement dimension
        z = np.array([z]).T # convert measurement to single column matrix

        # create matrix for sigma points in measurement space
        # from the x and y sigma point predictions
        Z_sigma = self.X_sigma[:nz]

        # mean predicted measurement
        z_pred = np.sum(Z_sigma * self.weights[:nz], axis=1).reshape(-1, 1)

        # measurement sigma and state residuals
        z_sig_diff = Z_sigma - z_pred
        x_diff = self.X_sigma - self.x
        # x_diff[3] -= TWO_PI * np.floor((x_diff[3] + np.pi) / TWO_PI)

        # compute cross correlation matrix
        Tc = (self.weights * x_diff).dot(z_sig_diff.T)

        # compute measurement covariance matrix + measurement noise
        S  = (self.weights * z_sig_diff).dot(z_sig_diff.T) + self.R
        S_inv = np.linalg.inv(S)

        # compute the Kalman gain
        K = Tc.dot(S_inv)

        # measurement residual
        z_diff = z - z_pred

        # update x and prediction uncertainty
        self.x += K.dot(z_diff)
        self.P -= K.dot(S).dot(K.T)

        # update NIS calculation
        self.NIS = z_diff.T.dot(S_inv).dot(z_diff)
