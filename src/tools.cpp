#include "./tools.h"
#include <iostream>
#include <vector>

using std::vector;

const float TWO_PI = 2 * M_PI;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::calculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd RMSE(4);
  RMSE << 0, 0, 0, 0;

  const std::size_t N = estimations.size();
  if (N == 0 || N != ground_truth.size()) {
    return RMSE;
  }

  // accumulate squared residuals
  for (std::size_t t = 0; t < N; ++t) {
    VectorXd res = estimations[t] - ground_truth[t];
    res = res.array().square();
    RMSE += res;
  }

  RMSE = (RMSE / N).array().sqrt();
  return RMSE;
}

float Tools::normalizeAngle(float angle) {
  if (angle < -M_PI || angle > M_PI) {
      angle -= TWO_PI * floor((angle + M_PI) / TWO_PI);
  }
  return angle;
}
