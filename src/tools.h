#ifndef SRC_TOOLS_H_
#define SRC_TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

class Tools {
 public:
  /**
   * Constructor.
   */
  Tools();

  /**
   * Destructor.
   */
  virtual ~Tools();

  /**
   * A helper method to calculate RMSE.
   */
  VectorXd calculateRMSE(const vector<VectorXd> &estimations,
                         const vector<VectorXd> &ground_truth);

  /**
   * normalize the given angle (in radians) to between -pi and pi
   */
  float normalizeAngle(float angle);

};

#endif  // SRC_TOOLS_H_
