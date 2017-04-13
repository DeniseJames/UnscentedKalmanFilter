
#ifndef UNSCENTED_KALMAN_FILTER_GROUND_TRUTH_PACKAGE_H
#define UNSCENTED_KALMAN_FILTER_GROUND_TRUTH_PACKAGE_H

//#pragma once
#include "Eigen/Dense"

class GroundTruthPackage {
public:
	long timestamp_;

	enum SensorType {
		LASER,
		RADAR
	} sensor_type_;

	Eigen::VectorXd gt_values_;

};

#endif //UNSCENTED_KALMAN_FILTER_GROUND_TRUTH_PACKAGE_H
