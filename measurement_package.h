//
// Created by denise on 3/22/17.
//

#ifndef UNSCENTED_KALMAN_FILTER_MEASUREMENT_PACKAGE_H
#define UNSCENTED_KALMAN_FILTER_MEASUREMENT_PACKAGE_H


#include "Eigen/Dense"

class MeasurementPackage {
public:
    long timestamp_;

    enum SensorType{
        LASER,
        RADAR
    } sensor_type_;

    Eigen::VectorXd raw_measurements_;

};


#endif //UNSCENTED_KALMAN_FILTER_MEASUREMENT_PACKAGE_H
