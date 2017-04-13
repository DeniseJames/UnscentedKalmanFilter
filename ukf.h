//
// Created by denise on 3/22/17.
//

#ifndef UNSCENTED_KALMAN_FILTER_UKF_H
#define UNSCENTED_KALMAN_FILTER_UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

    ///* initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_;

    // previous timestamp
    long long previous_timestamp_;

    //time difference
   

    ///* if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;

    ///* if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;

    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;

	VectorXd z_radar_;

    ///* state covariance matrix
    MatrixXd P_;

	// * laser state variable with vx and vy
	VectorXd x_vxvy_;

    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;

    /// sigma vector aug
    //create augmented mean vector
    VectorXd x_aug_;

    //create augmented state covariance
    MatrixXd P_aug_;

    //create sigma point matrix
    MatrixXd Xsig_aug_;

    //Measurement model radar vector
    MatrixXd Zsig_;

    //
    int n_z_;

    ///* time when the state is true, in us
    long time_us_;

    ///* Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a_;

    ///* Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd_;

    ///* Laser measurement noise standard deviation position1 in m
    double std_laspx_;

    ///* Laser measurement noise standard deviation position2 in m
    double std_laspy_;

    ///* Radar measurement noise standard deviation radius in m
    double std_radr_;

    ///* Radar measurement noise standard deviation angle in rad
    double std_radphi_;

    ///* Radar measurement noise standard deviation radius change in m/s
    double std_radrd_ ;

    ///* Weights of sigma points
    VectorXd weights_;

    // radar prediction vector
    VectorXd z_pred_;

    //create example matrix for predicted measurement covariance
    MatrixXd S_ ;

    //create matrix for cross correlation Tc
    MatrixXd Tc_;

    //Lidar H matrix
	// z is the measured value of the sensors
	 MatrixXd H_ ;
	
	 // lidar R function
	 MatrixXd R_laser_;
//    Need this to make delta_t_ accessible to other functions
    double delta_t_;

    ///* State dimension
    int n_x_;

    ///* Augmented state dimension
    int n_aug_;

    ///* Sigma point spreading parameter
    double lambda_;

    ///* the current NIS for radar
    double NIS_radar_;

    ///* the current NIS for laser
    double NIS_laser_;

    // incoming radar measurement
    VectorXd z_;

	VectorXd estm_;

    /**
     * Constructor
     */
    UKF();

    /**
     * Destructor
     */
    virtual ~UKF();

    /**
     * ProcessMeasurement
     * @param meas_package The latest measurement data of either radar or laser
     */
    void ProcessMeasurement(const MeasurementPackage meas_package);

    /**
     * Prediction Predicts sigma points, the state, and the state covariance
     * matrix
     * @param delta_t Time between k and k+1 in s
     */



    void AugmentedSigmaPoints();
    void SigmaPointPrediction(const double delta_t);
    void PredictMeanAndCovariance();

    /**
     * Updates the state and the state covariance matrix using a laser measurement
     * @param meas_package The measurement at k+1
     */

    void PredictRadarMeasurement();

	void UpdateLidar(MeasurementPackage meas_package);

    /**
     * Updates the state and the state covariance matrix using a radar measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateRadar(MeasurementPackage meas_package);
 
//	VectorXd UKF::CalculateRMSE(const std::vector<VectorXd> &estimations,
//		const std::vector<VectorXd> &ground_truth);
//
//VectorXd UKF::CalculateRMSE(const std::vector<VectorXd> &estimations,
//	const std::vector<VectorXd> &ground_truth);

    VectorXd CalculateRMSE(const std::vector<VectorXd> &estimations,
                                const std::vector<VectorXd> &ground_truth);
};


#endif //UNSCENTED_KALMAN_FILTER_UKF_H
