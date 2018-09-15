#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

//#define M_PI 3.1415926
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//#define DEBUG
/**
* Initializes Unscented Kalman filter
* This is scaffolding, do not modify
*/
UKF::UKF() {

  // Open NIS data files
  NISvals_radar_.open( "./NISvals_radar.txt", ios::out );
  NISvals_laser_.open( "./NISvals_laser.txt", ios::out );

  // Check for errors opening the files
  if( !NISvals_radar_.is_open() )
  {
    cout << "Error opening NISvals_radar.txt" << endl;
    exit(1);
  }

  if( !NISvals_laser_.is_open() )
  {
    cout << "Error opening NISvals_laser.txt" << endl;
    exit(1);
}
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;

  n_z_r = 3;
  n_z_l = 2;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    if (i == 0) weights_(i) = lambda_ / (lambda_ + n_aug_);
    else       weights_(i) = 0.5 / (lambda_ + n_aug_);
  }
  /*RADAR */
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;
#ifdef DEBUG
  cout << "initialization finished" << endl;
  cout << "n_x_: " << n_x_ << endl;
  cout << "n_aug_: " << n_aug_ << endl;
  cout << "lambda_: " << lambda_ << endl;
  cout << "Xsig_pred_ shape: " << Xsig_pred_.rows() << " , " << Xsig_pred_.cols() << endl;
  cout << "weights_: " << weights_ << endl;
  cout << "R_radar_: " << R_radar_ << endl;
  cout << "R_lidar_: " << R_lidar_ << endl;
#endif // DEBUG

}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} meas_package The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
  if (!is_initialized_)
  {
    x_.fill(0.0);
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    else
    {
      x_(0) = meas_package.raw_measurements_(0)*cos(meas_package.raw_measurements_(1));
      x_(1) = meas_package.raw_measurements_(0)*sin(meas_package.raw_measurements_(1));
      x_(2) = 0;
      x_(3) = meas_package.raw_measurements_(1);
      x_(4) = meas_package.raw_measurements_(2);
    }
    P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  /*****************************************************************************
  *  Prediction
  ****************************************************************************/
  //compute the time elapsed between the current and previous measurements
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;  //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (use_laser_  && meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
  else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else
  {
    cout << "--------------------------------ERROR------------------------------------" << endl;
  }
}

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} delta_t the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  /*****************************************************************************
  *  Generate sigma Point
  ****************************************************************************/
  /*Predict parameter*/
#ifdef DEBUG
  cout << "--------------------------------Predict Stage------------------------------------" << endl;
  cout << "n_x_" << n_x_ << endl;
  cout << "n_aug_" << n_aug_ << endl;
  cout << "lambda_" << lambda_ << endl;
#endif
  VectorXd x_aug = VectorXd(n_aug_);                //7
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);        // 7*7
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);//7*15

  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.fill(0.0);
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }
  /*****************************************************************************
  *  Predict sigma Point
  ****************************************************************************/
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double pred_x, pred_y;
    if (fabs(yawd) > 0.001)
    {
      pred_x = p_x + ((v / yawd)*(sin(yaw + yawd*delta_t) - sin(yaw))) + (0.5*delta_t*delta_t*cos(yaw)*nu_a);
      pred_y = p_y + ((v / yawd)*(-cos(yaw + yawd*delta_t) + cos(yaw))) + (0.5*delta_t*delta_t*sin(yaw)*nu_a);
    }
    else
    {
      pred_x = p_x + (v*cos(yaw)*delta_t) + (0.5*delta_t*delta_t*cos(yaw)*nu_a);
      pred_y = p_y + (v*sin(yaw)*delta_t) + (0.5*delta_t*delta_t*sin(yaw)*nu_a);
    }
    double pred_v = v + delta_t * nu_a;
    double pred_yaw = yaw + yawd*delta_t + 0.5*delta_t*delta_t*nu_yawdd;
    double pred_yawd = yawd + delta_t * nu_yawdd;

    Xsig_pred_(0, i) = pred_x;
    Xsig_pred_(1, i) = pred_y;
    Xsig_pred_(2, i) = pred_v;
    Xsig_pred_(3, i) = pred_yaw;
    Xsig_pred_(4, i) = pred_yawd;
  }
  /*****************************************************************************
  * Convert Predicted Sigma Points to Mean/Covariance
  ****************************************************************************/
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  }
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)  > M_PI)  x_diff(3) -= 2.*M_PI;
    while (x_diff(3)  < -M_PI)  x_diff(3) += 2.*M_PI;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
#ifdef DEBUG
  cout << "x_ shape : " << x_.rows() << " , " << x_.cols() << endl;
  cout << "P_ shape : " << P_.rows() << " , " << P_.cols() << endl;
  cout << "Xsig_pred_ shape : " << Xsig_pred_.rows() << " , " << Xsig_pred_.cols() << endl;

  cout << "x_" << x_ << endl;
  cout << "P_" << P_ << endl;
  cout << "+++++++++++++++++++++++++++finish prediction one time+++++++++++++++++++++++++++" << endl;
#endif
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  /*****************************************************************************
  * Predicted Radar Measurement
  ****************************************************************************/
#ifdef DEBUG
  cout << "--------------------------------UpdateLidar Stage------------------------------------" << endl;
  cout << "n_z_l" << n_x_ << endl;
  cout << "n_aug_" << n_aug_ << endl;
#endif
  MatrixXd Zsig_pred = MatrixXd(n_z_l, 2 * n_aug_ + 1);
  MatrixXd z_prd = VectorXd(n_z_l);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    Zsig_pred(0, i) = p_x;
    Zsig_pred(1, i) = p_y;
  }
#ifdef DEBUG
  cout << "Zsig_pred shape : " << Zsig_pred.rows() << " , " << Zsig_pred.cols() << endl;
#endif
  //calculate mean predicted measurement
  z_prd.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    z_prd = z_prd + weights_(i) * Zsig_pred.col(i);
  }
#ifdef DEBUG
  cout << "z_prd shape : " << z_prd.rows() << " , " << z_prd.cols() << endl;
#endif
  //calculate innovation covariance matrix S
  MatrixXd S_lasr_ = MatrixXd(n_z_l, n_z_l);
  S_lasr_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd z_diff = Zsig_pred.col(i) - z_prd;
    S_lasr_ = S_lasr_ + weights_(i) * z_diff * z_diff.transpose();
  }
  S_lasr_ = S_lasr_ + R_lidar_;
#ifdef DEBUG
  cout << "S_lasr_ shape : " << S_lasr_.rows() << " , " << S_lasr_.cols() << endl;
#endif
  /*****************************************************************************
  * UKF Update
  ****************************************************************************/
  //VectorXd z = meas_package.raw_measurements_;
  VectorXd z = VectorXd(n_z_l);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
#ifdef DEBUG
  cout << "z shape : " << z.rows() << " , " << z.cols() << endl;
#endif
  MatrixXd Tc = MatrixXd(n_x_, n_z_l);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd z_diff = Zsig_pred.col(i) - z_prd;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    /*
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    */
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K
  MatrixXd K = Tc * S_lasr_.inverse();
  //update state mean and covariance matrix
  VectorXd z_diff = z - z_prd;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_lasr_ * K.transpose();
  NIS_laser_ = z_diff.transpose()*S_lasr_.inverse()*z_diff;
  // cout << NIS_laser_ << endl;
  NISvals_laser_ << NIS_laser_ << endl;;
#ifdef DEBUG
  cout << "x_ shape : " << x_.rows() << " , " << x_.cols() << endl;
  cout << "P_ shape : " << P_.rows() << " , " << P_.cols() << endl;
  cout << "+++++++++++++++++++++++++++Finish UpdateLidar Stage+++++++++++++++++++++++++++" << endl;
#endif
}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  /*****************************************************************************
  * Predicted Radar Measurement
  ****************************************************************************/
#ifdef DEBUG
  cout << "--------------------------------UpdateRadar Stage------------------------------------" << endl;
  cout << "n_z_r" << n_x_ << endl;
  cout << "n_aug_" << n_aug_ << endl;
#endif
  MatrixXd Zsig_pred = MatrixXd(n_z_r, 2 * n_aug_ + 1);
  VectorXd z_prd = VectorXd(n_z_r);
  Zsig_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double rho = sqrt(p_x*p_x + p_y*p_y);
    double phi = atan2(p_y, p_x);
    double rho_dot = (p_x*v*cos(yaw) + p_y*v*sin(yaw)) / rho;
    if(fabs(rho) < 0.0001)
    {
      rho_dot = 0;
    }
    else
    {
    rho_dot = (p_x*v*cos(yaw) + p_y*v*sin(yaw))/rho;
    }
    Zsig_pred(0, i) = rho;
    Zsig_pred(1, i) = phi;
    Zsig_pred(2, i) = rho_dot;
  }
#ifdef DEBUG
  cout << "Zsig_pred shape : " << Zsig_pred.rows() << " , " << Zsig_pred.cols() << endl;
#endif
  //calculate mean predicted measurement
  z_prd.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    z_prd = z_prd + weights_(i)*Zsig_pred.col(i);
  }
#ifdef DEBUG
  cout << "z_prd shape : " << z_prd.rows() << " , " << z_prd.cols() << endl;
#endif

  while (z_prd(1) > M_PI) z_prd(1) -= 2.*M_PI;
  while (z_prd(1) < -M_PI) z_prd(1) += 2.*M_PI;
  //calculate innovation covariance matrix S
  MatrixXd S_radr = MatrixXd(n_z_r, n_z_r);
  S_radr.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd z_diff = Zsig_pred.col(i) - z_prd;
    while (z_diff(1) > M_PI)  z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S_radr = S_radr + weights_(i)*z_diff*z_diff.transpose();
  }
  S_radr = S_radr + R_radar_;
#ifdef DEBUG
  cout << "S_radr shape : " << S_radr.rows() << " , " << S_radr.cols() << endl;
#endif
  /*****************************************************************************
  * UKF Update
  ****************************************************************************/
  VectorXd z = VectorXd(n_z_r);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);
#ifdef DEBUG
  cout << "z shape : " << z.rows() << " , " << z.cols() << endl;
#endif
  MatrixXd Tc = MatrixXd(n_x_, n_z_r);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd z_diff = Zsig_pred.col(i) - z_prd;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
#ifdef DEBUG
  cout << "Tc shape(5, 3) : " << Tc.rows() << " , " << Tc.cols() << endl;
#endif
  //calculate Kalman gain K
  MatrixXd K = Tc * S_radr.inverse();
  //update state mean and covariance matrix
  VectorXd z_diff = z - z_prd;
  while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_radr * K.transpose();
  NIS_radar_ = z_diff.transpose()*S_radr.inverse()*z_diff;
  //cout << NIS_radar_ << endl;
  NISvals_radar_ << NIS_radar_ << endl;
#ifdef DEBUG
  cout << "x_ shape : " << x_.rows() << " , " << x_.cols() << endl;
  cout << "P_ shape : " << P_.rows() << " , " << P_.cols() << endl;
  cout << "+++++++++++++++++++++++++++Finish UpdateRadar Stage+++++++++++++++++++++++++++" << endl;
#endif
}
