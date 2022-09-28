
#pragma once

#include <Eigen/Dense>
#include <sensor_msgs/Imu.h>

namespace COMMON
{
  struct cmpVec3i
  {
    bool operator()(const Eigen::Vector3i &a, const Eigen::Vector3i &b) const
    {
      for (size_t i = 0; i < 3; ++i)
      {
        if (a[i] < b[i])
          return true;
        if (a[i] > b[i])
          return false;
      }
      return false;
    }
  };


  sensor_msgs::Imu imuCorrected(const sensor_msgs::Imu &imu, Eigen::Affine3d &T_ti)
  {
    //T_ti target frame -> imu frame
    //Get imu data in target frame.
    const Eigen::Vector3d acc(imu.linear_acceleration.x, imu.linear_acceleration.y, imu.linear_acceleration.z);
    const Eigen::Vector3d gyo(imu.angular_velocity.x, imu.angular_velocity.y, imu.angular_velocity.z);

    const Eigen::Matrix3d Rbi = T_ti.rotation().matrix();
    Eigen::Vector3d corrected_acc = Rbi * acc;
    Eigen::Vector3d corrected_gyo = Rbi * gyo;

    //const Eigen::Vector3d b_arm = T_ti.translation();
    //if (!b_arm.isZero())
    //{
    //  const Eigen::Matrix3d body_Omega_body = skew<double>(corrected_gyo);
    //  const Eigen::Vector3d b_velocity_bs = body_Omega_body * b_arm;
    //  corrected_acc -= body_Omega_body * b_velocity_bs;
    //}
    Eigen::Quaterniond q_from(imu.orientation.w, imu.orientation.x, imu.orientation.y, imu.orientation.z);
    Eigen::Quaterniond q_final = q_from * Eigen::Quaterniond(Rbi);

    sensor_msgs::Imu corrected_imu = imu;
    corrected_imu.linear_acceleration.x = corrected_acc.x();
    corrected_imu.linear_acceleration.y = corrected_acc.y();
    corrected_imu.linear_acceleration.z = corrected_acc.z();
    corrected_imu.angular_velocity.x = corrected_gyo.x();
    corrected_imu.angular_velocity.y = corrected_gyo.y();
    corrected_imu.angular_velocity.z = corrected_gyo.z();
    corrected_imu.orientation.x = q_final.x();
    corrected_imu.orientation.y = q_final.y();
    corrected_imu.orientation.z = q_final.z();
    corrected_imu.orientation.w = q_final.w();

    return corrected_imu;
  };

}
