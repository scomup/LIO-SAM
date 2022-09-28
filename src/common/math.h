/**
 * @file   math.h
 * @brief  演算処理をまとめたライブラリ
 * @author Panasonic Advanced Technology Development Co.,Ltd.
 * @date   2022/02/09
 */
#ifndef NDTMATCHING_MATH_H_
#define NDTMATCHING_MATH_H_

#include <Eigen/Dense>
#include <math.h>

namespace COMMON
{

#define EPSILON 1e-5


  Eigen::Vector2d gaussianKernel(double e2, double d = 1)
  {
    Eigen::Vector2d rho;
    double t = exp(-e2 * d);
    rho(0) = t;
    rho(1) = -d * t;
    //rho(2) = d * rho(1);
    return rho;
  }

  Eigen::Vector2d CauchyKernel(double e2, double delta = 2)
  {
    Eigen::Vector2d rho;
    double c2 = delta * delta;
    double c2_inv = 1. / c2;
    double aux = c2_inv * e2 + 1.0;
    rho(0) = (c2 * std::log(aux));
    rho(1) = (1. / aux);
    //rho(2) = -(c2_inv * rho[1] * rho[1]);
    return rho;
  }

  Eigen::Vector2d huberKernel(double e, double delta = 2)
  {
    Eigen::Vector2d rho;
    double dsqr = delta * delta;
    if (e <= dsqr)
    {
      rho[0] = e;
      rho[1] = 1.;
      //rho[2] = 0.;
    }
    else
    {
      double sqrte = sqrt(e);
      rho[0] =
          2 * sqrte * delta - dsqr;
      rho[1] = delta / sqrte;
      //rho[2] = -0.5 * rho[1] / e;
    }
    return rho;
  }

  template <typename T>
  Eigen::Matrix<T, 3, 1> R2rot(const Eigen::Matrix<T, 3, 3> &R)
  {
    const Eigen::Quaternion<T> q(R);
    const T squared_n = q.vec().dot(q.vec());
    T two_atan_nbyw_by_n;
    if (squared_n < EPSILON * EPSILON)
    {
      const T squared_w = q.w() * q.w();
      two_atan_nbyw_by_n = 2. / q.w() - (2.0 / 3.0) * (squared_n) / (q.w() * squared_w);
    }
    else
    {
      const T n = sqrt(squared_n);
      if (abs(q.w()) < EPSILON)
      {
        if (q.w() > 0.)
          two_atan_nbyw_by_n = M_PI / n;
        else
          two_atan_nbyw_by_n = -M_PI / n;
      }

      else
      {
        two_atan_nbyw_by_n = 2. * std::atan(n / q.w()) / n;
      }
    }
    return two_atan_nbyw_by_n * q.vec();
  }

  template <typename T>
  Eigen::Matrix<T, 3, 3> rot2R(const Eigen::Matrix<T, 3, 1> &v)
  {
    const T theta_sq = v.dot(v);

    T imag_factor = 0.;
    T real_factor = 0.;
    if (theta_sq < EPSILON * EPSILON)
    {
      const T theta_po4 = theta_sq * theta_sq;
      imag_factor = 0.5 - (1.0 / 48.0) * theta_sq + (1.0 / 3840.0) * theta_po4;
      real_factor = 1. - (1.0 / 8.0) * theta_sq + (1.0 / 384.0) * theta_po4;
    }
    else
    {
      const T theta = std::sqrt(theta_sq);
      const T half_theta = 0.5 * theta;
      const T sin_half_theta = std::sin(half_theta);
      imag_factor = sin_half_theta / theta;
      real_factor = std::cos(half_theta);
    }

    Eigen::Quaternion<T> q(real_factor, imag_factor * v.x(), imag_factor * v.y(), imag_factor * v.z());

    return q.toRotationMatrix();
  }

  template <typename T>
  Eigen::Matrix<T, 4, 4> param2T(const Eigen::Matrix<T, 6, 1> &p)
  {
    Eigen::Matrix<T, 4, 4> Trans = (Eigen::Translation<T, 3>(p(0), p(1), p(2)) * rot2R<T>(p.tail(3))).matrix();
    return Trans;
  }

  template <typename T>
  Eigen::Matrix<T, 6, 1> T2param(const Eigen::Matrix<T, 4, 4> &mat)
  {
    Eigen::Transform<T, 3, Eigen::Affine, Eigen::ColMajor> eig_transformation;
    eig_transformation.matrix() = mat;
    Eigen::Matrix<T, 6, 1> x;
    Eigen::Matrix<T, 3, 1> init_translation = eig_transformation.translation();
    Eigen::Matrix<T, 3, 1> init_rotation = R2rot<T>(eig_transformation.rotation());
    x << init_translation(0), init_translation(1), init_translation(2), init_rotation(0), init_rotation(1), init_rotation(2);
    return x;
  }

  template <typename T>
  inline Eigen::Matrix<T, 3, 3> skew(const Eigen::Matrix<T, 3, 1> &w)
  {
    return (Eigen::Matrix<T, 3, 3>() << 0.0, -w[2], +w[1], +w[2], 0.0, -w[0], -w[1], +w[0], 0.0).finished();
  }
} // namespace NDT

#endif
