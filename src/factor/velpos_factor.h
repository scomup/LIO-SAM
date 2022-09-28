#pragma once

#include <gtsam/geometry/Pose3.h>
#include <gtsam/nonlinear/NonlinearFactor.h>

namespace gtsam
{
  class VelposFactor : public NoiseModelFactor3<Pose3, Pose3, Eigen::Vector3d>
  {
  protected:
    double dt_;
    typedef gtsam::NoiseModelFactor3<Pose3, Pose3, Eigen::Vector3d> Base;

  public:
    virtual ~VelposFactor() {}

    VelposFactor(
        const Key &xi_key, const Key &xj_key, const Key &vi_key,
        const double &dt,
        const SharedNoiseModel &noiseModel)
        : Base(noiseModel, xi_key, xj_key, vi_key),
          dt_(dt)
    {
    }

    Vector evaluateError(const Pose3 &x_wbi, const Pose3 &x_wbj, const Eigen::Vector3d &vel,
                         boost::optional<Matrix &> H1 = boost::none,
                         boost::optional<Matrix &> H2 = boost::none,
                         boost::optional<Matrix &> H3 = boost::none) const
    {

      double dt_inv = 1. / dt_;
      Eigen::Vector3d r = vel - (x_wbj.translation() - x_wbi.translation()) * dt_inv;

      if (H1)
      {
        Matrix36 tmpH1 = Matrix36::Zero();
        tmpH1.block(0,3,3,3) = x_wbi.rotation().matrix() * dt_inv;

        *H1 = tmpH1;
      }
      if (H2)
      {
        Matrix36 tmpH2 = Matrix36::Zero();
        tmpH2.block(0,3,3,3) = -x_wbj.rotation().matrix() * dt_inv;
        *H2 = tmpH2;
      }

      if (H3)
      {
        *H3 = x_wbi.rotation().matrix();
      }
      return r;
    };
  };
}// gtsam
