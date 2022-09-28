#pragma once

#include <gtsam/nonlinear/NonlinearFactor.h>

namespace gtsam
{
  class GravityFactor : public NoiseModelFactor1<Pose3>
  {
  protected:
    Eigen::Vector3d acc_;
    typedef NoiseModelFactor1<Pose3> Base;

  public:
    virtual ~GravityFactor() {}

    GravityFactor(
        const Key &x_key,
        const Eigen::Vector3d &acc,
        const SharedNoiseModel &noiseModel)
        : Base(noiseModel, x_key),
        acc_(acc)
        {};

    Vector evaluateError(const Pose3 &x,
                         boost::optional<Matrix &> H = boost::none) const
    {
      //r = T(x_wb,acc_b) cross g
      Matrix33 daccwdx;
      Eigen::Vector3d acc_w = x.rotation().rotate(acc_, H ? &daccwdx : 0);
      Eigen::Vector3d g =  Eigen::Vector3d(0,0,9.8);
      Eigen::Vector3d r = acc_w.cross(g);
      if (H)
      {
        Matrix33 drdaccw = (Matrix33() << 0.0, -g(2), +g(1), +g(2), 0.0, -g(0), -g(1), +g(0), 0.0).finished();
        Matrix36 tmpH = Matrix36::Zero();
        tmpH.block(0,0,3,3) = -drdaccw*daccwdx;
        *H = tmpH;
      }
      return r;
    };

  };
}// gtsam
