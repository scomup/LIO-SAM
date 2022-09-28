#pragma once

#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/geometry/Pose3.h>

namespace gtsam
{
  class GPSFactor : public NoiseModelFactor1<Pose3>
  {
  protected:
    gtsam::Pose3 Tbg_; //b:body g:gps
    gtsam::Point3 twg_; //translation of gps in world
    typedef NoiseModelFactor1<Pose3> Base;

  public:
    virtual ~GPSFactor() {}

    GPSFactor(
        const Key &x_key,
        const gtsam::Pose3 &Tbg,
        const gtsam::Point3 &gps,
        const SharedNoiseModel &noiseModel)
        : Base(noiseModel, x_key),
        Tbg_(Tbg),
        twg_(gps)
        {};

    Vector evaluateError(const Pose3 &Twb,
                         boost::optional<Matrix &> H = boost::none) const
    {
      Matrix66 dTwgdTwb;
      gtsam::Pose3 Twg = Twb.compose(Tbg_, H ? &dTwgdTwb : 0);

      Matrix36 drdTwg;
      Eigen::Vector3d r = Twg.translation(drdTwg) - twg_;

      if (H)
      {
        *H = drdTwg * dTwgdTwb;
      }

     return r;
    };

  };
}// gtsam
