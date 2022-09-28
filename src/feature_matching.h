#pragma once

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/common.h>
#include <pcl/common/transforms.h>
#include <pcl/common/angles.h>
#include <pcl/filters/voxel_grid.h>

//#include <pcl_conversions/pcl_conversions.h>

namespace GRAPH_SLAM
{
  using PointType = pcl::PointXYZI;
  using PointCloud = pcl::PointCloud<PointType>;
  using PointCloudPtr = PointCloud::Ptr;
  using PointCloudConstPtr = PointCloud::ConstPtr;

  class FeatureMatching
  {
  public:
    FeatureMatching();

    void setInputTarget(const PointCloudPtr surf);

    void setInputSource(const PointCloudPtr surf);

    void align(const Eigen::Affine3d &guess);

    inline Eigen::Affine3d getFinalTransformation() const { return T_; }

  private:
    bool findPlane(PointCloudPtr scan,
                   std::vector<int> &index,
                   Eigen::Vector4d &plane) const;


    void reset();

    void surfGradientHessian(Eigen::Matrix<double, 6, 1> &g, Eigen::Matrix<double, 6, 6> &H);

    PointCloudPtr scan_surf_;
    PointCloudPtr tar_surf_;

    pcl::KdTreeFLANN<PointType>::Ptr surf_kdtree_;

    Eigen::Affine3d T_;
    //pcl::VoxelGrid<PointType> dsfilter_;

    Eigen::Matrix<double, 3, 6> calcdTadx(
        const Eigen::Matrix<double, 3, 3> &T,
        const Eigen::Vector3d &a) const;
  };

}