#include "feature_matching.h"
#include "common/math.h"

namespace GRAPH_SLAM
{
  Eigen::Matrix<double, 3, 6> FeatureMatching::calcdTadx(
      const Eigen::Matrix<double, 3, 3> &T,
      const Eigen::Vector3d &a) const
  {
    Eigen::Matrix<double, 3, 6> dTadx;
    dTadx << T, T * COMMON::skew<double>(-a);
    return dTadx;
  }

  bool FeatureMatching::findPlane(PointCloudPtr scan,
                                  std::vector<int> &index,
                                  Eigen::Vector4d &plane) const
  {
    int n = index.size();
    Eigen::Matrix<double, -1, 3> A(n, 3);
    Eigen::Matrix<double, -1, 1> b(n, 1);

    b.fill(-1);

    for (int j = 0; j < n; j++)
    {
      A.row(j) = scan->points[index[j]].getVector3fMap().template cast<double>();
    }

    Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);

    plane.head(3) = x;
    plane(3) = 1;
    plane /= x.norm();

    bool valid = true;
    for (int j = 0; j < n; j++)
    {
      Eigen::Vector3d pt = scan->points[index[j]].getVector3fMap().template cast<double>();

      double dist = fabs(pt.dot(plane.head(3)) + plane(3));

      if (dist > 0.2)
      {
        valid = false;
        break;
      }
    }

    return valid;
  }

  void FeatureMatching::setInputTarget(const PointCloudPtr surf)
  {
    reset();
    tar_surf_ = surf;
    surf_kdtree_->setInputCloud(tar_surf_);

  }

  void FeatureMatching::setInputSource(const PointCloudPtr surf)
  {
    scan_surf_ = surf;
  }

  FeatureMatching::FeatureMatching()
  {
  }

  void FeatureMatching::reset()
  {
    scan_surf_.reset(new pcl::PointCloud<PointType>());
    tar_surf_.reset(new pcl::PointCloud<PointType>());
    surf_kdtree_.reset(new pcl::KdTreeFLANN<PointType>());
  }

  void FeatureMatching::surfGradientHessian(Eigen::Matrix<double, 6, 1> &g, Eigen::Matrix<double, 6, 6> &H)
  {

    PointCloudPtr transformed_scan(new pcl::PointCloud<PointType>());
    pcl::transformPointCloud(*scan_surf_, *transformed_scan, T_);
    int n = scan_surf_->size();

    std::vector<Eigen::Matrix<double, 6, 1>> gs(n, Eigen::Matrix<double, 6, 1>::Zero());
    std::vector<Eigen::Matrix<double, 6, 6>> Hs(n, Eigen::Matrix<double, 6, 6>::Zero());

    #pragma omp parallel for num_threads(8)
    for (int i = 0; i < n; i++)
    {

      PointType &p = scan_surf_->points[i];
      PointType &p_star = transformed_scan->points[i];

      std::vector<int> index;
      std::vector<float> dist;
      surf_kdtree_->nearestKSearch(p_star, 5, index, dist);

      if (dist[4] >= 1.0)
        continue;

      Eigen::Vector4d plane;
      bool vaild = findPlane(tar_surf_, index, plane);

      if (!vaild)
        continue;

      Eigen::Vector3d a(p.x, p.y, p.z);

      double r = plane(0) * p_star.x + plane(1) * p_star.y + plane(2) * p_star.z + plane(3);

      double w = 1 - 0.9 * fabs(r) / sqrt(sqrt(a.dot(a)));

      if (w < 0.1)
        continue;

      r = r * w;

      double e = r * r;
      auto rho = COMMON::huberKernel(e, 0.3);

      Eigen::Matrix<double, 3, 6> dTadx = calcdTadx(T_.rotation(), a);
      Eigen::Matrix<double, 1, 6> J;
      Eigen::Vector3d dDdTa = plane.head(3) * w;

      J = dDdTa.transpose() * dTadx;
      Hs[i] = J.transpose() * J * rho(1);
      gs[i] = J.transpose() * r * rho(1);
    }
    for (int i = 0; i < n; i++)
    {
      H += Hs[i];
      g += gs[i];
    }
  }

  void FeatureMatching::align(const Eigen::Affine3d &guess)
  {
    //std::cout<<"scan num:"<< scan_->size()<<std::endl;
    //std::cout<<"tar num:"<<   tar_->size()<<std::endl;

    T_ = guess;

    for (int iterCount = 0; iterCount < 30; iterCount++)
    {
      Eigen::Matrix<double, 6, 1> g;
      Eigen::Matrix<double, 6, 6> H;

      g.setZero();
      H.setZero();

      surfGradientHessian(g, H);

      Eigen::Matrix<double, 6, 1> dx = H.colPivHouseholderQr().solve(-g);
      T_ = T_ * COMMON::param2T<double>(dx);

      float deltaR = sqrt(
          pow(pcl::rad2deg(dx(0)), 2) +
          pow(pcl::rad2deg(dx(1)), 2) +
          pow(pcl::rad2deg(dx(2)), 2));
      float deltaT = sqrt(
          pow(dx(3) * 100, 2) +
          pow(dx(4) * 100, 2) +
          pow(dx(5) * 100, 2));

      if (deltaR < 0.05 && deltaT < 0.05) // 0.05cm and 0.05 degree
      {
        break; // converged
      }
    }
  }

}
