
#ifndef _GEOMETRY_UTILS_ELLIPSOID_HPP_
#define _GEOMETRY_UTILS_ELLIPSOID_HPP_

#include "ros/ros.h"
/// ros msgs
#include "visualization_msgs/MarkerArray.h"

#include "Eigen/Dense"
#include "vector"

/// utils
#include "common_type_name.hpp"
#include "geo_utils.hpp"

class Ellipsoid {
 private:
  /// If the ellipsoid is empty
  bool undefined{true};

  /// The ellipsoid is defined by shape C and center d
  Mat3f C_, C_inv_;
  Mat3f R_;
  Vec3f d_, r_;

 public:
  Ellipsoid() {}

  Ellipsoid(const Mat3f &C, const Vec3f &d) : C_(C), d_(d) {
	undefined = false;
	C_inv_ = C_.inverse();

	Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::FullPivHouseholderQRPreconditioner> svd(C_, Eigen::ComputeFullU);
	Eigen::Matrix3d U = svd.matrixU();
	Eigen::Vector3d S = svd.singularValues();
	if (U.determinant() < 0.0) {
	  R_.col(0) = U.col(1);
	  R_.col(1) = U.col(0);
	  R_.col(2) = U.col(2);
	  r_(0) = S(1);
	  r_(1) = S(0);
	  r_(2) = S(2);
	} else {
	  R_ = U;
	  r_ = S;
	}
  }

  Ellipsoid(const Mat3f &R, const Vec3f &r, const Vec3f &d) : R_(R), r_(r), d_(d) {
	undefined = false;
	C_ = R_ * r_.asDiagonal() * R_.transpose();
	C_inv_ = C_.inverse();
  }

  /// If this ellipsoid is empty
  bool empty(){
	return undefined;
  }

  /// Find the closestPoint in a point set
  int nearestPointId(const Eigen::Matrix3Xd &pc) {
	Eigen::VectorXd dists = (C_inv_ * (pc.colwise() - d_)).colwise().norm();
	int np_id;
	double np_dist = dists.minCoeff(&np_id);
	return np_id;
  }

  /// Find the closestPoint in a point set
  Vec3f nearestPoint(const Eigen::Matrix3Xd &pc) {
	Eigen::VectorXd dists = (C_inv_ * (pc.colwise() - d_)).colwise().norm();
	int np_id;
	double np_dist = dists.minCoeff(&np_id);
	return pc.col(np_id);
  }

  /// Find the closestPoint in a point set
  double nearestPointDis(const Eigen::Matrix3Xd &pc, int &np_id) {
	Eigen::VectorXd dists = (C_inv_ * (pc.colwise() - d_)).colwise().norm();
	double np_dist = dists.minCoeff(&np_id);
	return np_dist;
  }

  /// Get the shape of the ellipsoid
  Mat3f C() {
	return C_;
  }

  /// Get the center of the ellipsoid
  Vec3f d() {
	return d_;
  }

  Mat3f R() {
	return R_;
  }

  Vec3f r() {
	return r_;
  }

  /// Convert a point to the ellipsoid frame
  Vec3f toEllipsoidFrame(const Vec3f &pt_w) {
	return C_inv_ * (pt_w - d_);
  }

  /// Convert a set of points to the ellipsoid frame
  Eigen::Matrix3Xd toEllipsoidFrame(const Eigen::Matrix3Xd &pc_w) {
	return C_inv_ * (pc_w.colwise() - d_);
  }

  /// Convert a point to the world frame
  Vec3f toWorldFrame(const Vec3f &pt_e) {
	return C_ * pt_e + d_;
  }

  /// Convert a set of points to the world frame
  Eigen::Matrix3Xd toWorldFrame(const Eigen::Matrix3Xd &pc_e) {
	return (C_ * pc_e).colwise() + d_;
  }

  /// Convert a plane to the ellipsoid frame
  Eigen::Vector4d toEllipsoidFrame(const Eigen::Vector4d &plane_w) {
	Eigen::Vector4d plane_e;
	plane_e.head(3) = plane_w.head(3).transpose() * C_;
	plane_e(3) = plane_w(3) + plane_w.head(3).dot(d_);
	return plane_e;
  }

  /// Convert a plane to the ellipsoid frame
  Eigen::Vector4d toWorldFrame(const Eigen::Vector4d &plane_e) {
	Eigen::Vector4d plane_w;
	plane_w.head(3) = plane_e.head(3).transpose() * C_inv_;
	plane_w(3) = plane_e(3) - plane_w.head(3).dot(d_);
	return plane_w;
  }

  /// Convert a set of planes to ellipsoid frame
  Eigen::MatrixX4d toEllipsoidFrame(const Eigen::MatrixX4d &planes_w) {
	Eigen::MatrixX4d planes_e(planes_w.rows(), planes_w.cols());
	planes_e.leftCols(3) = planes_w.leftCols(3) * C_;
	planes_e.rightCols(1) = planes_w.rightCols(1) + planes_w.leftCols(3) * d_;
	return planes_e;
  }

  /// Convert a set of planes to ellipsoid frame
  Eigen::MatrixX4d toWorldFrame(const Eigen::MatrixX4d &planes_e) {
	Eigen::MatrixX4d planes_w(planes_e.rows(), planes_e.cols());
	planes_w.leftCols(3) = planes_e.leftCols(3) * C_inv_;
	planes_w.rightCols(1) = planes_e.rightCols(1) + planes_w.leftCols(3) * d_;
	return planes_w;
  }

  /// Calculate the distance of a point in world frame
  double dist(const Vec3f &pt_w) const {
	return (C_inv_ * (pt_w - d_)).norm();
  }

  /// Calculate the distance of a point in world frame
  Eigen::VectorXd dist(const Eigen::Matrix3Xd &pc_w) const {
	return (C_inv_ * (pc_w.colwise() - d_)).colwise().norm();
  }

  bool noPointsInside(vector<Vec3f> &pc, const Eigen::Matrix3d R,
					  const Vec3f &r, const Vec3f &p) {
	Eigen::Matrix3d C_inv;
	C_inv = r.cwiseInverse().asDiagonal() * R.transpose();
	for (auto pt_w : pc) {
	  double d = (C_inv * (pt_w - p)).norm();
	  if (d <= 1) {
		return false;
	  }
	}
	return true;
  }

  bool pointsInside(const Eigen::Matrix3Xd &pc,
					Eigen::Matrix3Xd &out,
					int &min_pt_id) {
	Eigen::VectorXd vec = (C_inv_ * (pc.colwise() - d_)).colwise().norm();
	vector<Vec3f> pts;
	pts.reserve(pc.cols());
	int cnt = 0;
	min_pt_id = 0;
	double min_dis = numeric_limits<double>::max();
	for (int i = 0; i < vec.size(); i++) {
	  if (vec(i) <= 1) {
		pts.push_back(pc.col(i));
		if (vec(i) <= min_dis) {
		  min_pt_id = cnt;
		  min_dis = vec(i);
		}
		cnt++;
	  }
	}
	if(!pts.empty()){
	  out = Eigen::Map<const Eigen::Matrix<double, 3, -1, Eigen::ColMajor>>(pts[0].data(), 3, pts.size());
	  return true;
	}else{
	  return false;
	}
  }

  /// Check if the point is inside, non-exclusive
  bool inside(const Vec3f &pt) const {
	return dist(pt) <= 1;
  }

  void Visualize(const ros::Publisher &pub,
				 Color color = Color::Orange()) {
	visualization_msgs::Marker mkr;
	Eigen::Quaterniond q(R_);
	// static int id = 0;
	int id = 0;
	mkr.id = id++;
	mkr.type = visualization_msgs::Marker::SPHERE;
	mkr.header.frame_id = "world";
	mkr.header.stamp = ros::Time::now();
	mkr.ns = "ellp";
	mkr.id = id++;
	mkr.action = visualization_msgs::Marker::ADD;
	mkr.pose.orientation.w = q.w();
	mkr.pose.orientation.x = q.x();
	mkr.pose.orientation.y = q.y();
	mkr.pose.orientation.z = q.z();
	mkr.pose.position.x = d_.x();
	mkr.pose.position.y = d_.y();
	mkr.pose.position.z = d_.z();
	mkr.scale.x = r_.x() * 2;
	mkr.scale.y = r_.y() * 2;
	mkr.scale.z = r_.z() * 2;
	mkr.color = color;
	mkr.color.a = 0.3;
	visualization_msgs::MarkerArray mkr_arr;
	mkr_arr.markers.push_back(mkr);
	pub.publish(mkr_arr);
  }

};

#endif