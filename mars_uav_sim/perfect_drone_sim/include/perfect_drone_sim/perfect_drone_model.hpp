#ifndef _PERFECT_DRONE_SIM_HPP_
#define _PERFECT_DRONE_SIM_HPP_

#include "ros/ros.h"
#include "visualization_msgs/Marker.h"
#include "visualization_msgs/MarkerArray.h"
#include "quadrotor_msgs/PositionCommand.h"
#include "nav_msgs/Odometry.h"
#include "tf2_ros/transform_broadcaster.h"
#include "string"
#include "Eigen/Dense"

typedef Eigen::Matrix<double, 3, 1> Vec3;
typedef Eigen::Matrix<double, 3, 3> Mat33;

typedef Eigen::Matrix<double, 3, 3> StatePVA;
typedef Eigen::Matrix<double, 3, 4> StatePVAJ;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DynamicMat;
typedef Eigen::MatrixX4d MatX4;
typedef std::pair<double, Vec3> TimePosPair;

typedef Eigen::Matrix3Xd PolyhedronV;
typedef Eigen::MatrixX4d PolyhedronH;

class PerfectDrone {
 public:
  PerfectDrone(ros::NodeHandle &n) {
	nh_ = n;
	cmd_sub_ = nh_.subscribe("/planning/pos_cmd", 100, &PerfectDrone::cmdCallback, this);
	odom_pub_ = nh_.advertise<nav_msgs::Odometry>("/lidar_slam/odom", 100);
	robot_pub_ = nh_.advertise<visualization_msgs::Marker>("robot", 100);
	n.param("mesh_resource", mesh_resource_, std::string("package://perfect_drone_sim/meshes/f250.dae"));
	n.param("init_position/x", position_.x(), 0.0);
	n.param("init_position/y", position_.y(), 0.0);
	n.param("init_position/z", position_.z(), 1.0);
	q_ = Eigen::Quaterniond(Mat33::Identity());
	odom_.header.frame_id = "world";
	odom_pub_timer_ = nh_.createTimer(ros::Duration(0.01), &PerfectDrone::publishOdom, this);
  }
  ~PerfectDrone() {
  }

 private:
  ros::Subscriber cmd_sub_;
  ros::Publisher odom_pub_, robot_pub_;
  ros::Timer odom_pub_timer_;
  ros::NodeHandle nh_;
  Vec3 position_, vel_;
  Vec3 pos_r_, acc_r_;
  double yaw_r_;
  Eigen::Quaterniond q_;
  nav_msgs::Odometry odom_;
  std::string mesh_resource_;

	ros::Time last_cmd_time;
	Vec3 pos_last;

  void cmdCallback(const quadrotor_msgs::PositionCommandConstPtr &msg) {
	Vec3 pos(msg->position.x, msg->position.y, msg->position.z);
	Vec3 acc(msg->acceleration.x, msg->acceleration.y, msg->acceleration.z);
	double yaw = msg->yaw;
	// pos_r_ = pos;
	position_ = pos;
	acc_r_ = acc;
	yaw_r_ = yaw;
	// updateFlatness(pos, acc, yaw);
  }

  void publishOdom(const ros::TimerEvent &e) {
	updateFlatness(position_, acc_r_, yaw_r_);

	odom_.pose.pose.position.x = position_.x();
	odom_.pose.pose.position.y = position_.y();
	odom_.pose.pose.position.z = position_.z();
	odom_.twist.twist.linear.x = vel_.x();
	odom_.twist.twist.linear.y = vel_.y();
	odom_.twist.twist.linear.z = vel_.z();

	odom_.pose.pose.orientation.x = q_.x();
	odom_.pose.pose.orientation.y = q_.y();
	odom_.pose.pose.orientation.z = q_.z();
	odom_.pose.pose.orientation.w = q_.w();

	odom_.header.stamp = ros::Time::now();

	odom_pub_.publish(odom_);

	static tf2_ros::TransformBroadcaster br_map_ego;
	geometry_msgs::TransformStamped transformStamped;
	transformStamped.header.stamp = odom_.header.stamp;
	transformStamped.header.frame_id = "world";
	transformStamped.child_frame_id = "drone";
	transformStamped.transform.translation.x = odom_.pose.pose.position.x;
	transformStamped.transform.translation.y = odom_.pose.pose.position.y;
	transformStamped.transform.translation.z = odom_.pose.pose.position.z;
	transformStamped.transform.rotation.x = odom_.pose.pose.orientation.x;
	transformStamped.transform.rotation.y = odom_.pose.pose.orientation.y;
	transformStamped.transform.rotation.z = odom_.pose.pose.orientation.z;
	transformStamped.transform.rotation.w = odom_.pose.pose.orientation.w;
	br_map_ego.sendTransform(transformStamped);

	visualization_msgs::Marker meshROS;
	meshROS.header.frame_id = "world";
	meshROS.header.stamp = odom_.header.stamp;
	meshROS.ns = "mesh";
	meshROS.id = 0;
	meshROS.type = visualization_msgs::Marker::MESH_RESOURCE;
	meshROS.action = visualization_msgs::Marker::ADD;
	meshROS.pose.position= odom_.pose.pose.position;
	meshROS.pose.orientation = odom_.pose.pose.orientation;
	meshROS.scale.x = 1;
	meshROS.scale.y = 1;
	meshROS.scale.z = 1;
	meshROS.mesh_resource = mesh_resource_;
	meshROS.mesh_use_embedded_materials = true;
	meshROS.color.a=1.0;
	meshROS.color.r=1.0;
	meshROS.color.g=1.0;
	meshROS.color.b=1.0;
	robot_pub_.publish(meshROS);

  }

  void updateFlatness(const Vec3 &pos, const Vec3 &acc, const double yaw) {
	Vec3 gravity_ = 9.80 * Eigen::Vector3d(0, 0, 1);
	position_ = pos;
	// double dt = 0.01;
	double dt = (ros::Time::now() - last_cmd_time).toSec();
	last_cmd_time = ros::Time::now();
	vel_ = (position_ - pos_last) / dt;
	pos_last = position_;
	double a_T = (gravity_ + acc).norm();
	Eigen::Vector3d xB, yB, zB;
	Eigen::Vector3d xC(cos(yaw), sin(yaw), 0);

	zB = (gravity_ + acc).normalized();
	yB = ((zB).cross(xC)).normalized();
	xB = yB.cross(zB);
	Eigen::Matrix3d R;
	R << xB, yB, zB;
	q_ = Eigen::Quaterniond(R);
  }

};

#endif