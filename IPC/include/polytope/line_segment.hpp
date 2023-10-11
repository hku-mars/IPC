#ifndef _GEOMETRY_UTILS_LINE_SEGMENT_
#define _GEOMETRY_UTILS_LINE_SEGMENT_


/// For utils
#include "common_type_name.hpp"
#include "geo_utils.hpp"

/// For ROS
#include "ros/ros.h"
#include "visualization_msgs/MarkerArray.h"

class LineSegment {
 private:
	bool undefined_ = true;
	Vec3f p1_, p2_, center_;
	double length_;
 public:

  LineSegment() {}

  LineSegment(const Vec3f & p1, const Vec3f & p2): p1_(p1),p2_(p2) {
	undefined_ = false;
  }

  bool empty(){
	return undefined_;
  }

  double getLength(){
	length_ = (p1_-p2_).norm();
	return length_;
  }

  Vec3f getCenter(){
	center_ = (p1_+p2_)/2;
	return center_;
  }


  void Visualize(const ros::Publisher &pub,
				 string ns = "line_segment",
				 Color point_color = Color::Blue(),
				 decimal_t pt_size = 0.1,
				 Color line_color = Color::Blue(),
				 decimal_t line_size = 0.05) {
	pt_size = max(pt_size, 0.05);
	visualization_msgs::MarkerArray mkr_arr;
	// static int point_id = 0;
	// static int line_cnt = 0;
	int point_id = 0;
	int line_cnt = 0;
	/* Publish point */
	visualization_msgs::Marker point;
	point.header.frame_id = "world";
	point.header.stamp = ros::Time::now();
	point.ns = ns.c_str();
	point.id = point_id++;
	point.action = visualization_msgs::Marker::ADD;
	point.pose.orientation.w = 1.0;
	point.type = visualization_msgs::Marker::SPHERE;
	// LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
	point.scale.x = pt_size;
	point.scale.y = pt_size;
	point.scale.z = pt_size;
	// Line list is blue
	point.color = point_color;
	point.color.a = 0.3;
	// Create the vertices for the points and lines
	geometry_msgs::Point p;
	p.x = p1_.x();
	p.y = p1_.y();
	p.z = p1_.z();
	point.pose.position = p;
	mkr_arr.markers.push_back(point);

	p.x = p2_.x();
	p.y = p2_.y();
	p.z = p2_.z();
	point.pose.position = p;
	point.id = point_id++;
	mkr_arr.markers.push_back(point);

// publish cylinder
	visualization_msgs::Marker line_list;
	line_list.header.frame_id = "world";
	line_list.header.stamp = ros::Time::now();
	line_list.ns = ns + "line";
	line_list.id = line_cnt++;
	line_list.action = visualization_msgs::Marker::ADD;
	line_list.type = visualization_msgs::Marker::CYLINDER;

	line_list.scale.x = line_size;
	line_list.scale.y = line_size;
	line_list.scale.z = getLength();
	// Line list is blue
	line_list.color = line_color;
	line_list.color.a = 0.3;
	Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), (p2_ - p1_));
	// Create the vertices for the points and lines
	line_list.pose.orientation.w = q.w();
	line_list.pose.orientation.x = q.x();
	line_list.pose.orientation.y = q.y();
	line_list.pose.orientation.z = q.z();
	Vec3f center = getCenter();
	line_list.pose.position.x = center.x();
	line_list.pose.position.y = center.y();
	line_list.pose.position.z = center.z();
	line_list.points.push_back(p);
	mkr_arr.markers.push_back(line_list);
	pub.publish(mkr_arr);
  }

};

typedef vector<LineSegment> LineSegments;

#endif