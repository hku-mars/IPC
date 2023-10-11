#pragma once

#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/kdtree.h>
#include <pcl/search/impl/kdtree.hpp>

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>
#include <geometry_msgs/Vector3.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <Eigen/Eigen>
#include <math.h>
#include <random>
#include "fmt/color.h"
#include "racing_map_generator/config.hpp"


using namespace fmt;
using namespace std;
namespace racing_map_generator {
    class RacingMapGenerator {
    private:
        ros::Subscriber click_sub_;
        ros::NodeHandle nh_;

        ros::Timer pub_timer_;
        ros::Publisher pc_pub_, pc_no_ground_pub_;
        pcl::PointCloud<pcl::PointXYZ> pcl_map_, pcl_no_ground_map_;
        pcl::search::KdTree<pcl::PointXYZ> kdtreeMap;
        sensor_msgs::PointCloud2 pc2_map_, pc2_no_ground_map_;
        vector<int> pointIdxSearch;
        vector<float> pointSquaredDistance;
        bool map_is_ready{false};

        RacingMapConfig cfg_;


        template<class T>
        void LoadParam(string param_name, T &param_value, T default_value) {
            if (nh_.getParam(param_name, param_value)) {
                print("\033[0;32m Load param {} succes: \033[0;0m", param_name.c_str());
                cout << param_value << endl;
            } else {
                print("\033[0;33m Load param {} failed, use default value: \033[0;0m", param_name);
                param_value = default_value;
                cout << param_value << endl;
            }
        }

        template<class T>
        void LoadParam(string param_name, vector<T> &param_value, vector<T> default_value) {
            if (nh_.getParam(param_name, param_value)) {
                print("\033[0;32m Load param {} succes: \033[0;0m", param_name.c_str());
                for (int i = 0; i < param_value.size(); i++) {
                    cout << param_value[i] << " ";
                }
                cout << endl;
            } else {
                print("\033[0;33m Load param {} failed, use default value: \033[0;0m", param_name);
                param_value = default_value;
                for (int i = 0; i < param_value.size(); i++) {
                    cout << param_value[i] << " ";
                }
                cout << endl;
            }
        }

    public:
        RacingMapGenerator(ros::NodeHandle &nh) {
            nh_ = nh;
            cfg_ = RacingMapConfig(nh);
            pc_pub_ = nh_.advertise<sensor_msgs::PointCloud2>("/map_generator/global_cloud", 1);
            pub_timer_ = nh_.createTimer(ros::Duration(1.0 /cfg_.pub_freq), &RacingMapGenerator::PubTimerCallback, this);
            GenerateCircles();
            pcl_map_.width = pcl_map_.points.size();
            pcl_map_.height = 1;
            pcl_map_.is_dense = true;
            pcl::toROSMsg(pcl_map_, pc2_map_);
            pc2_map_.header.frame_id = "world";
            map_is_ready = true;
            map_is_ready = true;
        }

        ~RacingMapGenerator() {}

    private:

        void PubTimerCallback(const ros::TimerEvent &e) {
            if (!map_is_ready) {
                return;
            }
            pc2_map_.header.stamp = ros::Time::now();
            pc_pub_.publish(pc2_map_);
        }

        void GenerateCircles() {
            pcl::PointXYZ pt_random;
            for (auto py: cfg_.pos_yaw_pair) {
                double x0, y0, z0;
                std::vector<Eigen::Vector3d> circle_set;

                x0 = py[0];
                y0 = py[1];
                z0 = py[2];

                double x, y, z;
                Eigen::Vector3d pt3, pt3_rot;
                double rmin = cfg_.size - 0.3;
                double rmax = cfg_.size + 0.3;
                for (double theta = -M_PI; theta < M_PI; theta += 0.025) {
                    for (double R = rmin; R < rmax; R += cfg_.resolution) {
                        x = cos(theta) * R;
                        y = sin(theta) * R;
                        z = 0;
                        pt3 << x, y, z;
                        circle_set.push_back(pt3);
                    }
                }
                // Define a random 3d rotation matrix
                Eigen::Matrix3d Rot = (Eigen::AngleAxisd(py[3], Eigen::Vector3d::UnitZ())
                                      * Eigen::AngleAxisd(-1.5741, Eigen::Vector3d::UnitY())
                                      * Eigen::AngleAxisd(0, Eigen::Vector3d::UnitX())).toRotationMatrix();

                for (auto pt: circle_set) {
                    pt3_rot = Rot * pt;
                    pt_random.x = pt3_rot(0) + x0 + 0.001;
                    pt_random.y = pt3_rot(1) + y0 + 0.001;
                    pt_random.z = pt3_rot(2) + z0 + 0.001;
                    pcl_map_.points.push_back(pt_random);
                }
            }
            print(fg(color::light_sea_green), " [Rings] Generate success.\n");
        }


    };
}
