//
// Created by yunfan on 2022/2/17.
//

#ifndef _DEPTH_GRID_HPP_
#define _DEPTH_GRID_HPP_

#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl_conversions/pcl_conversions.h>
#include "pcl_ros/point_cloud.h"
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/kdtree.h>
#include <pcl/search/impl/kdtree.hpp>
#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>
#include <nav_msgs/Odometry.h>
#include <Eigen/Eigen>
#include "aros_util/transport_util.h"
#include <opencv2/opencv.hpp>
#include "pthread.h"
#include<pcl/common/transforms.h>
#define TINYCOLORMAP_WITH_EIGEN
#include "tinycolormap.hpp"
#include "fmt/color.h"
using namespace std;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Matrix3d Mat33;
using namespace Eigen;
using namespace fmt;

class DepthGrid {

private:
    template<class T>
    bool LoadParam(string param_name, T &param_value, T default_value) {
        if (nh_.getParam(param_name, param_value)) {
            printf("\033[0;32m Load param %s succes: \033[0;0m", param_name.c_str());
            cout << param_value << endl;
            return true;
        } else {
            printf("\033[0;33m Load param %s failed, use default value: \033[0;0m", param_name.c_str());
            param_value = default_value;
            cout << param_value << endl;
            return false;
        }
    }

public:
    DepthGrid() {}

    DepthGrid(ros::NodeHandle nh) {
        init(nh);
    }

    bool IsFreeInKnownSpace(const Vector3d & world_pt){
        Vec3 fov_pt = toFovFrame(world_pt);
        double pt_depth;
        Vector2i grid_id = toGridIndex(fov_pt, pt_depth);
        if(!IsGridIdValid(grid_id)){
            // not in FOV
            return false;
        }
        double map_dep = depth_map_[grid_id.x()][grid_id.y()];
        if(map_dep == 0){map_dep = cfg_.sensing_range;}
        if(map_dep > pt_depth || pt_depth > cfg_.sensing_range){
            return  true;
        }
        return false;
    }

    typedef shared_ptr<DepthGrid> Ptr;

    void init(ros::NodeHandle &nh) {
        nh_ = nh;
        print(fg(color::sea_shell), " -- [DepthGrid] Init process\n");
        // Default resolution
        LoadParam("/depth_grid/sensing_range", cfg_.sensing_range, 15.0);
        LoadParam("/depth_grid/resolution", cfg_.resolution, 0.1);
        LoadParam("/depth_grid/vertical_fov", cfg_.vertical_fov, 70.0);
        LoadParam("/depth_grid/horizon_fov", cfg_.horizon_fov, 70.0);
        string cloud_topic, odom_topic;
        LoadParam("/depth_grid/cloud_topic", cloud_topic, string("/pcl_render_node/cloud"));
        LoadParam("/depth_grid/odom_topic", odom_topic, string("/lidar_slam/odom"));
        LoadParam("/depth_grid/min_radius", cfg_.min_radius, 0.3);
        LoadParam("/depth_grid/max_radius", cfg_.max_radius, 5.0);
        LoadParam("/depth_grid/point_filt_num", cfg_.point_filt_num, 1);
        LoadParam("/depth_grid/bisec_num", cfg_.bisec_num, 10);
        if (LoadParam("/depth_grid/angular_reso", cfg_.angular_reso, 0.1)) {
            print(fg(color::gold), " -- [DepthGrid] Use angular resolution\n");
        } else {
            cfg_.angular_reso = cfg_.resolution / cfg_.sensing_range;
            print(fg(color::sea_green), " -- [DepthGrid] Use distance resolution, angle resolution: {}.\n",
                  cfg_.angular_reso);
        }
        // Conver deg FOV to rad.
        cfg_.vertical_fov = cfg_.vertical_fov / 2 * M_PI / 180;
        cfg_.horizon_fov = cfg_.horizon_fov / 2 * M_PI / 180;

        // Calculate Depth Grid Size
        cfg_.grid_size_hor = static_cast<int>(2 * cfg_.horizon_fov / cfg_.angular_reso);
        cfg_.grid_size_ver = static_cast<int>(2 * cfg_.vertical_fov / cfg_.angular_reso);

        //        cfg_.bisec_num = static_cast<int>(log(1 / cfg_.resolution) / log(2) * 1.5);
        print(fg(color::lime_green), " -- [DepthGrid] Init success.\n");
        print(fg(color::azure), "\tgrid_size_hor = {}\n", cfg_.grid_size_hor);
        print(fg(color::azure), "\tgrid_size_ver = {}\n", cfg_.grid_size_ver);
        print(fg(color::azure), "\tbisec_num = {}\n", cfg_.bisec_num);


        depth_map_ = new double *[cfg_.grid_size_ver];
        for (int i = 0; i < cfg_.grid_size_ver; i++) {
            depth_map_[i] = new double[cfg_.grid_size_hor];
        }




        local_cloud_sub_ = nh_.subscribe(cloud_topic, 10, &DepthGrid::LocalCloudCallback, this);
        odom_sub_ = nh.subscribe(odom_topic, 10, &DepthGrid::OdomCallback, this);
        depth_img_pub_ = nh.advertise<sensor_msgs::Image>("depth", 10);
        depth_pc_pub_ = nh_.advertise<sensor_msgs::PointCloud2>("depth_pc", 10);
    }

    ~DepthGrid() {
        delete[] depth_map_;

    }

    bool IsGridIdValid(Vector2i &id) {
        if (id.x() < 0 || id.x() >= cfg_.grid_size_ver) return false;
        if (id.y() < 0 || id.y() >= cfg_.grid_size_hor) return false;
        return true;
    }


private:
    // ros
    ros::NodeHandle nh_;
    ros::Subscriber odom_sub_;
    ros::Subscriber local_cloud_sub_;
    ros::Publisher depth_img_pub_, mkr_pub_, mkr_arr_pub_, depth_pc_pub_;
    // Data
    double **depth_map_;
    bool has_odom_;
    struct Config {
        int grid_size_ver, grid_size_hor;
        int bisec_num, point_filt_num;
        double resolution, angular_reso;
        double sensing_range;
        double vertical_fov;
        double horizon_fov;
        double min_radius, max_radius;
    } cfg_;

    struct CamPose {
        Vec3 position;
        Mat33 rotation;
        nav_msgs::Odometry odom;
    } cam_pose_;

private:
    template<typename T>
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    Vector3d toFovFrame(const Vector3d &world) {
        Vec3 fov = cam_pose_.rotation.inverse() * (world - cam_pose_.position);
        return fov;
    }

    Vector3d toWorldFrame(const Vector3d &fov) {
        Vec3 world = cam_pose_.rotation * fov + cam_pose_.position;
        return world;
    }

    Vector2i toGridIndex(const Vector3d &fov_pt, double &depth) {
        depth = fov_pt.norm();
        double theta_ver = -atan2(fov_pt.z(), fov_pt.x());
        double theta_hor = -atan2(fov_pt.y(), fov_pt.x());
        int ver = floor((theta_ver + cfg_.vertical_fov) / cfg_.angular_reso);
        int hor = floor((theta_hor + cfg_.horizon_fov) / cfg_.angular_reso);
        Vector2i idx(ver, hor);
        return idx;
    }

    Vector3d toFovPoint(const Vector2i &idx) {
        double depth = depth_map_[idx.x()][idx.y()];
        Vec3 fov_pt;
        double ver = idx.x() * cfg_.angular_reso - cfg_.vertical_fov;
        double hor = idx.y() * cfg_.angular_reso - cfg_.horizon_fov;

        fov_pt.x() = sqrt(depth*depth / (1 + tan(ver)* tan(ver) + tan(hor)* tan(hor)));
        fov_pt.z() = -tan(ver) * fov_pt.x();
        fov_pt.y() = -tan(hor) * fov_pt.x();
        return fov_pt;
    }

    // 问题化简为搜索圆锥而不是托圆锥，影响不大。
    double GetMaxRadiusInFov(const Vector3d &seed) {
        Vec3 seed_in_fov = toFovFrame(seed);
        double l = seed_in_fov.norm();
        Vec3 v1 = seed_in_fov.normalized();
        Vec3 v2(1, 0, 0);
        double sensin_ang = min(cfg_.vertical_fov, cfg_.horizon_fov);
        double theta_1 = acos(v1.dot(v2));
        double theta_2 = sensin_ang - theta_1;
        double max_r;
        if (theta_2 < 0) {
            max_r = 0;
            return 0;
        }
        max_r = l * sin(theta_2);
        return max_r;
    }

    void LocalCloudCallback(const sensor_msgs::PointCloud2ConstPtr &msg) {
        if (!has_odom_) return;
        pcl::PointCloud<pcl::PointXYZ> cloud_input;
        pcl::fromROSMsg(*msg, (cloud_input));
        // Update camera pose
        UpdateCamPose();
        ClearDepthGrid();

        pcl::PointCloud<pcl::PointXYZ> pc_in_fov;
//        TimeConsuming t__("to fov");
//        t__.start();
//        print("{} points in cur frame.\n", cloud_input.points.size()*1.0/cfg_.point_filt_num);

        Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
        transform.block<3,3>(0,0)= cam_pose_.rotation.inverse();
        transform.col(3).head(3) =- cam_pose_.position;
        Eigen::Affine3d transform_to_fov( transform);
        pcl::transformPointCloud(cloud_input,pc_in_fov,transform);
        // Update depth map
        for (size_t i = 0; i < (pc_in_fov.points.size()); i+= cfg_.point_filt_num) {
            pcl::PointXYZ pt_in = pc_in_fov.points[i];
            Eigen::Vector3d pose_pt(pt_in.x, pt_in.y, pt_in.z);
            if(pose_pt.norm()>cfg_.sensing_range){
                continue;
            }
            double depth;
            Vector2i idx = toGridIndex(pose_pt, depth);
            if (!IsGridIdValid(idx)) continue;
            if ((depth > depth_map_[idx.x()][idx.y()] && depth_map_[idx.x()][idx.y()] != 0)) {
                continue;
            }
            if(depth_map_ [idx.x()][idx.y()] == 0){
                depth_map_ [idx.x()][idx.y()] = cfg_.sensing_range;
            }
            depth_map_[idx.x()][idx.y()] = min(depth,depth_map_[idx.x()][idx.y()]);
//            if (!IsGridIdValid(idx)) continue;
//            for(float dx = -0.1; dx <=0.1; dx+=0.1){
//                for(float dy = -0.1; dy <=0.1; dy+=0.1){
//                    for(float dz = -0.1; dz <=0.1; dz+=0.1)
//                    {
//                        double fake_dep;
//                        Vec3 pt = pose_pt + Vec3(dx,dy,dz);
//                        Vector2i idx = toGridIndex(pt, fake_dep);
//                        if (!IsGridIdValid(idx)) continue;
//                        if ((depth > depth_map_[idx.x()][idx.y()] && depth_map_[idx.x()][idx.y()] != 0)) {
//                            continue;
//                        }
//                        if(depth_map_ [idx.x()][idx.y()] == 0){
//                            depth_map_ [idx.x()][idx.y()] = cfg_.sensing_range;
//                        }
//                        depth_map_[idx.x()][idx.y()] = min(depth,depth_map_[idx.x()][idx.y()]);
//                    }
//                }
//            }


        }
//        t__.stop();
        ConvertDepthGridToImage();
        ConvertDepthGridToPointCloud();
    }

    void ConvertDepthGridToPointCloud(){
        pcl::PointCloud<pcl::PointXYZ> pc;
        for(int ver = 0 ; ver < cfg_.grid_size_ver ; ver ++){
            for(int hor = 0 ; hor < cfg_.grid_size_hor; hor ++){
                Vector2i idx(ver,hor);
                Vec3 pt = toFovPoint(idx);
                if(pt.norm() < 0.1){
                    continue;
                }
                pt = toWorldFrame(pt);
                pcl::PointXYZ pcl_pt(pt.x(), pt.y(),pt.z());
                pc.points.push_back(pcl_pt);
            }
        }
        sensor_msgs::PointCloud2 pc2;
        pcl::toROSMsg(pc,pc2);
        pc2.header.frame_id = "world";
        pc2.header.stamp = ros::Time::now();
        depth_pc_pub_.publish(pc2);

    }

    void ConvertDepthGridToImage() {
        tinycolormap::ColormapType ct = tinycolormap::ColormapType::Jet;
        cv::Mat img(cfg_.grid_size_ver, cfg_.grid_size_hor, CV_8UC3, cv::Vec3b(0, 0, 0));  //创建一副图像
        for (int i = 0; i < cfg_.grid_size_ver; i++) {
            for (int j = 0; j < cfg_.grid_size_hor; j++) {
                double depth = depth_map_[i][j];
                if (depth == 0) {
                    continue;
                }
                double color_id = depth / cfg_.sensing_range;
                if (color_id > 1.0) { color_id = 1.0; }

                Eigen::Vector3d color_mag = tinycolormap::GetColor(color_id, ct).ConvertToEigen();
                img.at<cv::Vec3b>(i, j)[0] = static_cast<uint8_t>((color_mag[0] * 255));
                img.at<cv::Vec3b>(i, j)[1] = static_cast<uint8_t>((color_mag[1] * 255));
                img.at<cv::Vec3b>(i, j)[2] = static_cast<uint8_t>((color_mag[2] * 255));
            }
        }

        sensor_msgs::Image img_ros;
        img_ros.encoding = sensor_msgs::image_encodings::BGR8;
        img_ros.header.stamp = ros::Time::now();
        img_ros.header.frame_id = "world";
        toImageMsg(img_ros, img);
        depth_img_pub_.publish(img_ros);
    }

    void ConvertDepthGridToImage(const Vector2i & center,const int& radius) {
        tinycolormap::ColormapType ct = tinycolormap::ColormapType::Jet;
        cv::Mat img(cfg_.grid_size_ver, cfg_.grid_size_hor, CV_8UC3, cv::Vec3b(0, 0, 0));  //创建一副图像
        for (int i = 0; i < cfg_.grid_size_ver; i++) {
            for (int j = 0; j < cfg_.grid_size_hor; j++) {
                double depth = depth_map_[i][j];
                if (depth == 0) {
                    continue;
                }
                double color_id = depth / cfg_.sensing_range;
                if (color_id > 1.0) { color_id = 1.0; }

                Eigen::Vector3d color_mag = tinycolormap::GetColor(color_id, ct).ConvertToEigen();
                img.at<cv::Vec3b>(i, j)[0] = static_cast<uint8_t>((color_mag[0] * 255));
                img.at<cv::Vec3b>(i, j)[1] = static_cast<uint8_t>((color_mag[1] * 255));
                img.at<cv::Vec3b>(i, j)[2] = static_cast<uint8_t>((color_mag[2] * 255));
            }
        }

        cv::circle(img, cv::Point(center.y(), center.x()), radius,cv::Scalar(0,255,0),1);

        sensor_msgs::Image img_ros;
        img_ros.encoding = sensor_msgs::image_encodings::BGR8;
        img_ros.header.stamp = ros::Time::now();
        img_ros.header.frame_id = "world";
        toImageMsg(img_ros, img);
        depth_img_pub_.publish(img_ros);
    }



    void ClearDepthGrid() {
        for (int i = 0; i < cfg_.grid_size_ver; i++) {
            memset(depth_map_[i], 0, cfg_.grid_size_hor * sizeof(double));
        }
    }

    void UpdateCamPose() {
        Eigen::Quaterniond q;
        q.w() = cam_pose_.odom.pose.pose.orientation.w;
        q.x() = cam_pose_.odom.pose.pose.orientation.x;
        q.y() = cam_pose_.odom.pose.pose.orientation.y;
        q.z() = cam_pose_.odom.pose.pose.orientation.z;

        cam_pose_.rotation = q.toRotationMatrix();
        cam_pose_.position = Vec3(cam_pose_.odom.pose.pose.position.x,
                                  cam_pose_.odom.pose.pose.position.y,
                                  cam_pose_.odom.pose.pose.position.z);
    }

    void OdomCallback(const nav_msgs::OdometryConstPtr &msg) {
        has_odom_ = true;
        cam_pose_.odom = *msg;
    }


};

#endif //_DEPTH_GRID_HPP_
