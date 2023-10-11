#ifndef PLANNER_H
#define PLANNER_H

#include <ros/ros.h>
#include <ros/package.h>

#include <fstream>
#include <mutex>
#include <tf/tf.h>

#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/Imu.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <quadrotor_msgs/PositionCommand.h>
#include <mavros_msgs/AttitudeTarget.h>
#include <mavros_msgs/RCIn.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/crop_box.h>

#include "../include/astar.h"
#include "../include/mpc.h"
#include "../include/local_astar.h"

enum UAVMode_e {
    Manual  = 0,
    Hover   = 1,
    Takeoff = 2,
    Land    = 3,
    Command = 4
};

class PlannerClass {
public:
    // PlannerClass() {}
    PlannerClass(ros::NodeHandle &nh) {
        // load param
        double freq;
        nh.param("/ipc_node/simulation", simu_flag_, true);
        nh.param("/ipc_node/perfect_simu", perfect_simu_flag_, false);
        nh.param("/ipc_node/frequency", freq, 100.0);
        nh.param("/ipc_node/ctrl_delay", ctrl_delay_, 0.1);
        nh.param("/ipc_node/sfc_dis", sfc_dis_, 0.1);
        nh.param("/ipc_node/thrust_limit", thrust_limit_, 0.5);
        nh.param("/ipc_node/hover_esti", hover_esti_flag_, true);
        nh.param("/ipc_node/hover_perc", hover_perc_, 0.23);
        nh.param("/ipc_node/yaw_ctrl_flag", yaw_ctrl_flag_, false);
        nh.param("/ipc_node/yaw_gain", yaw_gain_, 0.1);
        nh.param("/ipc_node/goal_x", goal_p_.x(), 0.0);
        nh.param("/ipc_node/goal_y", goal_p_.y(), 0.0);
        nh.param("/ipc_node/goal_z", goal_p_.z(), 1.0);
        Eigen::Vector3d map_size, map_low, map_upp;
        nh.param("/ipc_node/astar/resolution", resolution_, 0.1);
        nh.param("/ipc_node/astar/map_x_size", map_size.x(), 10.0);
        nh.param("/ipc_node/astar/map_y_size", map_size.y(), 10.0);
        nh.param("/ipc_node/astar/map_z_size", map_size.z(), 5.0);
        nh.param("/ipc_node/astar/expand_dyn", expand_dyn_, 0.25);
        nh.param("/ipc_node/astar/expand_fix", expand_fix_, 0.25);
        map_low << -map_size.x()/2.0, -map_size.y()/2.0, 0;
        map_upp <<  map_size.x()/2.0,  map_size.y()/2.0, map_size.z();
        map_upp_ = map_upp;
        nh.param("/ipc_node/fsm/ref_dis", ref_dis_, 1);
        nh.param("/ipc_node/fsm/path_dis", path_dis_, 0.1);

        // instantiation
        odom_p_ << 0, 0, 1;
        Gravity_ << 0, 0, 9.81;
        odom_v_.setZero();
        odom_a_.setZero();
        imu_a_.setZero();
        if (simu_flag_) {
            thrust_ = 0.7;
            mode_ = Command;
        } else {
            thrust_ = hover_perc_;
            mode_ = Manual;
        }
        thr2acc_ = 9.81 / thrust_;

        mpc_   = std::make_shared<MPCPlannerClass>(nh);
        local_astar_ = std::make_shared<LoaclAstarClass>();
        local_astar_->InitMap(resolution_, map_low, map_upp);
        // local_astar_->InitMap(resolution_, Eigen::Vector3d(-5, -5, 0), Eigen::Vector3d(5, 5, 3));

        // ros topic
        astar_pub_ = nh.advertise<visualization_msgs::Marker>("astar_path", 1);
        gird_map_pub_ = nh.advertise<sensor_msgs::PointCloud2>("grid_map", 1);
        cmd_pub_ = nh.advertise<quadrotor_msgs::PositionCommand>("cmd", 1);
        mpc_path_pub_ = nh.advertise<nav_msgs::Path>("mpc_path", 1);
        sfc_pub_ = nh.advertise<visualization_msgs::MarkerArray>("sfc", 1);
        px4ctrl_pub_ = nh.advertise<mavros_msgs::AttitudeTarget>("px4ctrl", 1);
        goal_pub_ = nh.advertise<geometry_msgs::PoseStamped>("goal_pub", 1);

        local_pc_sub_ = nh.subscribe<sensor_msgs::PointCloud2>("local_pc", 1, &PlannerClass::LocalPcCallback, this, ros::TransportHints().tcpNoDelay());
        odom_sub_ = nh.subscribe<nav_msgs::Odometry>("odom", 1, &PlannerClass::OdomCallback, this, ros::TransportHints().tcpNoDelay());
        goal_sub_ = nh.subscribe<geometry_msgs::PoseStamped>("goal", 1, &PlannerClass::GoalCallback, this, ros::TransportHints().tcpNoDelay());
        imu_sub_ = nh.subscribe<sensor_msgs::Imu>("imu", 1, &PlannerClass::IMUCallback, this, ros::TransportHints().tcpNoDelay());
        rc_sub_ = nh.subscribe<mavros_msgs::RCIn>("rc_in", 1, &PlannerClass::RCInCallback, this, ros::TransportHints().tcpNoDelay());

        timer_ = nh.createTimer(ros::Duration(1.0/freq), &PlannerClass::TimerCallback, this, false, true);

        std::string file = ros::package::getPath("ipc") + "/config";
        write_time_.open((file+"/time_consuming.csv"), std::ios::out | std::ios::trunc);
        log_times_.resize(5, 1);
        write_time_ << "mapping" << ", " << "replan" << ", " << "sfc" << ", " << "mpc" << ", " << "df" << ", " << std::endl;
    }
    ~PlannerClass() {}

private:
    void AstarPublish(std::vector<Eigen::Vector3d>& nodes, uint8_t type, double scale) {
        visualization_msgs::Marker node_vis; 
        node_vis.header.frame_id = "world";
        node_vis.header.stamp = ros::Time::now();

        if (type == 0) {
            node_vis.ns = "astar_path";
            node_vis.color.a = 1.0;
            node_vis.color.r = 0.0;
            node_vis.color.g = 0.0;
            node_vis.color.b = 0.0;
        } else if (type == 1) {
            node_vis.ns = "floyd_path";
            node_vis.color.a = 1.0;
            node_vis.color.r = 1.0;
            node_vis.color.g = 0.0;
            node_vis.color.b = 0.0;
        } else if (type == 2) {
            node_vis.ns = "short_path";
            node_vis.color.a = 1.0;
            node_vis.color.r = 0.0;
            node_vis.color.g = 0.0;
            node_vis.color.b = 1.0;
        } else if (type == 3) {
            node_vis.ns = "set_points";
            node_vis.color.a = 1.0;
            node_vis.color.r = 0.0;
            node_vis.color.g = 1.0;
            node_vis.color.b = 0.0;
        }

        node_vis.type = visualization_msgs::Marker::CUBE_LIST;
        node_vis.action = visualization_msgs::Marker::ADD;
        node_vis.id = 0;
        node_vis.pose.orientation.x = 0.0;
        node_vis.pose.orientation.y = 0.0;
        node_vis.pose.orientation.z = 0.0;
        node_vis.pose.orientation.w = 1.0;
        
        node_vis.scale.x = scale;
        node_vis.scale.y = scale;
        node_vis.scale.z = scale;

        geometry_msgs::Point pt;
        for (int i = 0; i < int(nodes.size()); i++) {
            Eigen::Vector3d coord = nodes[i];
            pt.x = coord(0);
            pt.y = coord(1);
            pt.z = coord(2);
            node_vis.points.push_back(pt);
        }
        astar_pub_.publish(node_vis);
    }
    void CmdPublish(Eigen::Vector3d p_r, Eigen::Vector3d v_r, Eigen::Vector3d a_r, Eigen::Vector3d j_r) {
        quadrotor_msgs::PositionCommand msg;
        msg.header.frame_id = "world";
        msg.header.stamp    = ros::Time::now();
        msg.position.x      = p_r.x();
        msg.position.y      = p_r.y();
        msg.position.z      = p_r.z();
        msg.velocity.x      = v_r.x();
        msg.velocity.y      = v_r.y();
        msg.velocity.z      = v_r.z();
        msg.acceleration.x  = a_r.x();
        msg.acceleration.y  = a_r.y();
        msg.acceleration.z  = a_r.z();
        msg.jerk.x          = j_r.x();
        msg.jerk.y          = j_r.y();
        msg.jerk.z          = j_r.z();
        if (yaw_ctrl_flag_) {
            double yaw_error = yaw_r_ - yaw_;
            if (yaw_error >  M_PI) yaw_error -= M_PI * 2;
            if (yaw_error < -M_PI) yaw_error += M_PI * 2;
            msg.yaw     = yaw_ + yaw_error * 0.1;
            msg.yaw_dot = 0;
        } else {
            msg.yaw     = 0;
            msg.yaw_dot = 0;
        }
        cmd_pub_.publish(msg);
    }
    void MPCPathPublish(std::vector<Eigen::Vector3d> &pt) {
        nav_msgs::Path msg;
        msg.header.frame_id = "world";
        msg.header.stamp = ros::Time::now();
        for (int i = 0; i < pt.size(); i++) {
            geometry_msgs::PoseStamped pose;
            pose.pose.position.x = pt[i].x();
            pose.pose.position.y = pt[i].y();
            pose.pose.position.z = pt[i].z();
            msg.poses.push_back(pose);
            // std::cout << "mpc path " << i << ": " << pt[i].transpose() << std::endl;
        }
        mpc_path_pub_.publish(msg);
    }
    void WriteLogTime(void) {
        for (int i = 0; i < log_times_.size(); i++) {
            write_time_ << log_times_[i] << ", ";
            log_times_[i] = 0.0;
        }
        write_time_ << std::endl;
    }
    
    void ComputeThrust(Eigen::Vector3d acc) {
        const Eigen::Vector3d zB = odom_q_ * Eigen::Vector3d::UnitZ();
        double des_acc_norm = acc.dot(zB);
        thrust_ = des_acc_norm / thr2acc_;
    }
    void ConvertCommand(Eigen::Vector3d acc, Eigen::Vector3d jerk) {
        Eigen::Vector3d xB, yB, zB, xC;
        if (yaw_ctrl_flag_) {
            double yaw_error = yaw_r_ - yaw_;
            if (yaw_error >  M_PI) yaw_error -= M_PI * 2;
            if (yaw_error < -M_PI) yaw_error += M_PI * 2;
            yaw_dot_r_ = yaw_error * yaw_gain_;
        } else {
            yaw_dot_r_ = (0 - yaw_) * yaw_gain_;
        }
        xC << std::cos(yaw_), std::sin(yaw_), 0;

        zB = acc.normalized();
        yB = (zB.cross(xC)).normalized();
        xB = yB.cross(zB);
        Eigen::Matrix3d R;
        R << xB, yB, zB;
        u_q_ = R;

        Eigen::Vector3d hw = (jerk - (zB.dot(jerk) * zB)) / acc.norm();
        rate_.x() = -hw.dot(yB);
        rate_.y() = hw.dot(xB);
        rate_.z() = yaw_dot_r_ * zB.dot(Eigen::Vector3d(0, 0, 1));
    }
    bool estimateThrustModel(const Eigen::Vector3d &est_a)
    {
        if (hover_esti_flag_ == false) {
            thr2acc_ = 9.81 / hover_perc_;
            return true;
        }
        if (mode_ != Command) {
            P_ = 100.0;
            thr2acc_ = 9.81 / hover_perc_;
            return true;
        }
        ros::Time t_now = ros::Time::now();
        if (timed_thrust_.size() == 0) return false;
        std::pair<ros::Time, double> t_t = timed_thrust_.front();

        while (timed_thrust_.size() >= 1) {
            double delta_t = (t_now - t_t.first).toSec();
            if (delta_t > 1.0) {
                timed_thrust_.pop();
                continue;
            } 
            if (delta_t < 0.035) {
                return false;
            }

            /* Recursive least squares algorithm with vanishing memory */
            double thr = t_t.second;
            timed_thrust_.pop();
            /* Model: est_a(2) = thr2acc * thr */
            double R = 0.3; // using Kalman filter
            double K = P_ / (P_ + R);
            thr2acc_ = thr2acc_ + K * (est_a(2) - thr * thr2acc_);
            P_ = (1 - K * thr) * P_;
            double hover_percentage = 9.81 / thr2acc_;
            if (hover_percentage > 0.8 || hover_percentage < 0.1) {
                // ROS_INFO_THROTTLE(1, "Estimated hover_percentage >0.8 or <0.1! Perhaps the accel vibration is too high!");
                thr2acc_ = hover_percentage > 0.8 ? 9.81 / 0.8 : thr2acc_;
                thr2acc_ = hover_percentage < 0.1 ? 9.81 / 0.1 : thr2acc_;
            }
            // ROS_WARN("[PX4CTRL] hover_percentage = %f", hover_percentage);
            return true;
        }
        return false;
    }
    void AttitudeCtrlPub(const Eigen::Quaterniond &q, const double thrust, const ros::Time &stamp) {
        mavros_msgs::AttitudeTarget msg;
        msg.header.stamp = stamp;
        msg.header.frame_id = std::string("FCU");
        msg.type_mask = mavros_msgs::AttitudeTarget::IGNORE_ROLL_RATE |
                        mavros_msgs::AttitudeTarget::IGNORE_PITCH_RATE |
                        mavros_msgs::AttitudeTarget::IGNORE_YAW_RATE;
        msg.orientation.x = q.x();
        msg.orientation.y = q.y();
        msg.orientation.z = q.z();
        msg.orientation.w = q.w();
        msg.thrust = thrust;
        if (msg.thrust < 0.1) msg.thrust = 0.1;
        if (msg.thrust > 0.9) msg.thrust = 0.9;
        if (mode_ == Manual) msg.thrust = 0.05;
        if (!simu_flag_ && msg.thrust > thrust_limit_) msg.thrust = thrust_limit_;
        px4ctrl_pub_.publish(msg);
    }
    void BodyrateCtrlPub(const Eigen::Vector3d &rate, const double thrust, const ros::Time &stamp) {
        mavros_msgs::AttitudeTarget msg;
        msg.header.stamp = stamp;
        msg.header.frame_id = std::string("FCU");
        msg.type_mask = mavros_msgs::AttitudeTarget::IGNORE_ATTITUDE;
        msg.body_rate.x = rate.x();
        msg.body_rate.y = rate.y();
        msg.body_rate.z = rate.z();
        msg.thrust = thrust;
        if (msg.thrust < 0.08) msg.thrust = 0.08;
        if (msg.thrust > 0.9) msg.thrust = 0.9;
        if (mode_ == Manual) msg.thrust = 0.05;
        if (!simu_flag_ && msg.thrust > thrust_limit_) msg.thrust = thrust_limit_;
        px4ctrl_pub_.publish(msg);
    }

    void OdomCallback(const nav_msgs::OdometryConstPtr& msg);
    void GoalCallback(const geometry_msgs::PoseStampedConstPtr& msg);
    void IMUCallback(const sensor_msgs::ImuConstPtr& msg);
    void LocalPcCallback(const sensor_msgs::PointCloud2ConstPtr& msg);
    void RCInCallback(const mavros_msgs::RCInConstPtr& msg);

    void TimerCallback(const ros::TimerEvent &);

    void PathReplan(bool extend);
    void GeneratePolyOnPath();
    void GenerateAPolytope(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Matrix<double, Eigen::Dynamic, 4>& planes, uint8_t index);

    ros::Timer      timer_;
    ros::Subscriber odom_sub_, goal_sub_, imu_sub_, local_pc_sub_, rc_sub_;
    ros::Publisher  gird_map_pub_, astar_pub_, cmd_pub_, sfc_pub_, mpc_path_pub_, px4ctrl_pub_, goal_pub_;
    ros::Time       odom_time_, last_mpc_time_;
    std::mutex  odom_mutex_, goal_mutex_, cloud_mutex_, timer_mutex_, imu_mutex_, 
                local_pc_mutex_, rc_mutex_;

    UAVMode_e mode_;

    bool simu_flag_, perfect_simu_flag_, pc_ctrl_flag_, hover_esti_flag_, yaw_ctrl_flag_;
    bool has_map_flag_{false}, has_odom_flag_{false}, replan_flag_{false}, new_goal_flag_{false};
    double resolution_;
    double ctrl_delay_;
    double thrust_limit_, hover_perc_;
    double sfc_dis_, path_dis_, expand_dyn_, expand_fix_;
    int ref_dis_;
    int mpc_ctrl_index_;
    std::vector<Eigen::Vector3d> remain_nodes_;

    std::ofstream write_time_;
    std::vector<double> log_times_;

    Eigen::Vector3d goal_p_, map_upp_;
    Eigen::Vector3d odom_p_, odom_v_, odom_a_, imu_a_;
    Eigen::Quaterniond odom_q_, u_q_;
    Eigen::Vector3d rate_;
    double yaw_{0}, yaw_r_{0}, yaw_dot_r_{0}, yaw_gain_;

    int astar_index_{0};
    std::vector<Eigen::Vector3d> astar_path_;
    std::vector<Eigen::Vector3d> waypoints_;
    std::vector<Eigen::Vector3d> follow_path_;
    std::vector<Eigen::Vector3d> replan_path_;
    std::vector<Eigen::Vector3d> local_pc_;
    std::vector<Eigen::Vector3d> local_pc_buffer_[10];
    std::vector<Eigen::Vector3d> mpc_goals_;

    double thr2acc_;
    double thrust_;
    double P_{100.0};
    Eigen::Vector3d Gravity_;
    std::queue<std::pair<ros::Time, double>> timed_thrust_;

    std::deque<pcl::PointCloud<pcl::PointXYZ>> vec_cloud_;
    pcl::PointCloud<pcl::PointXYZ> static_map_;
    pcl::PointCloud<pcl::PointXYZ>::Ptr static_cloud_;

    std::shared_ptr<LoaclAstarClass> local_astar_;
    std::shared_ptr<MPCPlannerClass> mpc_;
};

#endif
