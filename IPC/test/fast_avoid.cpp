#include <ros/ros.h>
#include <ros/package.h>

#include <Eigen/Eigen>

#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <sensor_msgs/PointCloud2.h>
#include <quadrotor_msgs/PositionCommand.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <fstream>

ros::Publisher goal_pub, local_pc_pub, map_pub, obs_pub;
Eigen::Vector3d odom_p, odom_v, odom_p_last;
Eigen::Vector3d obs_pos, obs_pos_r, obs_size;
double obs_dis = 1.0, min_dis = 10000.0, obs_v, resolusion, odom_p_dt, v_f;
bool obs_flag = false, first_flag = false;
pcl::PointCloud<pcl::PointXYZ> cloud, map_cloud;
ros::Time odom_time;

void OdomCallback(const nav_msgs::OdometryConstPtr& msg)
{
    odom_p_last = odom_p;
    odom_p << msg->pose.pose.position.x, msg->pose.pose.position.y, msg->pose.pose.position.z;
    odom_v << msg->twist.twist.linear.x, msg->twist.twist.linear.y, msg->twist.twist.linear.z;
    for (int i = 0; i < cloud.points.size(); i++) {
        Eigen::Vector3d pt;
        pt << cloud.points[i].x, cloud.points[i].y, cloud.points[i].z;
        // double dis = (odom_p - pt).norm();
        double dis = std::sqrt((odom_p.x()-pt.x())*(odom_p.x()-pt.x()) + (odom_p.y()-pt.y())*(odom_p.y()-pt.y()));
        if (dis < min_dis) {
            min_dis = dis;
        }
    }

    static bool print_flag = false;
    if (std::fabs(odom_p.x() - 10.0) < 0.5 && print_flag == false) {
        print_flag = true;
        min_dis -= 0.25;
        if (min_dis > 0) ROS_INFO("\033[1;32m minimum distance: %f \033[0m", min_dis);
        else ROS_ERROR("minimum distance: %f", min_dis);
    }

    odom_p_dt = (msg->header.stamp - odom_time).toSec();
    odom_time = msg->header.stamp;
}

void GoalPathPublish(void)
{
    // nav_msgs::Path msg;
    // msg.header.frame_id = "world";
    // msg.header.stamp = ros::Time::now();
    // geometry_msgs::PoseStamped point;
    // point.pose.position.x = 0.0;
    // point.pose.position.y = 0.0;
    // point.pose.position.z = 1.0;
    // msg.poses.push_back(point);
    // point.pose.position.x = 10.0;
    // point.pose.position.y = 0.0;
    // point.pose.position.z = 1.0;
    // msg.poses.push_back(point);
    // goal_pub.publish(msg);
    
    geometry_msgs::PoseStamped msg;
    msg.header.frame_id = "world";
    msg.header.stamp = ros::Time::now();
    msg.pose.position.x = 10.0;
    msg.pose.position.y = 0.0;
    msg.pose.position.z = 1.0;
    goal_pub.publish(msg);
}

void LocalPCPublish(bool flag)
{
    pcl::PointCloud<pcl::PointXYZ> pub_cloud;
    pub_cloud += map_cloud;

    if (flag) {
        cloud.clear();
        obs_pos.y() -= obs_v * 0.1;
        if (obs_pos.y() < obs_pos_r.y()) obs_pos.y() = obs_pos_r.y();
        for (double i = obs_pos.x()-obs_size.x()/2; i <= obs_pos.x()+obs_size.x()/2; i += resolusion) {
            for (double j = obs_pos.y()-obs_size.y()/2; j <= obs_pos.y()+obs_size.y()/2; j += resolusion) {
                for (double k = obs_pos.z()-obs_size.z()/2; k <= obs_pos.z()+obs_size.z()/2; k += resolusion) {
                    cloud.push_back(pcl::PointXYZ(i, j, k));
                }
            }
        }
    }

    pub_cloud += cloud;
    pub_cloud.width    = pub_cloud.points.size();
    pub_cloud.height   = 1;
    pub_cloud.is_dense = true;

    sensor_msgs::PointCloud2 msg;
    pcl::toROSMsg(pub_cloud, msg);
    msg.header.frame_id = "world";
    msg.header.stamp = ros::Time::now();
    local_pc_pub.publish(msg);

    sensor_msgs::PointCloud2 msg2;
    pcl::toROSMsg(cloud, msg2);
    msg2.header.frame_id = "world";
    msg2.header.stamp = ros::Time::now();
    if (flag) obs_pub.publish(msg2);
}

void PVAJCallbacl(const quadrotor_msgs::PositionCommandConstPtr& msg)
{
    if (std::fabs(msg->acceleration.x) > 20) min_dis = -0.25;
    if (std::fabs(msg->acceleration.y) > 20) min_dis = -0.25;
    if (msg->acceleration.z > 20 || msg->acceleration.z < -10) min_dis = -0.25;
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "fast_avoid_node");
    ros::NodeHandle nh;

    nh.param("/fast_avoid/resolusion", resolusion, 0.1);
    nh.param("/fast_avoid/obs_dis", obs_dis, 1.0);
    nh.param("/fast_avoid/obs_v", obs_v, 1.0);
    nh.param("/fast_avoid/obs_x", obs_pos_r.x(), 0.0);
    nh.param("/fast_avoid/obs_y", obs_pos_r.y(), 0.0);
    nh.param("/fast_avoid/obs_z", obs_pos_r.z(), 0.0);
    nh.param("/fast_avoid/obs_l", obs_size.x(), 0.0);
    nh.param("/fast_avoid/obs_w", obs_size.y(), 0.0);
    nh.param("/fast_avoid/obs_h", obs_size.z(), 0.0);

    obs_pos << obs_pos_r.x(), obs_pos_r.y() + 1.5, obs_pos_r.z();
    std::cout << "obs_dis:   " << obs_dis << std::endl;
    std::cout << "obs_pos:   " << obs_pos.transpose() << std::endl;
    std::cout << "obs_pos_r: " << obs_pos_r.transpose() << std::endl;
    std::cout << "obs_size:  " << obs_size.transpose() << std::endl;

    ros::Subscriber odom_sub = nh.subscribe<nav_msgs::Odometry>("odom", 1, &OdomCallback, ros::TransportHints().tcpNoDelay());
    ros::Subscriber pvaj_sub = nh.subscribe<quadrotor_msgs::PositionCommand>("pvaj", 1, &PVAJCallbacl, ros::TransportHints().tcpNoDelay());
    
    goal_pub = nh.advertise<geometry_msgs::PoseStamped>("goal", 1);
    // goal_pub = nh.advertise<nav_msgs::Path>("goal", 1);
    local_pc_pub = nh.advertise<sensor_msgs::PointCloud2>("local_pc", 1);
    map_pub = nh.advertise<sensor_msgs::PointCloud2>("map", 1);
    obs_pub = nh.advertise<sensor_msgs::PointCloud2>("obs", 1);

    ros::Rate(1).sleep();

    pcl::PointCloud<pcl::PointXYZ> map_cloud_1, obs_cloud;
    for (double x = 0.0; x < 10.0; x += resolusion) {
        if (std::fabs(x - obs_pos_r.x()) < obs_size.x() + 0.0) continue;
        for (double y = 0.8; y < 1.0; y += resolusion) {
            for (double z = 0.0; z < 3.0; z += resolusion) {
                map_cloud_1.push_back(pcl::PointXYZ(x, y, z));
            }
        }
    }
    map_cloud_1.width    = map_cloud_1.points.size();
    map_cloud_1.height   = 1;
    map_cloud_1.is_dense = true;
    sensor_msgs::PointCloud2 msg;
    pcl::toROSMsg(map_cloud_1, msg);
    msg.header.frame_id = "world";
    msg.header.stamp = ros::Time::now();
    map_pub.publish(msg);

    for (double i = obs_pos.x()-obs_size.x()/2; i <= obs_pos.x()+obs_size.x()/2; i += resolusion) {
        for (double j = obs_pos.y()-obs_size.y()/2; j <= obs_pos.y()+obs_size.y()/2; j += resolusion) {
            for (double k = obs_pos.z()-obs_size.z()/2; k <= obs_pos.z()+obs_size.z()/2; k += resolusion) {
                obs_cloud.push_back(pcl::PointXYZ(i, j, k));
            }
        }
    }
    sensor_msgs::PointCloud2 msg2;
    pcl::toROSMsg(obs_cloud, msg2);
    msg2.header.frame_id = "world";
    msg2.header.stamp = ros::Time::now();
    obs_pub.publish(msg2);

    for (double x = 0.0; x < 10.0; x += resolusion) {
        // if (std::fabs(x - obs_pos_r.x()) < obs_size.x() + 0.0) continue;
        for (double y = 0.8; y < 1.0; y += resolusion) {
            for (double z = 0.0; z < 3.0; z += resolusion) {
                map_cloud.push_back(pcl::PointXYZ(x, y, z));
            }
        }
    }
    map_cloud.width    = map_cloud.points.size();
    map_cloud.height   = 1;
    map_cloud.is_dense = true;
    // std::string file = ros::package::getPath("ipc") + "/config/map.pcd";
    // pcl::io::savePCDFileASCII(file, map_cloud);

    ros::Rate(1.0).sleep();
    odom_p << 0, 0, 1;
    odom_v = Eigen::Vector3d::Zero();
    GoalPathPublish();

    ros::Rate rate(1000);
    while (ros::ok()) {
        ros::spinOnce();
        static int count = 0;
        count++;
        double dis = std::fabs(odom_p.x() - (obs_pos_r.x()-obs_size.x()/2));
        if (count > 100) { // 10Hz
            count = 0;
            // GoalPathPublish();
            static bool new_flag = false;
            if (dis <= obs_dis && new_flag == false) {
                new_flag = true;
                obs_flag = true;
                obs_pos_r.y() = odom_p.y();
                obs_pos.y() = odom_p.y() + 1.0;
            }
            LocalPCPublish(obs_flag);
        }
        if (first_flag == false && dis <= obs_dis) {
            LocalPCPublish(true);
            first_flag = true;
            rate.sleep();
            // ROS_WARN("obstacle distance: %f, velocity now: %f", obs_dis, odom_v.norm());
            v_f = (odom_p-odom_p_last).norm()/odom_p_dt;
            ROS_WARN("obstacle distance: %f, velocity now: %f", obs_dis, v_f);
        }

        rate.sleep();
    }

    std::ofstream outfile;
    std::string file = ros::package::getPath("mpc_planner") + "/config";
    outfile.open((file+"/mpc.csv"), std::ios::out | std::ios::app);
    if (std::fabs(odom_p.x() - 10.0) > 0.5) min_dis = -0.25;
    outfile << obs_dis << ", " << v_f << ", " << min_dis << ", ";
    outfile << std::endl;

    return 0;
}