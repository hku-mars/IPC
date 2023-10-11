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

#define TINYCOLORMAP_WITH_EIGEN

#include "tinycolormap.hpp"

using namespace std;


ros::Publisher color_pc_pub;
typedef pcl::PointXYZRGB PointType;
pcl::PointCloud<PointType> color_pc;
struct ColorMapParams {
    double max_height;
    double min_height;
    double range;
    string map_type_name, frame_id;
    bool init_map_success{false};
    bool static_pc_mode,inverse_color{false};
    double publish_rate;
} cmp_;

void PubColorCloudCallback(const ros::TimerEvent &e) {
    if (!cmp_.init_map_success) {
        return;
    }
    color_pc_pub.publish(color_pc);
}

void DoubleCalp(double &data_in, const double min_value, const double max_value) {
    data_in = data_in > max_value ? max_value : data_in;
    data_in = data_in < min_value ? min_value : data_in;
}

void PointCloudCallback(const sensor_msgs::PointCloud2ConstPtr &msg) {
    if (cmp_.init_map_success) {
        return;
    }
    tinycolormap::ColormapType ct;
    if (cmp_.map_type_name == "Parula") {
        ct = tinycolormap::ColormapType::Parula;
    }
    if (cmp_.map_type_name == "Heat") {
        ct = tinycolormap::ColormapType::Heat;
    }
    if (cmp_.map_type_name == "Jet") {
        ct = tinycolormap::ColormapType::Jet;
    }
    if (cmp_.map_type_name == "Turbo") {
        ct = tinycolormap::ColormapType::Turbo;
    }
    if (cmp_.map_type_name == "Hot") {
        ct = tinycolormap::ColormapType::Hot;
    }
    if (cmp_.map_type_name == "Gray") {
        ct = tinycolormap::ColormapType::Gray;
    }
    if (cmp_.map_type_name == "Magma") {
        ct = tinycolormap::ColormapType::Magma;
    }
    if (cmp_.map_type_name == "Inferno") {
        ct = tinycolormap::ColormapType::Inferno;
    }
    if (cmp_.map_type_name == "Plasma") {
        ct = tinycolormap::ColormapType::Plasma;
    }
    if (cmp_.map_type_name == "Viridis") {
        ct = tinycolormap::ColormapType::Viridis;
    }
    if (cmp_.map_type_name == "Cividis") {
        ct = tinycolormap::ColormapType::Cividis;
    }
    if (cmp_.map_type_name == "Github") {
        ct = tinycolormap::ColormapType::Github;
    }
    pcl::PCLPointCloud2::Ptr cloud(new pcl::PCLPointCloud2);
    pcl_conversions::toPCL(*msg, *cloud);
    pcl::fromPCLPointCloud2( *cloud, color_pc);

    for (size_t i = 0; i < color_pc.size(); i++) {
        double color_id = color_pc.points[i].z;
        DoubleCalp(color_id, cmp_.min_height, cmp_.max_height);
        color_id -= cmp_.min_height;
        color_id /= cmp_.range;
        if(cmp_.inverse_color){
            color_id = 1.0 - color_id;
        }
        Eigen::Vector3d color_mag = tinycolormap::GetColor(color_id, ct).ConvertToEigen();
        color_pc.points[i].r = static_cast<uint8_t>((color_mag[0] * 255));
        color_pc.points[i].g = static_cast<uint8_t>((color_mag[1] * 255));;
        color_pc.points[i].b = static_cast<uint8_t>((color_mag[2] * 255));;
    }

    std_msgs::Header hd;
    hd.frame_id = cmp_.frame_id;
    pcl_conversions::toPCL(hd, color_pc.header);

    cmp_.init_map_success = true;
    printf("111111111111111111111.\n");
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "color_map_pc2");
    ros::NodeHandle n("~");

    ros::Subscriber pc_sub = n.subscribe("cloud", 1000000, PointCloudCallback);

    color_pc_pub = n.advertise<pcl::PointCloud<pcl::PointXYZRGB>>("color_cloud", 100000);
    n.param("color/min_height", cmp_.min_height, 0.0);
    n.param("color/max_height", cmp_.max_height, 5.0);
    n.param("color/map_type_name", cmp_.map_type_name, std::string("Turbo"));
    n.param("color/frame_id", cmp_.frame_id, std::string("world"));
    n.param("color/publish_rate", cmp_.publish_rate, 1.0);
    n.param("color/inverse_color", cmp_.inverse_color, false);
    cmp_.range = cmp_.max_height - cmp_.min_height;
    double dt = 1.0 / cmp_.publish_rate;
    ros::Timer pub_pc_timer = n.createTimer(ros::Duration(dt), &PubColorCloudCallback);

    ros::AsyncSpinner spinner(0);
    spinner.start();

    ros::waitForShutdown();
}
