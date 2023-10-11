#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>
#include <pcl_conversions/pcl_conversions.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <pcl/search/impl/kdtree.hpp>
#include <pcl/features/normal_3d.h>
#include <pcl/common/common.h>
#include <vector>
#include <cmath>
#include <tr1/unordered_map>
#include "FOV_Checker/FOV_Checker.h"
#include <ros/package.h>
#include <stdlib.h>
#include <time.h>
#include <tf2_ros/transform_broadcaster.h>
#include <geometry_msgs/TransformStamped.h>
#include <pcl/common/transforms.h>

// #include <decay_map/raycast.h>
// #include "hash.h"
// #include <hash_map>
using namespace std;
using namespace Eigen;

#define MAX_INTENSITY 200
#define MIN_INTENSITY 1

std::string pkg_path;
std::ofstream myfile;

std::string quad_name;
// int uav_num;

// hash_map<int, std::vector<pcl::PointXYZI>> point_hashmap;
std::tr1::unordered_map<int, std::vector<pcl::PointXYZI>> point_hashmap;
std::tr1::unordered_map<int, std::vector<int>> pointindex_hashmap;

pcl::PointCloud<pcl::PointXYZI> cloud_explored;

FOV_Checker fov_checker;
BoxPointType env_box;

PointType min_center, max_center;

int cube_numx, cube_numy, cube_numz;

struct polar3D {
    int theta;
    int fi;
    float r;
};

ros::Publisher pub_cloud, pub_pose, pub_intercloud, pub_dyncloud, pub_uavcloud;

sensor_msgs::PointCloud2 local_map_pcl;
sensor_msgs::PointCloud2 local_depth_pcl;

ros::Subscriber odom_sub;
ros::Subscriber global_map_sub, local_map_sub;
ros::Subscriber quad1_odom_sub, quad2_odom_sub;

ros::Timer local_sensing_timer, pose_timer, dynobj_timer;

bool has_global_map(false);
bool has_local_map(false);
bool has_odom(false);
bool has_dyn_map(false);

nav_msgs::Odometry odom_, quad1_odom, quad2_odom;
Eigen::Matrix4d sensor2body, sensor2world;

int output_pcd;
int collisioncheck_enable;
int is_360lidar;
int use_avia_pattern, use_vlp32_pattern, use_minicf_pattern;
int livox_linestep;
double sensing_horizon, sensing_rate, estimation_rate, polar_resolution, yaw_fov, vertical_fov, min_raylength, downsample_res, curvature_limit, hash_cubesize, collision_range;
double x_size, y_size, z_size;
double gl_xl, gl_yl, gl_zl;
double resolution, inv_resolution;
int GLX_SIZE, GLY_SIZE, GLZ_SIZE;

//multi uav variables
int drone_id = 0;
vector<Eigen::Vector3d> other_uav_pos;
vector<vector<pcl::PointXYZI>> otheruav_points, otheruav_points_inrender;
vector<vector<int>> otheruav_pointsindex, otheruav_pointsindex_inrender;
pcl::PointCloud<pcl::PointXYZI> otheruav_points_vis;
double uav_size[3];
int uav_points_num;
int drone_num = 0;

//dynamic objects variables
int dynobj_enable;
double dynobject_size;
vector<Eigen::Vector3d> dynobj_poss;
vector<Eigen::Vector3d> dynobj_dir;
vector<pcl::PointXYZI> dynobj_points;
pcl::PointCloud<pcl::PointXYZI> dynobj_points_vis;
vector<int> dynobj_pointsindex;
double dyn_velocity;
int dynobject_num;
int dyn_mode;
Eigen::Vector3d map_min, map_max;

int origin_mapptcount = 0;

pcl::PointCloud<pcl::Normal>::Ptr all_normals(new pcl::PointCloud<pcl::Normal>);
pcl::NormalEstimation<pcl::PointXYZI, pcl::Normal> normalEstimation;

ros::Time last_odom_stamp = ros::TIME_MAX;

ros::Time start_time, dyn_start_time;

ros::Timer drone_cloud_pub;

pcl::PointCloud<pcl::PointXYZI> cloud_all_map, local_map, local_map_filled, point_in_sensor;
pcl::VoxelGrid<pcl::PointXYZI> _voxel_sampler;
sensor_msgs::PointCloud2 local_map_pcd, sensor_map_pcd;

pcl::search::KdTree<pcl::PointXYZI> _kdtreeLocalMap, kdtree_dyn;
vector<int> pointIdxRadiusSearch;
vector<float> pointRadiusSquaredDistance;

inline Eigen::Vector3d gridIndex2coord(const Eigen::Vector3i &index) {
    Eigen::Vector3d pt;
    pt(0) = ((double) index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double) index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double) index(2) + 0.5) * resolution + gl_zl;

    return pt;
};

inline Eigen::Vector3i coord2gridIndex(const Eigen::Vector3d &pt) {
    Eigen::Vector3i idx;
    idx(0) = std::min(std::max(int((pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1);
    idx(1) = std::min(std::max(int((pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1);
    idx(2) = std::min(std::max(int((pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);

    return idx;
};

sensor_msgs::PointCloud2 dynobj_points_pcd;

void dynobjGenerate(const ros::TimerEvent &event) {
    if (has_global_map == true && dynobj_enable == 1) {
        dynobj_points.clear();
        dynobj_points_vis.clear();
        dynobj_pointsindex.clear();
        pcl::PointXYZI temp_point;
        int dynpt_count = 0;

        //! you can add mode here
        switch (dyn_mode) {
            case 0:
                double fly_time = (ros::Time::now() - dyn_start_time).toSec();
                for (int n = 0; n < dynobject_num; n++) {
                    dynobj_poss[n] = dynobj_poss[n] + fly_time * dynobj_dir[n];
                    if (dynobj_poss[n](0) < map_min(0) || dynobj_poss[n](0) > map_max(0) \
 || dynobj_poss[n](1) < map_min(1) || dynobj_poss[n](1) > map_max(1) \
 || dynobj_poss[n](2) < map_min(2) || dynobj_poss[n](2) > map_max(2)) {
                        dynobj_poss[n](0) = rand() / double(RAND_MAX) * (map_max(0) - map_min(0)) + map_min(0);
                        dynobj_poss[n](1) = rand() / double(RAND_MAX) * (map_max(1) - map_min(1)) + map_min(1);
                        dynobj_poss[n](2) = rand() / double(RAND_MAX) * (map_max(2) - map_min(2)) + map_min(2);
                        Eigen::Vector3d dyntemp_dir_polar;
                        dyntemp_dir_polar(0) = rand() / double(RAND_MAX) * 3.1415926;
                        dyntemp_dir_polar(1) = (rand() / double(RAND_MAX) - 0.5) * 3.1415926;
                        dynobj_dir[n](0) = dyn_velocity * sin(dyntemp_dir_polar(1));
                        dynobj_dir[n](1) = dyn_velocity * cos(dyntemp_dir_polar(1)) * sin(dyntemp_dir_polar(0));
                        dynobj_dir[n](2) = dyn_velocity * cos(dyntemp_dir_polar(1)) * cos(dyntemp_dir_polar(0));
                    }
                }
                break;

        }

        // #pragma omp parallel default (none) \
    //                  shared (dynobject_num, dynobject_size, dynobj_poss, temp_point, dynobj_points_vis,dynobj_points,downsample_res)
        // {
        // #pragma omp for
        for (int n = 0; n < dynobject_num; n++) {
            for (double i = dynobj_poss[n](0) - 0.5 * dynobject_size;
                 i < dynobj_poss[n](0) + 0.5 * dynobject_size; i = i + downsample_res) {
                for (double j = dynobj_poss[n](1) - 0.5 * dynobject_size;
                     j < dynobj_poss[n](1) + 0.5 * dynobject_size; j = j + downsample_res) {
                    for (double k = dynobj_poss[n](2) - 0.5 * dynobject_size;
                         k < dynobj_poss[n](2) + 0.5 * dynobject_size; k = k + downsample_res) {
                        if (sqrt((i - dynobj_poss[n](0)) * (i - dynobj_poss[n](0)) +
                                 (j - dynobj_poss[n](1)) * (j - dynobj_poss[n](1)) +
                                 (k - dynobj_poss[n](2)) * (k - dynobj_poss[n](2))) > 0.5 * dynobject_size) {
                            continue;
                        }
                        temp_point.x = i;
                        temp_point.y = j;
                        temp_point.z = k;
                        dynobj_points_vis.push_back(temp_point);
                        dynobj_points.push_back(temp_point);

                        dynobj_pointsindex.push_back(dynpt_count + origin_mapptcount);
                        dynpt_count++;
                    }
                }
            }
        }
        // }

        dynobj_points_vis.width = dynobj_points_vis.points.size();
        dynobj_points_vis.height = 1;
        dynobj_points_vis.is_dense = true;

        pcl::toROSMsg(dynobj_points_vis, dynobj_points_pcd);
        dynobj_points_pcd.header = odom_.header;
        dynobj_points_pcd.header.frame_id = "world";
        pub_dyncloud.publish(dynobj_points_pcd);

        kdtree_dyn.setInputCloud(dynobj_points_vis.makeShared());
        has_dyn_map = true;

        dyn_start_time = ros::Time::now();
    }
}

//得到机身点云
void multiOdometryCallbck(const nav_msgs::OdometryConstPtr & msg, int drone_id){
    Eigen::Vector3d uav_pos;
    uav_pos(0) = msg->pose.pose.position.x;
    uav_pos(1) = msg->pose.pose.position.y;
    uav_pos(2) = msg->pose.pose.position.z;
    other_uav_pos[drone_id] = uav_pos;

    // ROS_INFO("multi uav variable init success");

    otheruav_points[drone_id].clear();
    otheruav_pointsindex[drone_id].clear();
    // otheruav_points_vis[drone_id].clear();

    pcl::PointXYZI temp_point;
    int uavpt_count = 0;

    //make uav point clouds
    for (double i = uav_pos(0) - 0.5 * uav_size[0];
         i < uav_pos(0) + 0.5 * uav_size[0]; i = i + downsample_res) {
        for (double j = uav_pos(1) - 0.5 * uav_size[1];
             j < uav_pos(1) + 0.5 * uav_size[1]; j = j + downsample_res) {
            for (double k = uav_pos(2) - 0.5 * uav_size[2];
                 k < uav_pos(2) + 0.5 * uav_size[2]; k = k + downsample_res) {
                if (sqrt((i - uav_pos(0)) * (i - uav_pos(0)) +
                         (j - uav_pos(1)) * (j - uav_pos(1)) +
                         (k - uav_pos(2)) * (k - uav_pos(2))) > 0.5 * 1.732 * uav_size[0]) {
                    continue;
                }
                temp_point.x = i;
                temp_point.y = j;
                temp_point.z = k;
                otheruav_points[drone_id].push_back(temp_point);
                otheruav_pointsindex[drone_id].push_back(uavpt_count + origin_mapptcount + 100000 + drone_id*uav_points_num );
                // otheruav_points_vis[drone_id].push_back(temp_point);
                uavpt_count++;

                // ROS_INFO("UAV points = %f,%f,%f ", otheruav_points[drone_id][uavpt_count-1].x, otheruav_points[drone_id][uavpt_count-1].y, otheruav_points[drone_id][uavpt_count-1].z);

            }
        }
    }

    uav_points_num = uavpt_count;

    // otheruav_points_vis.width = otheruav_points_vis.points.size();
    // otheruav_points_vis.height = 1;
    // otheruav_points_vis.is_dense = true;

    // sensor_msgs::PointCloud2 otheruav_points_pcd;
    // pcl::toROSMsg(otheruav_points_vis, otheruav_points_pcd);
    // otheruav_points_pcd.header = odom_.header;
    // otheruav_points_pcd.header.frame_id = "world";
    // pub_uavcloud.publish(otheruav_points_pcd);
    // ROS_INFO("UAV Point counts = %d",uavpt_count);
}

void rcvOdometryCallbck(const nav_msgs::Odometry &odom) {
    /*if(!has_global_map)
      return;*/
    has_odom = true;
    odom_ = odom;

    Matrix4d body2world = Matrix4d::Identity();

    Eigen::Vector3d request_position;
    Eigen::Quaterniond pose;
    pose.x() = odom.pose.pose.orientation.x;
    pose.y() = odom.pose.pose.orientation.y;
    pose.z() = odom.pose.pose.orientation.z;
    pose.w() = odom.pose.pose.orientation.w;
    body2world.block<3, 3>(0, 0) = pose.toRotationMatrix();
    body2world(0, 3) = odom.pose.pose.position.x;
    body2world(1, 3) = odom.pose.pose.position.y;
    body2world(2, 3) = odom.pose.pose.position.z;

    // convert to cam pose
    sensor2world = body2world * sensor2body;

    //collision check
    if (has_global_map && collisioncheck_enable) {
        pcl::PointXYZI searchPoint;
        searchPoint.x = odom.pose.pose.position.x;
        searchPoint.y = odom.pose.pose.position.y;
        searchPoint.z = odom.pose.pose.position.z;
        if (_kdtreeLocalMap.radiusSearch(searchPoint, collision_range, pointIdxRadiusSearch,
                                         pointRadiusSquaredDistance) > 0) {
            ROS_ERROR("ENVIRONMENT COLLISION DETECTED!!!");
        }
        if (dynobj_enable && has_dyn_map) {
            if (kdtree_dyn.radiusSearch(searchPoint, collision_range, pointIdxRadiusSearch,
                                        pointRadiusSquaredDistance) > 0) {
                ROS_ERROR("DYNAMIC OBJECTS COLLISION DETECTED!!!");
            }
        }

    }

}

void rcvGlobalPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map) {
    if (has_global_map)
        return;

    ROS_WARN("Global Pointcloud received..");

    pcl::PointCloud<pcl::PointXYZI> cloud_input;
    pcl::fromROSMsg(pointcloud_map, cloud_input);

    _voxel_sampler.setLeafSize(downsample_res, downsample_res, downsample_res);
    _voxel_sampler.setInputCloud(cloud_input.makeShared());
    _voxel_sampler.filter(cloud_all_map);

    origin_mapptcount = cloud_all_map.points.size();

    _kdtreeLocalMap.setInputCloud(cloud_all_map.makeShared());

    normalEstimation.setInputCloud(cloud_all_map.makeShared());
    //对于每一个点都用半径为3cm的近邻搜索方式
    normalEstimation.setRadiusSearch(2.0 * downsample_res);
    //Kd_tree是一种数据结构便于管理点云以及搜索点云，法线估计对象会使用这种结构来找到最近邻点
    pcl::search::KdTree<pcl::PointXYZI>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZI>);
    normalEstimation.setSearchMethod(kdtree);
    //计算法线
    normalEstimation.compute(*all_normals);

    ROS_WARN("Normal compute finished.., mapsize = %d", origin_mapptcount);

    //trans the map into hash map
    pcl::PointXYZI pt_in, center;
    long int ind_x, ind_y, ind_z;

    //get max xyz
    pcl::PointXYZI global_mapmin;
    pcl::PointXYZI global_mapmax;
    pcl::getMinMax3D(cloud_all_map, global_mapmin, global_mapmax);
    map_min(0) = global_mapmin.x;
    map_min(1) = global_mapmin.y;
    map_min(2) = global_mapmin.z;
    map_max(0) = global_mapmax.x;
    map_max(1) = global_mapmax.y;
    map_max(2) = global_mapmax.z;

    cube_numx = floor((global_mapmax.x - global_mapmin.x) / hash_cubesize) + 20;
    cube_numy = floor((global_mapmax.y - global_mapmin.y) / hash_cubesize) + 20;
    cube_numz = floor((global_mapmax.z - global_mapmin.z) / hash_cubesize) + 200;

    env_box.vertex_min[0] = round(global_mapmin.x / hash_cubesize) * hash_cubesize - 10 * hash_cubesize;
    env_box.vertex_min[1] = round(global_mapmin.y / hash_cubesize) * hash_cubesize - 10 * hash_cubesize;
    env_box.vertex_min[2] = round(global_mapmin.z / hash_cubesize) * hash_cubesize - 100 * hash_cubesize;
    env_box.vertex_max[0] = round(global_mapmax.x / hash_cubesize) * hash_cubesize + 10 * hash_cubesize;
    env_box.vertex_max[1] = round(global_mapmax.y / hash_cubesize) * hash_cubesize + 10 * hash_cubesize;
    env_box.vertex_max[2] = round(global_mapmax.z / hash_cubesize) * hash_cubesize + 100 * hash_cubesize;
    fov_checker.Set_Env(env_box);
    fov_checker.Set_BoxLength(hash_cubesize);

    for (int i = 0; i < int(cloud_all_map.points.size()); i++) {
        pt_in = cloud_all_map.points[i];
        ind_x = (round((pt_in.x - env_box.vertex_min[0] + EPSS) / hash_cubesize));
        ind_y = (round((pt_in.y - env_box.vertex_min[1] + EPSS) / hash_cubesize));
        ind_z = (round((pt_in.z - env_box.vertex_min[2] + EPSS) / hash_cubesize));

        long int box_index = ind_x + ind_y * cube_numx + ind_z * cube_numx * cube_numy;
        point_hashmap[box_index].push_back(pt_in);
        pointindex_hashmap[box_index].push_back(i);
    }


    if (dynobj_enable == 1) {
        //dynamic objects initial pos generate
        srand((unsigned) time(NULL));
        for (int i = 0; i < dynobject_num; i++) {
            Eigen::Vector3d dyntemp_pos;
            dyntemp_pos(0) = rand() / double(RAND_MAX) * (global_mapmax.x - global_mapmin.x) + global_mapmin.x;
            dyntemp_pos(1) = rand() / double(RAND_MAX) * (global_mapmax.y - global_mapmin.y) + global_mapmin.y;
            dyntemp_pos(2) = rand() / double(RAND_MAX) * (global_mapmax.z - global_mapmin.z) + global_mapmin.z;
            dynobj_poss.push_back(dyntemp_pos);
            Eigen::Vector3d dyntemp_dir, dyntemp_dir_polar;
            dyntemp_dir_polar(0) = rand() / double(RAND_MAX) * 3.1415926;
            dyntemp_dir_polar(1) = (rand() / double(RAND_MAX) - 0.5) * 3.1415926;
            dyntemp_dir[2] = dyn_velocity * sin(dyntemp_dir_polar(1));
            dyntemp_dir[1] = dyn_velocity * cos(dyntemp_dir_polar(1)) * sin(dyntemp_dir_polar(0));
            dyntemp_dir[0] = dyn_velocity * cos(dyntemp_dir_polar(1)) * cos(dyntemp_dir_polar(0));
            dynobj_dir.push_back(dyntemp_dir);
        }
        dyn_start_time = ros::Time::now();
    }


    has_global_map = true;
}

inline void euc2polar(Eigen::Vector3d &euc_pt, float length, polar3D *polar_pt) {
    polar_pt->theta = round((atan2(euc_pt[1], euc_pt[0])) / M_PI * 180.0 / polar_resolution);
    polar_pt->fi = round((atan2(euc_pt[2], euc_pt.head<2>().norm()) / M_PI * 180.0 / polar_resolution));
    polar_pt->r = length;

}

inline void polar2euc(polar3D *polar_pt, Eigen::Vector3d &euc_pt) {
    euc_pt[2] = polar_pt->r * sin(polar_pt->fi * polar_resolution / 180.0 * M_PI);
    euc_pt[1] = polar_pt->r * cos(polar_pt->fi * polar_resolution / 180.0 * M_PI) *
                sin(polar_pt->theta * polar_resolution / 180.0 * M_PI);
    euc_pt[0] = polar_pt->r * cos(polar_pt->fi * polar_resolution / 180.0 * M_PI) *
                cos(polar_pt->theta * polar_resolution / 180.0 * M_PI);
}


int sense_count = 0;
double duration1 = 0.0;
double duration2 = 0.0;
double duration3 = 0.0;
double duration4 = 0.0;
double d5 = 0.0;
double duration_pub = 0.0;

void renderSensedPoints(const ros::TimerEvent &event) {

    ros::Time t1 = ros::Time::now();

    if (!has_global_map || !has_odom)
        return;

    Eigen::Quaterniond q;
    q.x() = odom_.pose.pose.orientation.x;
    q.y() = odom_.pose.pose.orientation.y;
    q.z() = odom_.pose.pose.orientation.z;
    q.w() = odom_.pose.pose.orientation.w;

    //trans other uav point in
    if(otheruav_points.size() > 0)
    {
        otheruav_points_inrender = otheruav_points;
        otheruav_pointsindex_inrender = otheruav_pointsindex;
    }

    Eigen::Vector3d pos(odom_.pose.pose.position.x, odom_.pose.pose.position.y, odom_.pose.pose.position.z);
    // pos << odom_.pose.pose.position.x, odom_.pose.pose.position.y, odom_.pose.pose.position.z;

    Eigen::Matrix3d rot(q.toRotationMatrix());
    // rot = q.toRotationMatrix();
    Eigen::Vector3d yaw_x(1, 0, 0);
    // yaw_x << 1,0,0;

    Eigen::Vector3d rotmyaw = rot * yaw_x;
    Eigen::Vector3d yaw_vec(rotmyaw(0), rotmyaw(1), 0);
    // yaw_vec(2) = 0;

    local_map.points.clear();
    local_map_filled.points.clear();

    Eigen::Vector3d searchpoint_pos;
    Eigen::Vector3d sensorrange_vec;
    sensorrange_vec << 0.5 * sensing_horizon, 0, 0;
    searchpoint_pos = pos + rot * sensorrange_vec;

    pcl::PointXYZ searchPoint(searchpoint_pos(0), searchpoint_pos(1), searchpoint_pos(2));
    // pcl::PointXYZI searchPoint(odom_.pose.pose.position.x, odom_.pose.pose.position.y, odom_.pose.pose.position.z);

    pointIdxRadiusSearch.clear();
    pointRadiusSquaredDistance.clear();

    //polar coordinate matrix
    // add margin
    double margin_ = 8.0;
    double vertical_fov_margined_ = vertical_fov + margin_;
    int polar_height = ceil(vertical_fov_margined_ / polar_resolution);
    double temp_yaw;
    int temp_width;
    if (is_360lidar == 1) {
        temp_yaw = 360.0;
        temp_width = ceil(360.0 / polar_resolution);
    } else {
        temp_yaw = yaw_fov + margin_;
        temp_width = ceil(temp_yaw / polar_resolution);
    }

    double yaw_fov_margined_ = temp_yaw;
    int polar_width = temp_width;

    Eigen::MatrixXf polar_matrix = Eigen::MatrixXf::Zero(polar_width, polar_height);
    polar_matrix.setConstant(2 * sensing_horizon);
    Eigen::MatrixXi polarindex_matrix = Eigen::MatrixXi::Zero(polar_width, polar_height);
    Eigen::MatrixXi polarabnormal_matrix = Eigen::MatrixXi::Zero(polar_width, polar_height);
    Eigen::MatrixXi polarpointtype_matrix = Eigen::MatrixXi::Zero(polar_width, polar_height);
    polarpointtype_matrix.setConstant(-1);
    int original_pointcount = 0;
    int pointcount = 0;
    int changepointcount = 0;

    double tan_ver_ = tan(M_PI * (vertical_fov_margined_ / 2.0) / 180.0);
    double yaw_ = (M_PI * (yaw_fov_margined_ / 2.0) / 180.0);
    double ver_threshold = cos(M_PI * (vertical_fov_margined_ / 2.0) / 180.0);
    double cover_dis = 0.55 * 1.7321 * downsample_res;// 0.707
    //compute effective range
    double effect_range = cover_dis / sin(0.5 * polar_resolution * M_PI / 180.0);
    // ROS_INFO("POLAR R = %f, EFFECT RANGE = %f",polar_pt.r,effect_range);

    double max_fov_angle;
    if (yaw_fov_margined_ > vertical_fov_margined_) {
        max_fov_angle = yaw_fov_margined_;
    } else {
        max_fov_angle = vertical_fov_margined_;
    }

    //save the points need to interline
    // vector<Eigen::Vector2i> culling_vec;
    vector<int> culling_vectheta;
    vector<int> culling_vecfi;
    vector<int> culling_kdindex;

    sense_count++;
    ros::Time t2 = ros::Time::now();
    duration1 = (duration1 + (t2 - t1).toSec());

    //pattern
    Eigen::MatrixXf pattern_matrix = Eigen::MatrixXf::Zero(polar_width, polar_height);
    pattern_matrix.setConstant(1);

    //avia pattern
    if (use_avia_pattern == 1) {
        float w1 = 763.82589;//7294.0*2.0*3.1415926/60.0
        float w2 = -488.41293788;// -4664.0*2.0*3.1415926/60.0
        // int linestep = 2;
        // double point_duration = 0.000025;
        double point_duration = 0.000004167 * 6;

        float t_duration = 1.0 / sensing_rate;
        double t_start = (t2 - start_time).toSec();

        double scale_x = 0.48 * polar_width / 2.0;
        double scale_y = 0.43 * polar_height / 2.0;

        int linestep = ceil(livox_linestep / polar_resolution);

        for (double t_i = t_start; t_i < t_start + t_duration; t_i = t_i + point_duration) {
            int x = round(scale_x * (cos(w1 * t_i) + cos(w2 * t_i))) + round(0.5 * polar_width);
            int y = round(scale_y * (sin(w1 * t_i) + sin(w2 * t_i))) + round(0.5 * polar_height);

            if (x > (polar_width - 1)) {
                x = (polar_width - 1);
            } else if (x < 0) {
                x = 0;
            }

            if (y > (polar_height - 1)) {
                y = (polar_height - 1);
            } else if (y < 0) {
                y = 0;
            }

            pattern_matrix(x, y) = 2;//real pattern
            pattern_matrix(x, y + linestep) = 2;
            pattern_matrix(x, y + 2 * linestep) = 2;
            pattern_matrix(x, y + 3 * linestep) = 2;
            pattern_matrix(x, y - linestep) = 2;
            pattern_matrix(x, y - 2 * linestep) = 2;
        }
    }

    //vlp32 pattern
    if (use_vlp32_pattern == 1) {
        int step = round((40.0) / 32 / polar_resolution);
        int bottom_pattern = -ceil(25.0 / polar_resolution);
        for (int i = 0; i < 32; i++) {
            int y = bottom_pattern + i * step + round(0.5 * polar_height);

            for (int j = 0; j < polar_width; j++) {
                if (y > (polar_height - 1)) {
                    y = (polar_height - 1);
                } else if (y < 0) {
                    y = 0;
                }
                pattern_matrix(j, y) = 2;
            }
        }
    }

    if (use_minicf_pattern == 1) {
        double point_duration = 1.0 / 200000.0;

        float t_duration = 1.0 / sensing_rate;
        double t_start = (t2 - start_time).toSec();

        double scale_x = 0.48 * polar_width / 2.0;
        double scale_y = 0.43 * polar_height / 2.0;
        double PI = 3.141519265357;

        for (double t_i = t_start; t_i < t_start + t_duration; t_i = t_i + point_duration) {
            int x = (int(-round(-62050.63 * t_i + 3.11 * cos(314159.2 * t_i) * sin(628.318 * 2 * t_i))) % 360) /
                    polar_resolution;
            int y = round(22.5 * cos(20 * PI * t_i) + 4 * cos(2 * PI / 0.006 * t_i) * cos(10000 * PI * t_i) + 22.5) /
                    polar_resolution + round(0.5 * polar_height);

            // ROS_INFO("X = %d, Y = %d",x,y);
            if (x > (polar_width - 1)) {
                x = (polar_width - 1);
            } else if (x < 0) {
                x = 0;
            }

            if (y > (polar_height - 1)) {
                y = (polar_height - 1);
            } else if (y < 0) {
                y = 0;
            }

            pattern_matrix(x, y) = 2;//real pattern
        }
    }

    ros::Time t8 = ros::Time::now();

    //try hashmap with fov checker
    vector<BoxPointType> fov_boxes;
    // rotmyaw << 0,1,0;
    // ROS_INFO("POS = %f,%f,%f, rotmyaw = %f,%f,%f", pos(0), pos(1),pos(2), rotmyaw(0),rotmyaw(1),rotmyaw(2));
    vector<pcl::PointXYZI> fov_points;
    vector<int> fov_pointsindex;

    if (is_360lidar == 1) {
        //box expand to search
        int expand_num = ceil(sensing_horizon / hash_cubesize);
        int cen_x = (round((pos(0) - env_box.vertex_min[0] + EPSS) / hash_cubesize));
        int cen_y = (round((pos(1) - env_box.vertex_min[1] + EPSS) / hash_cubesize));
        int cen_z = (round((pos(2) - env_box.vertex_min[2] + EPSS) / hash_cubesize));
        long int ind_x, ind_y, ind_z;
        for (int i = -expand_num; i <= expand_num; i++) {
            for (int j = -expand_num; j <= expand_num; j++) {
                for (int k = -expand_num; k <= expand_num; k++) {
                    ind_x = cen_x + i;
                    ind_y = cen_y + j;
                    ind_z = cen_z + k;

                    long int box_index = ind_x + ind_y * cube_numx + ind_z * cube_numx * cube_numy;
                    fov_points.insert(fov_points.end(), point_hashmap[box_index].begin(),
                                      point_hashmap[box_index].end());
                    fov_pointsindex.insert(fov_pointsindex.end(), pointindex_hashmap[box_index].begin(),
                                           pointindex_hashmap[box_index].end());
                }
            }
        }
        ROS_INFO("POINT SIZE = %d, box size = %d", fov_points.size(), fov_boxes.size());

    } else {
        fov_checker.check_fov(pos, rotmyaw, M_PI * (max_fov_angle / 2.0 + 10) / 180.0, sensing_horizon, fov_boxes);
        Eigen::Vector3d box_center_temp;
        long int ind_x, ind_y, ind_z;
        for (int i = 0; i < fov_boxes.size(); i++) {
            box_center_temp(0) =
                    (fov_boxes[i].vertex_max[0] - fov_boxes[i].vertex_min[0]) * 0.5 + fov_boxes[i].vertex_min[0];
            box_center_temp(1) =
                    (fov_boxes[i].vertex_max[1] - fov_boxes[i].vertex_min[1]) * 0.5 + fov_boxes[i].vertex_min[1];
            box_center_temp(2) =
                    (fov_boxes[i].vertex_max[2] - fov_boxes[i].vertex_min[2]) * 0.5 + fov_boxes[i].vertex_min[2];

            // if (i == 0) {
            //    ROS_INFO("Box pos = %f,%f,%f",box_center_temp(0),box_center_temp(1),box_center_temp(2));
            //    printf("box points: (%0.2f,%0.2f), (%0.2f,%0.2f), (%0.2f,%0.2f)", fov_boxes[i].vertex_min[0],fov_boxes[i].vertex_max[0],fov_boxes[i].vertex_min[1],fov_boxes[i].vertex_max[1],fov_boxes[i].vertex_min[2],fov_boxes[i].vertex_max[2]);
            // }

            // ind_x = (round((box_center_temp(0) - map_min(0) + EPSS)/hash_cubesize));
            // ind_y = (round((box_center_temp(1) - map_min(1)+ EPSS)/hash_cubesize));
            // ind_z = (round((box_center_temp(2) - map_min(2) + EPSS)/hash_cubesize));

            ind_x = (round((box_center_temp(0) - env_box.vertex_min[0] + EPSS) / hash_cubesize));
            ind_y = (round((box_center_temp(1) - env_box.vertex_min[1] + EPSS) / hash_cubesize));
            ind_z = (round((box_center_temp(2) - env_box.vertex_min[2] + EPSS) / hash_cubesize));

            pcl::PointXYZI box_pt;
            box_pt.x = box_center_temp(0);
            box_pt.y = box_center_temp(1);
            box_pt.z = box_center_temp(2);
            local_map.push_back(box_pt);

            long int box_index = ind_x + ind_y * cube_numx + ind_z * cube_numx * cube_numy;

            fov_points.insert(fov_points.end(), point_hashmap[box_index].begin(), point_hashmap[box_index].end());
            fov_pointsindex.insert(fov_pointsindex.end(), pointindex_hashmap[box_index].begin(),
                                   pointindex_hashmap[box_index].end());

            // ROS_INFO("SEE BOX: %f,%f,%f", box_center_temp(0),box_center_temp(1),box_center_temp(2));
            // std::cout << point_hashmap[box_index][0] << endl;
        }
        // ROS_INFO("POINT SIZE = %d, box size = %d", fov_points.size(), fov_boxes.size());
    }


    //add dyanmic objects
    // Step 1: add dynpoints into fov points
    if (dynobj_enable == 1) {
        fov_points.insert(fov_points.end(), dynobj_points.begin(), dynobj_points.end());
        fov_pointsindex.insert(fov_pointsindex.end(), dynobj_pointsindex.begin(), dynobj_pointsindex.end());
    }

    if (drone_num > 1) {
        // fov_points.insert(fov_points.end(), otheruav_points.begin(), otheruav_points.end());
        // fov_pointsindex.insert(fov_pointsindex.end(), otheruav_pointsindex.begin(), otheruav_pointsindex.end());

        for(int i = 0;i<drone_num;i++)
        {
            fov_points.insert(fov_points.end(), otheruav_points_inrender[i].begin(), otheruav_points_inrender[i].end());
            fov_pointsindex.insert(fov_pointsindex.end(), otheruav_pointsindex_inrender[i].begin(), otheruav_pointsindex_inrender[i].end());
        }
    }

    // ROS_INFO("Pass insert points");

    // if (_kdtreeLocalMap.radiusSearch(searchPoint, 0.6*sensing_horizon, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
    // {
    ros::Time t3 = ros::Time::now();
    duration2 = (duration2 + (t3 - t8).toSec());
    duration4 = 0.0;
//    double top_dura= 0.0;
//    double down_dura= 0.0;

    // const size_t size_ = pointIdxRadiusSearch.size();
    size_t size_ = fov_points.size();

    omp_set_num_threads(32);

//rot, polar_width, polar_height, size_, pos, vertical_fov_margined_,, yaw_vec, yaw_fov_margined_
#pragma omp parallel default (none) \
                     shared (pattern_matrix, polar_matrix, cloud_all_map, pointIdxRadiusSearch, \
                     min_raylength, downsample_res, polar_resolution, use_avia_pattern, _kdtreeLocalMap, \
                     polarindex_matrix, local_map_filled, local_map, \
                     all_normals, vertical_fov, yaw_fov, sensing_horizon, fov_points, fov_pointsindex, \
                     polarabnormal_matrix, original_pointcount, size_, pos, rot, polar_width, polar_height, cover_dis, drone_num, polarpointtype_matrix, uav_points_num, \
                     effect_range, origin_mapptcount)
    {
#pragma omp for
        for (size_t i = 0; i < size_; ++i) {
            // auto pt = cloud_all_map.points[pointIdxRadiusSearch[i]];
            auto pt = fov_points[i];
            Eigen::Vector3d pt3;
            pt3[0] = pt.x;
            pt3[1] = pt.y;
            pt3[2] = pt.z;
            auto dir = pt3 - pos;

            if (dir.norm() > sensing_horizon) {
                continue;
            }

            //this is a simple 2D realize
            // if (fabs(dir[2]) > dir.head<2>().norm() * tan_ver_)
            //     continue;

            // Eigen::Vector3d dir_norm = dir;
            // dir_norm.normalize();
            // if (dir.dot(yaw_vec) < 0)
            //     continue;

            // Vector3d dir_xy = dir;
            // dir_xy(2) = 0;
            // // ROS_INFO("YAWVEC = %f,%f, dir = %f,%f",yaw_vec(0),yaw_vec(1),dir(0),dir(1));
            // if (acos(dir.dot(yaw_vec) / (dir.norm() * yaw_vec.norm())) > yaw_)//add yaw fov  0.1944
            //     continue;

            //dead zone
            // if (dir.norm() < min_raylength)
            //     continue;

            // ROS_INFO("After min ray length");

            original_pointcount++;

            //use polar to filter the right pointclouds
            polar3D polar_pt;
            Eigen::Vector3d dir_vec;
            // dir_vec = pt3 - pos;

            //trans coordinate to lidar coordinate
            dir_vec = rot.transpose() * dir;

            euc2polar(dir_vec, dir_vec.norm(), &polar_pt);
            // ROS_INFO("dir_vec = %f,%f,%f, polar = %d,%d",dir(0),dir(1),dir(2), polar_pt.theta,polar_pt.fi);
            int cen_theta_index = polar_pt.theta + round(0.5 * polar_width);
            int cen_fi_index = polar_pt.fi + round(0.5 * polar_height);

            int half_cover_angle = ceil(
                    (asin(cover_dis / dir_vec.norm()) / (M_PI * polar_resolution / 180.0)));
            // ROS_INFO("half cover angle = %d",half_cover_angle);
            // int half_cover_angle = 1;

            if (polar_pt.r > effect_range) {

                if (polar_matrix(cen_theta_index, cen_fi_index) > polar_pt.r) {

                    polar_matrix(cen_theta_index, cen_fi_index) = polar_pt.r;
                    polarabnormal_matrix(cen_theta_index, cen_fi_index) = 1;
                    // polarindex_matrix(cen_theta_index,cen_fi_index) = pointIdxRadiusSearch[i];

                    if ((fov_pointsindex[i] >= origin_mapptcount+100000)&&drone_num > 1) {
                        int point_temp_index = fov_pointsindex[i] - origin_mapptcount - 100000;
                        int droneid_temp = (point_temp_index/(uav_points_num));
                        polarpointtype_matrix(cen_theta_index, cen_fi_index) = droneid_temp;
                    }else{
                        polarpointtype_matrix(cen_theta_index, cen_fi_index) = -1;
                    }
                }

            } else {

                for (int theta_index_o = cen_theta_index - half_cover_angle;
                     theta_index_o <= cen_theta_index + half_cover_angle; theta_index_o++) {
                    for (int fi_index_o = cen_fi_index - half_cover_angle;
                         fi_index_o <= cen_fi_index + half_cover_angle; fi_index_o++) {

                        // ros::Time t8 = ros::Time::now();

                        int theta_index = theta_index_o;
                        int fi_index = fi_index_o;
                        if (theta_index > (polar_width - 1)) {
                            continue;
                            theta_index = (polar_width - 1);
                        } else if (theta_index < 0) {
                            continue;
                            theta_index = 0;
                        }

                        if (fi_index > (polar_height - 1)) {
                            continue;
                            fi_index = (polar_height - 1);
                        } else if (fi_index < 0) {
                            continue;
                            fi_index = 0;
                        }

                        // double ray_depth = inter_point.norm();
                        double ray_depth = polar_pt.r;

                        // output_pt.x = inter_point_world(0);
                        // output_pt.y = inter_point_world(1);
                        // output_pt.z = inter_point_world(2);

                        if (polar_matrix(theta_index, fi_index) > ray_depth) {//+ sensing_horizon

                            polar_matrix(theta_index, fi_index) = ray_depth;//+ sensing_horizon
                            // polarindex_matrix(theta_index,fi_index) = pointIdxRadiusSearch[i];
                            polarindex_matrix(theta_index, fi_index) = fov_pointsindex[i];
                            // polar_matrix(theta_index, fi_index) = polar_pt.r;

                            polarabnormal_matrix(theta_index, fi_index) = -1;

                            if ((fov_pointsindex[i] >= origin_mapptcount+100000)&&drone_num > 1) {
                                int point_temp_index = fov_pointsindex[i] - origin_mapptcount - 100000;
                                int droneid_temp = (point_temp_index/(uav_points_num));
                                polarpointtype_matrix(theta_index, fi_index) = droneid_temp;
                            }else{
                                polarpointtype_matrix(theta_index, fi_index) = -1;
                            }

                        }

                    }
                }

            }   // end of for loop
        }
    }

    ros::Time t5 = ros::Time::now();
    d5 = (t5 - t3).toSec();

    // ROS_INFO("After first filter");


#pragma omp parallel default (none) \
                     shared (polarindex_matrix, culling_kdindex, polar_width, polar_height, cover_dis, rot, pos)
    {
        std::vector<int> vec_private;

#pragma omp for nowait collapse(2)
        //compute plane interline
        for (int i = 0; i < polar_width; i++) {
            for (int j = 0; j < polar_height; j++) {
                if (polarindex_matrix(i, j) != 0) {
                    //check if the index in the vector
                    if (std::find(vec_private.begin(), vec_private.end(), polarindex_matrix(i, j)) !=
                        vec_private.end()) {
                        continue;
                    } else {
                        vec_private.push_back(polarindex_matrix(i, j));
                    }

                }
            }
        }
#pragma omp critical
        culling_kdindex.insert(culling_kdindex.end(), vec_private.begin(), vec_private.end());
    }
    // ROS_INFO("CULLING COUNT = %d",culling_kdindex.size());

    //print other uav points
    // for(int i = 0;i<otheruav_points_inrender[0].size();i++)
    // {
    //     ROS_INFO("Uav 0 POINTS in render is %f,%f,%f ", otheruav_points_inrender[0][i].x, otheruav_points_inrender[0][i].y,otheruav_points_inrender[0][i].z);
    // }

// omp_set_num_threads(32);
    std::vector<int>::iterator iter;
#pragma omp parallel default (none) \
                     shared (culling_kdindex, \
                     polar_matrix, cloud_all_map, \
                     downsample_res, polar_resolution, \
                     all_normals, curvature_limit, fov_points, origin_mapptcount, dynobj_points, \
                     polarabnormal_matrix, sensing_horizon, dynobj_enable, \
                     otheruav_points, drone_num, polarpointtype_matrix, polar_width, polar_height, cover_dis, rot, pos, uav_points_num,otheruav_points_inrender\
                     )
    {
#pragma omp for

        // for(iter = culling_kdindex.begin();iter!=culling_kdindex.end();++iter)
        for (int i = 0; i < culling_kdindex.size(); i++) {

            //multi uav add
            int droneid_temp,point_temp_index;

            int point_index = culling_kdindex[i];

            //!check if it is dyn points
            pcl::PointXYZI pt;
            if (point_index >= origin_mapptcount) {
                if (dynobj_enable == 1) {
                    pt = dynobj_points[point_index - origin_mapptcount];
                } else if (drone_num > 1) {
                    point_temp_index = point_index - origin_mapptcount - 100000;
                    droneid_temp = (point_temp_index/(uav_points_num));
                    // ROS_INFO("point_temp_index = %d,uav_points_num = %d, droneid_temp = %d",point_temp_index,uav_points_num,droneid_temp);
                    // std::cout<<"droneid_temp  =  " << droneid_temp<<std::endl;
                    if((point_temp_index - droneid_temp*uav_points_num) > (otheruav_points_inrender[droneid_temp].size()-1))
                    {
                        // point_temp_index = drone_num-1;
                        continue;
                    }
                    pt = otheruav_points_inrender[droneid_temp][point_temp_index - droneid_temp*uav_points_num];
                    // ROS_INFO("Pts = %f,%f,%f",pt.x,pt.y,pt.z);
                } else {
                    // ROS_INFO("POINT INDEX = %d",point_index);
                    continue;
                }

            } else {
                pt = cloud_all_map.points[point_index];
            }

            // ROS_INFO("Take points");

            // auto pt = fov_points[point_index];
            Eigen::Vector3d pt3;
            pt3[0] = pt.x;
            pt3[1] = pt.y;
            pt3[2] = pt.z;
            auto dir = pt3 - pos;

            polar3D polar_pt;
            Eigen::Vector3d dir_vec;
            // dir_vec = pt3 - pos;

            //trans coordinate to lidar coordinate
            dir_vec = rot.transpose() * dir;

            euc2polar(dir_vec, dir_vec.norm(), &polar_pt);
            // ROS_INFO("dir_vec = %f,%f,%f, polar = %d,%d",dir(0),dir(1),dir(2), polar_pt.theta,polar_pt.fi);
            int cen_theta_index = polar_pt.theta + round(0.5 * polar_width);
            int cen_fi_index = polar_pt.fi + round(0.5 * polar_height);

            int half_cover_angle = ceil(
                    (asin(cover_dis / dir_vec.norm()) / (M_PI * polar_resolution / 180.0)));

            // ROS_INFO("Before interline");

            //!check if it is dyn points
            if (point_index >= origin_mapptcount) {
                for (int theta_index_o = cen_theta_index - half_cover_angle;
                     theta_index_o <= cen_theta_index + half_cover_angle; theta_index_o++) {
                    for (int fi_index_o = cen_fi_index - half_cover_angle;
                         fi_index_o <= cen_fi_index + half_cover_angle; fi_index_o++) {

                        // ros::Time t8 = ros::Time::now();

                        int theta_index = theta_index_o;
                        int fi_index = fi_index_o;
                        if (theta_index > (polar_width - 1)) {
                            continue;
                            theta_index = (polar_width - 1);
                        } else if (theta_index < 0) {
                            continue;
                            theta_index = 0;
                        }

                        if (fi_index > (polar_height - 1)) {
                            continue;
                            fi_index = (polar_height - 1);
                        } else if (fi_index < 0) {
                            continue;
                            fi_index = 0;
                        }

                        double ray_depth = polar_pt.r;

                        if (polar_matrix(theta_index, fi_index) > ray_depth) {

                            polar_matrix(theta_index, fi_index) = ray_depth;
                            polarabnormal_matrix(theta_index, fi_index) = -1;

                        }

                        if ((point_index >= origin_mapptcount+100000)&&drone_num > 1) {
                            polarpointtype_matrix(theta_index, fi_index) = droneid_temp;
                            // std::cout<<"droneid_temp  =  " << droneid_temp<<std::endl;
                        }else{
                            polarpointtype_matrix(theta_index, fi_index) = -1;
                        }


                    }
                }
            } else {

                Eigen::Vector3d plane_normal;
                plane_normal(0) = all_normals->points[point_index].normal_x;
                plane_normal(1) = all_normals->points[point_index].normal_y;
                plane_normal(2) = all_normals->points[point_index].normal_z;
                float curvature;
                curvature = all_normals->points[point_index].curvature;

                // ROS_INFO("Get normal");

// #pragma omp for nowait collapse(2)
                for (int theta_index_o = cen_theta_index - half_cover_angle;
                     theta_index_o <= cen_theta_index + half_cover_angle; theta_index_o++) {
                    for (int fi_index_o = cen_fi_index - half_cover_angle;
                         fi_index_o <= cen_fi_index + half_cover_angle; fi_index_o++) {

                        // ros::Time t8 = ros::Time::now();

                        int theta_index = theta_index_o;
                        int fi_index = fi_index_o;
                        if (theta_index > (polar_width - 1)) {
                            continue;
                            theta_index = (polar_width - 1);
                        } else if (theta_index < 0) {
                            continue;
                            theta_index = 0;
                        }

                        if (fi_index > (polar_height - 1)) {
                            continue;
                            fi_index = (polar_height - 1);
                        } else if (fi_index < 0) {
                            continue;
                            fi_index = 0;
                        }

                        //compute plane interline
                        polar3D cur_polarpt;
                        cur_polarpt.theta = theta_index - round(0.5 * polar_width);
                        cur_polarpt.fi = fi_index - round(0.5 * polar_height);
                        cur_polarpt.r = 1.0;
                        Eigen::Vector3d ray, inter_point, inter_point_world;
                        polar2euc(&cur_polarpt, ray);
                        ray = rot * ray;
                        ray.normalize();

                        //METHOD2 compute line and plane's cross point directly
                        double vpt = ray.dot(plane_normal);
                        // if(vpt == 0.0)// || curvature > curvature_limit
                        // {
                        //   //do not fill in free
                        //   // polar_matrix(theta_index, fi_index) = 3*sensing_horizon;
                        //   // polar_matrix(theta_index, fi_index) = polar_pt.r;
                        //   continue;
                        // }
                        double line_t = dir.dot(plane_normal) / vpt;
                        inter_point_world = pos + line_t * ray;

                        if ((inter_point_world - pt3).norm() > 0.5 * 1.7321 * downsample_res)//sqrt 3
                        {
                            // polar_matrix(theta_index, fi_index) = 3*sensing_horizon;
                            // polar_matrix(theta_index, fi_index) = polar_pt.r;
                            polarabnormal_matrix(theta_index, fi_index) = -1;
                            continue;
                        }

                        inter_point = rot.transpose() * (inter_point_world - pos);
                        double ray_depth = inter_point.norm();

                        if (polar_matrix(theta_index, fi_index) > ray_depth) {

                            polar_matrix(theta_index, fi_index) = ray_depth;
                            polarabnormal_matrix(theta_index, fi_index) = -1;
                            polarpointtype_matrix(theta_index, fi_index) = -1;

                        }

                    }
                }

                //  ROS_INFO("End one");

            }
// #pragma omp critical
        }
    }


    ros::Time t4 = ros::Time::now();
    duration3 = (t4 - t5).toSec();

    Eigen::Vector3d odom_xyz;
    odom_xyz(0) = odom_.pose.pose.position.x;
    odom_xyz(1) = odom_.pose.pose.position.y;
    odom_xyz(2) = odom_.pose.pose.position.z;

    // pcl::PointXYZI searchPoint;
    searchPoint.x = pos(0);
    searchPoint.y = pos(1);
    searchPoint.z = pos(2);

    int free_flag = 0;

    // trans polar matrix to point clouds
    for (int i = 0; i < polar_width; i++) {
        for (int j = 0; j < polar_height; j++) {
            if (use_avia_pattern == 1 || use_vlp32_pattern || use_minicf_pattern) {
                if (pattern_matrix(i, j) == 0 || pattern_matrix(i, j) == 1) {
                    continue;
                }
            }
            if (polar_matrix(i, j) < sensing_horizon) {
                Eigen::Vector3d addeuc_pt;
                polar3D polarindex_pt;
                polarindex_pt.theta = i - round(0.5 * polar_width);
                polarindex_pt.fi = j - round(0.5 * polar_height);
                polarindex_pt.r = polar_matrix(i, j);
                polar2euc(&polarindex_pt, addeuc_pt);
                addeuc_pt = rot * addeuc_pt + pos;
                pcl::PointXYZI add_pcl_pt;
                add_pcl_pt.x = addeuc_pt(0);
                add_pcl_pt.y = addeuc_pt(1);
                add_pcl_pt.z = addeuc_pt(2);

                if (polarpointtype_matrix(i, j) > -1) {
                    add_pcl_pt.intensity = (MAX_INTENSITY-MIN_INTENSITY)*(polarpointtype_matrix(i, j)+1.0)/(float(drone_num)) + MIN_INTENSITY;
                    // std::cout << "Other uav intensity is " <<  add_pcl_pt.intensity << " Point type is " << polarpointtype_matrix(i, j) << std::endl;
                } else {
                    add_pcl_pt.intensity = MIN_INTENSITY;
                }

                // filter the margin
                Eigen::Vector3d pt3;
                pt3[0] = add_pcl_pt.x;
                pt3[1] = add_pcl_pt.y;
                pt3[2] = add_pcl_pt.z;
                auto dir = pt3 - pos;

                if (dir.norm() < min_raylength)
                    continue;

                if (is_360lidar == 1) {

                } else {
                    if (acos(dir.dot(rotmyaw) / (dir.norm() * rotmyaw.norm()) > M_PI * (vertical_fov / 2.0) / 180.0)) {
                        continue;
                    }
                }


                // push back the filtered points
                local_map_filled.push_back(add_pcl_pt);
            } else if (polar_matrix(i, j) >= 1.95 * sensing_horizon) {
                Eigen::Vector3d addeuc_pt;
                polar3D polarindex_pt;
                polarindex_pt.theta = i - round(0.5 * polar_width);
                polarindex_pt.fi = j - round(0.5 * polar_height);
                polarindex_pt.r = 2 * sensing_horizon;//polar_matrix(i,j)
                polar2euc(&polarindex_pt, addeuc_pt);
                addeuc_pt = rot * addeuc_pt + odom_xyz;//.transpose()*
                pcl::PointXYZI add_pcl_pt;
                add_pcl_pt.x = addeuc_pt(0);
                add_pcl_pt.y = addeuc_pt(1);
                add_pcl_pt.z = addeuc_pt(2);
                add_pcl_pt.intensity = MIN_INTENSITY;
                // filter the margin
                Eigen::Vector3d pt3;
                pt3[0] = add_pcl_pt.x;
                pt3[1] = add_pcl_pt.y;
                pt3[2] = add_pcl_pt.z;
                auto dir = pt3 - pos;

                if (dir.norm() < min_raylength)
                    continue;

                if (is_360lidar == 1) {

                } else {
                    if (acos(dir.dot(rotmyaw) / (dir.norm() * rotmyaw.norm()) > M_PI * (vertical_fov / 2.0) / 180.0)) {
                        continue;
                    }
                }

                free_flag = 1;
//                local_map_filled.push_back(add_pcl_pt);
            }
        }
    }

    if (free_flag == 1) {
//      ROS_ERROR("Give free far !!!!!!!!!!!!!!!!!!!!!!!");
    }

    ROS_INFO("GET OUT OF LOOP, pointcount = %d, origin pointcount = %d, change point = %d",
             local_map_filled.points.size(), original_pointcount, changepointcount);

    local_map.width = local_map.points.size();
    local_map.height = 1;
    local_map.is_dense = true;

    local_map_filled.width = local_map_filled.points.size();
    local_map_filled.height = 1;
    local_map_filled.is_dense = true;

    pcl::toROSMsg(local_map_filled, local_map_pcd);

    local_map_pcd.header = odom_.header;
    local_map_pcd.header.frame_id = "world";
    pub_cloud.publish(local_map_pcd);


    // transform
    std::string sensor_frame_id_ = "sensor";
    static tf2_ros::TransformBroadcaster br;
    geometry_msgs::TransformStamped transform;
    ros::Time time_stamp_ = odom_.header.stamp;
    transform.header.stamp = time_stamp_;
    transform.header.frame_id = "world";
    transform.child_frame_id = sensor_frame_id_;

    transform.transform.translation.x = pos.x();
    transform.transform.translation.y = pos.y();
    transform.transform.translation.z = pos.z();
    transform.transform.rotation.x = q.x();
    transform.transform.rotation.y = q.y();
    transform.transform.rotation.z = q.z();
    transform.transform.rotation.w = q.w();

    br.sendTransform(transform);

    Eigen::Matrix4d sensor2world;
    sensor2world <<
                 rot(0, 0), rot(0, 1), rot(0, 2), pos.x(),
            rot(1, 0), rot(1, 1), rot(1, 2), pos.y(),
            rot(2, 0), rot(2, 1), rot(2, 2), pos.z(),
            0, 0, 0, 1;
    Eigen::Matrix4d world2sensor;
    world2sensor = sensor2world.inverse();

    // pointcloud
    point_in_sensor.points.clear();
    // pcl::copyPointCloud(local_map_filled, point_in_sensor);

    //将观测到队友的点云添加进sensor_cloud
//    pcl::PointCloud<pcl::PointXYZI> local_map_uav;
//    pcl::copyPointCloud(local_map_filled, local_map_uav);
//    for (int i = 0; i < otheruav_points_vis.size(); i++) {
//        local_map_uav.push_back(otheruav_points_vis.points[i]);
//    }

    pcl::transformPointCloud(local_map_filled, point_in_sensor, world2sensor);
    point_in_sensor.width = point_in_sensor.points.size();
    point_in_sensor.height = 1;
    point_in_sensor.is_dense = true;
    pcl::toROSMsg(point_in_sensor, sensor_map_pcd);
    sensor_map_pcd.header.frame_id = sensor_frame_id_;
    sensor_map_pcd.header.stamp = time_stamp_;
    pub_intercloud.publish(sensor_map_pcd);

    ros::Time t7 = ros::Time::now();
    duration_pub = (duration_pub + (t7 - t4).toSec());

    // }
    ros::Time t_total = ros::Time::now();
    // std::cout << "ONE FOV GENERATE:"<<  (t_total-t1).toSec()<<std::endl;;
    // ROS_INFO("Time statics: %f,%f,%f,%f,%f total_time = %f",duration1/sense_count,duration2/sense_count,duration3/sense_count,duration4/pointIdxRadiusSearch.size(),duration_pub/sense_count,(t_total-t1).toSec());
    ROS_INFO("Time statics: %f,%f,%f,%f total_time = %f", duration1 / sense_count, duration2 / sense_count, duration3,
             duration_pub / sense_count, (t_total - t1).toSec());

}


void pubSensorPose(const ros::TimerEvent &e) {
    Eigen::Quaterniond q;
    q = sensor2world.block<3, 3>(0, 0);

    geometry_msgs::PoseStamped sensor_pose;
    sensor_pose.header = odom_.header;
    sensor_pose.header.frame_id = "/map";
    sensor_pose.pose.position.x = sensor2world(0, 3);
    sensor_pose.pose.position.y = sensor2world(1, 3);
    sensor_pose.pose.position.z = sensor2world(2, 3);
    sensor_pose.pose.orientation.w = q.w();
    sensor_pose.pose.orientation.x = q.x();
    sensor_pose.pose.orientation.y = q.y();
    sensor_pose.pose.orientation.z = q.z();
    pub_pose.publish(sensor_pose);
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "pcl_render");
    ros::NodeHandle nh("~");

    nh.param("quadrotor_name", quad_name, std::string("quadrotor"));

    nh.getParam("is_360lidar", is_360lidar);
    nh.getParam("sensing_horizon", sensing_horizon);
    nh.getParam("sensing_rate", sensing_rate);
    nh.getParam("estimation_rate", estimation_rate);
    nh.getParam("polar_resolution", polar_resolution);
    nh.getParam("yaw_fov", yaw_fov);
    nh.getParam("vertical_fov", vertical_fov);
    nh.getParam("min_raylength", min_raylength);
    nh.getParam("downsample_res", downsample_res);
    nh.getParam("livox_linestep", livox_linestep);
    nh.getParam("use_avia_pattern", use_avia_pattern);
    nh.getParam("curvature_limit", curvature_limit);
    nh.getParam("hash_cubesize", hash_cubesize);
    nh.getParam("use_vlp32_pattern", use_vlp32_pattern);
    nh.getParam("use_minicf_pattern", use_minicf_pattern);

    //dyn parameters
    nh.getParam("dynobj_enable", dynobj_enable);
    nh.getParam("dynobject_size", dynobject_size);
    nh.getParam("dynobject_num", dynobject_num);
    nh.getParam("dyn_mode", dyn_mode);
    nh.getParam("dyn_velocity", dyn_velocity);

    nh.getParam("collisioncheck_enable", collisioncheck_enable);
    nh.getParam("collision_range", collision_range);

    nh.getParam("output_pcd", output_pcd);

    nh.getParam("map/x_size", x_size);
    nh.getParam("map/y_size", y_size);
    nh.getParam("map/z_size", z_size);

    //subscribe other uav pos
    // int drone_num = 3;
    nh.param("uav_num", drone_num, 1);
    nh.param("drone_id", drone_id, 0);
    other_uav_pos.resize(drone_num);
    otheruav_points.resize(drone_num);
    otheruav_pointsindex.resize(drone_num);
    // uav_num = drone_num;
    // /quad_1/lidar_slam/odom
    ros::Subscriber* subs = new ros::Subscriber[drone_num];
    for(int i = 0 ; i < drone_num ; i++){
        if(i == drone_id){continue;}
        string topic = "/quad_";
        topic+= to_string(i);
        topic+= "/lidar_slam/odom";
        cout<<topic<<endl;
        subs[i] = nh.subscribe<nav_msgs::Odometry>(topic, 50, boost::bind(&multiOdometryCallbck, _1, i));
    }
    uav_size[0] = 0.5;
    uav_size[1] = 0.5;
    uav_size[2] = 0.5;
    uav_points_num = 1*ceil(uav_size[0] * uav_size[1] * uav_size[2] / (downsample_res*downsample_res*downsample_res));
    ROS_INFO("Uav points num = %d", uav_points_num);


    // subscribe point cloud
    global_map_sub = nh.subscribe("global_map", 1, rcvGlobalPointCloudCallBack);
    odom_sub = nh.subscribe("odometry", 50, rcvOdometryCallbck);

    // publisher depth image and color image
    pub_dyncloud = nh.advertise<sensor_msgs::PointCloud2>("dyn_cloud", 10);
    pub_intercloud = nh.advertise<sensor_msgs::PointCloud2>("sensor_cloud", 10);
    pub_cloud = nh.advertise<sensor_msgs::PointCloud2>("cloud", 10);
    pub_pose = nh.advertise<geometry_msgs::PoseStamped>("sensor_pose", 10);
    pub_uavcloud = nh.advertise<sensor_msgs::PointCloud2>("uav_cloud", 10); //扫描机身的点云
    double sensing_duration = 1.0 / sensing_rate;
    double estimate_duration = 1.0 / estimation_rate;

    start_time = ros::Time::now();

    local_sensing_timer = nh.createTimer(ros::Duration(sensing_duration), renderSensedPoints);
    pose_timer = nh.createTimer(ros::Duration(estimate_duration), pubSensorPose);
    dynobj_timer = nh.createTimer(ros::Duration(sensing_duration), dynobjGenerate);

    pkg_path = ros::package::getPath("octomap_server");
    pkg_path.append("/data/" + quad_name + "_explore_persentage.txt");
    std::cout << "\nFound pkg_path = " << pkg_path << std::endl;
    myfile.open(pkg_path.c_str(), std::ios_base::out);//, std::ios_base::out

    inv_resolution = 1.0 / resolution;
    gl_xl = -x_size / 2.0;
    gl_yl = -y_size / 2.0;
    gl_zl = 0.0;
    GLX_SIZE = (int) (x_size * inv_resolution);
    GLY_SIZE = (int) (y_size * inv_resolution);
    GLZ_SIZE = (int) (z_size * inv_resolution);

    sensor2body << 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    ros::Rate rate(100);
    bool status = ros::ok();
    while (status) {
        ros::spinOnce();
        status = ros::ok();
        rate.sleep();
    }

    if (output_pcd == 1) {
        //write pcd
        pkg_path = ros::package::getPath("octomap_server");
        std::string pcd_name(quad_name + "cloud_explored");
        if (pcl::io::savePCDFileASCII(pkg_path + "/data/" + pcd_name + ".pcd", local_map) >= 0) {
            std::cerr << "Saved  " << pcd_name << ".pcd" << std::endl;
        }
    }

}