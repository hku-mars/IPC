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
#include "bvp/obvp_solver.hpp"

using namespace fmt;
using namespace std;
class RandomMap {
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

    bvp::ObvpSolver::Ptr bvp_;

    struct RandomUtils {
        default_random_engine eng;
        uniform_real_distribution<double> unif_pm_1;
        uniform_real_distribution<double> unif_0_1;
    } rand_;

    struct ConfigParams {
        int cylinder_num, circle_num, spline_num, pillar_num;
        double cylinder_density;
        double map_size_x, map_size_y, map_size_z, resolution;
        double init_x, init_y, init_z;
        string frame_id;
        double ground_height;

        vector<Eigen::Vector3d> clear_area;
        vector<double> clear_radius;
        // ground params
        double ground_size_x,ground_size_y,ground_size_z;
        double ground_resolution;
        bool use_ground{false};
        bool visual_no_groud_map{false};

        // all params has [0] as min and [1] as max
        // cylinder params
        double cylinder_radius[2];
        double cylinder_height[2];
        double cylinder_angle_resolution;
        vector<Eigen::Vector3d> fix_cylinder;

        // pillar params
        double pillar_width[2];
        double pillar_height[2];
        double pillar_length[2];

        // circle parames
        double circle_height[2];
        double circle_radius[2];
        double circle_thickness[2];


        // spline params
        double bending_rate;
        double height_max, height_min;
        double radius_max, radius_min;
        double sample_theta_reso;
        double sample_t_reso;
        int max_branch_num;
        int total_num;
        double tree_width;

    } cfg_;

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
            for(int i = 0 ; i< param_value.size();i++){
                cout << param_value[i] <<" ";
            }cout<<endl;
        } else {
            print("\033[0;33m Load param {} failed, use default value: \033[0;0m", param_name);
            param_value = default_value;
            for(int i = 0 ; i< param_value.size();i++){
                cout << param_value[i] <<" ";
            }cout<<endl;
        }
    }

public:
    RandomMap(ros::NodeHandle &nh) {
        nh_ = nh;

        rand_.unif_pm_1 = uniform_real_distribution<double>(-1, 1);
        rand_.unif_0_1 = uniform_real_distribution<double>(0, 1);

        LoadParam("map/x_size", cfg_.map_size_x, 20.0);
        LoadParam("map/y_size", cfg_.map_size_y, 20.0);
        LoadParam("map/z_size", cfg_.map_size_z, 5.0);
        LoadParam("map/ground_height", cfg_.ground_height, -1.0);
        LoadParam("map/init_x", cfg_.init_x, 0.0);
        LoadParam("map/init_y", cfg_.init_y, 0.0);
        LoadParam("map/init_z", cfg_.init_z, 0.0);
        LoadParam("map/resolution", cfg_.resolution, 0.1);
        double hz;
        LoadParam("map/pub_freq", hz, 1.0);
        LoadParam<string>("map/frame_id", cfg_.frame_id, "world");
        int seed;
        LoadParam("map/random_seed", seed, -1);
        rand_.eng = default_random_engine (seed);
        // Load ground parameters
        LoadParam("ground/enable", cfg_.use_ground, true);
        LoadParam("ground/visual_no_ground", cfg_.visual_no_groud_map, true);
        LoadParam("ground/x_size", cfg_.ground_size_x, cfg_.map_size_x);
        LoadParam("ground/y_size", cfg_.ground_size_y, cfg_.map_size_y);
        LoadParam("ground/z_size", cfg_.ground_size_z, cfg_.map_size_z);
        LoadParam("ground/resolution", cfg_.ground_resolution, cfg_.resolution);

        // Load cylinder parameters
        LoadParam("cylinder/num", cfg_.cylinder_num, 20);
        LoadParam("cylinder/density", cfg_.cylinder_density, -1.0);

        LoadParam("cylinder/angle_resolution", cfg_.cylinder_angle_resolution, 0.1);
        LoadParam("cylinder/height_min", cfg_.cylinder_height[0], 1.0);
        LoadParam("cylinder/height_max", cfg_.cylinder_height[1], 3.0);
        LoadParam("cylinder/radius_min", cfg_.cylinder_radius[0], 0.3);
        LoadParam("cylinder/radius_max", cfg_.cylinder_radius[1], 0.5);


        // Load circle parameters
        LoadParam("circle/num", cfg_.circle_num, 20);
        LoadParam("circle/height_min", cfg_.circle_height[0], 1.0);
        LoadParam("circle/height_max", cfg_.circle_height[1], 3.0);
        LoadParam("circle/radius_min", cfg_.circle_radius[0], 0.3);
        LoadParam("circle/radius_max", cfg_.circle_radius[1], 0.5);
        LoadParam("circle/thickness_min", cfg_.circle_thickness[0], 0.1);
        LoadParam("circle/thickness_max", cfg_.circle_thickness[1], 0.1);

        // Load pillar parameters
        LoadParam("pillar/num", cfg_.pillar_num, 20);
        LoadParam("pillar/height_min", cfg_.pillar_height[0], 1.0);
        LoadParam("pillar/height_max", cfg_.pillar_height[1], 3.0);
        LoadParam("pillar/width_min", cfg_.pillar_width[0], 0.3);
        LoadParam("pillar/width_max", cfg_.pillar_width[1], 0.5);
        LoadParam("pillar/length_min", cfg_.pillar_length[0], 0.3);
        LoadParam("pillar/length_max", cfg_.pillar_length[1], 0.5);

        //Load spline parameters
//        LoadParam("spline/num", cfg_.spline_num, 5);

        LoadParam("spline/num", cfg_.spline_num, 5);
        LoadParam("spline/radius_min", cfg_.radius_min, 0.1);
        LoadParam("spline/radius_max", cfg_.radius_max, 0.2);
        LoadParam("spline/bending_rate", cfg_.bending_rate, 0.5);
        LoadParam("spline/height_max", cfg_.height_max, 4.0);
        LoadParam("spline/height_min", cfg_.height_min, 3.0);
        LoadParam("spline/width", cfg_.tree_width, 3.0);
        LoadParam("spline/sample_t_reso", cfg_.sample_t_reso, 0.05);
        LoadParam("spline/sample_theta_reso", cfg_.sample_theta_reso, 0.5);
        LoadParam("spline/max_branch_num", cfg_.max_branch_num, 3);

        vector<double> clear_x, clear_y;
        LoadParam("map/clear_radius",cfg_.clear_radius,cfg_.clear_radius);
        LoadParam("map/clear_area/x",clear_x,clear_x);
        LoadParam("map/clear_area/y",clear_y,clear_y);
        for(size_t i = 0; i < clear_x.size();i++){
            cfg_.clear_area.push_back(Eigen::Vector3d(clear_x[i], clear_y[i],0));
        }

        vector<double> fix_x, fix_y;
        LoadParam("cylinder/fix/x",fix_x,fix_x);
        LoadParam("cylinder/fix/y",fix_y,fix_y);
        for(size_t i = 0; i < fix_x.size();i++){
            cfg_.fix_cylinder.push_back(Eigen::Vector3d(fix_x[i], fix_y[i],0));
        }

        pc_pub_ = nh_.advertise<sensor_msgs::PointCloud2>("/map_generator/global_cloud", 1);
        pc_no_ground_pub_ = nh_.advertise<sensor_msgs::PointCloud2>("/map_generator/visual_cloud", 1);
        pub_timer_ = nh_.createTimer(ros::Duration(1.0 / hz), &RandomMap::PubTimerCallback, this);
        GenerateRandomMap();
    }

    ~RandomMap() {}

private:

    void PubTimerCallback(const ros::TimerEvent &e) {
        if (!map_is_ready) {
            return;
        }
        pc2_map_.header.stamp = ros::Time::now();
        pc_pub_.publish(pc2_map_);

        if(cfg_.visual_no_groud_map){
            pc2_no_ground_map_.header.stamp = ros::Time::now();
            pc_no_ground_pub_.publish(pc2_no_ground_map_);
        }
    }

    void GenerateRandomMap() {
        GeneratePillars();
        GenerateCylinders();
        GenerateCircles();
        GenerateSpline();
        if(cfg_.visual_no_groud_map){
            pcl_no_ground_map_ = pcl_map_;
            pcl_no_ground_map_.width = pcl_no_ground_map_.points.size();
            pcl_no_ground_map_.height = 1;
            pcl_no_ground_map_.is_dense = true;
            pcl::toROSMsg(pcl_no_ground_map_, pc2_no_ground_map_);
            pc2_no_ground_map_.header.frame_id = cfg_.frame_id;
        }

        GenerateGround();
        pcl_map_.width = pcl_map_.points.size();
        pcl_map_.height = 1;
        pcl_map_.is_dense = true;
        pcl::toROSMsg(pcl_map_, pc2_map_);
        pc2_map_.header.frame_id = cfg_.frame_id;
        map_is_ready = true;
    }

    void GenerateCircles() {
        pcl::PointXYZ pt_random;
        uniform_real_distribution<double> rand_roll = uniform_real_distribution<double>(-M_PI, +M_PI);
        uniform_real_distribution<double> rand_pitch = uniform_real_distribution<double>(+M_PI / 4.0, +M_PI / 2.0);
        uniform_real_distribution<double> rand_yaw = uniform_real_distribution<double>(+M_PI / 4.0, +M_PI / 2.0);
        uniform_real_distribution<double> rand_ellipse_c = uniform_real_distribution<double>(0.5, 2.0);
        for (int i = 0; i < cfg_.circle_num; i++) {
            double x0, y0, z0, rmin, rmax;
            std::vector<Eigen::Vector3d> circle_set;

            x0 = rand_.unif_pm_1(rand_.eng) * cfg_.map_size_x / 2;
            y0 = rand_.unif_pm_1(rand_.eng) * cfg_.map_size_y / 2;
            bool conti{false};
            for(size_t i = 0 ; i< cfg_.clear_area.size();i++){
                Eigen::Vector3d cur(x0,y0,0);
                if((cur-cfg_.clear_area[i]).norm() < cfg_.clear_radius[i]){
                    conti = true;
                    break;
                }
            }
            if(conti){
                i--;
                continue;
            }

            z0 = rand_.unif_0_1(rand_.eng) * (cfg_.circle_height[1] - cfg_.circle_height[0]) + cfg_.circle_height[0];
            rmin = rand_.unif_0_1(rand_.eng) * (cfg_.circle_radius[1] - cfg_.circle_radius[0]) +
                    cfg_.circle_radius[0];
            rmax = rmin + rand_.unif_0_1(rand_.eng) * (cfg_.circle_thickness[1] - cfg_.circle_thickness[0]) +
                    cfg_.circle_thickness[0];


            double a, b;
            a = rand_ellipse_c(rand_.eng);
            b = rand_ellipse_c(rand_.eng);

            double x, y, z;
            Eigen::Vector3d pt3, pt3_rot;
            for (double theta = -M_PI; theta < M_PI; theta += 0.025) {
                for(double R = rmin; R < rmax; R += cfg_.resolution){
                    x = a * cos(theta) * R;
                    y = b * sin(theta) * R;
                    z = 0;
                    pt3 << x, y, z;
                    circle_set.push_back(pt3);
                }
            }
            // Define a random 3d rotation matrix
            Eigen::Matrix3d Rot;
            double roll, pitch, yaw;
            double alpha, beta, gama;
            roll = rand_roll(rand_.eng); // alpha
            pitch = rand_pitch(rand_.eng); // beta
            yaw = rand_yaw(rand_.eng); // gama

            alpha = roll;
            beta = pitch;
            gama = yaw;

            double p = rand_.unif_0_1(rand_.eng);
            if (p < 0.5) {
                beta = M_PI / 2.0;
                gama = M_PI / 2.0;
            }

            Rot << cos(alpha) * cos(gama) - cos(beta) * sin(alpha) * sin(gama), -cos(beta) * cos(gama) * sin(alpha) -
                                                                                cos(alpha) * sin(gama), sin(alpha) *
                                                                                                        sin(beta),
                    cos(gama) * sin(alpha) + cos(alpha) * cos(beta) * sin(gama), cos(alpha) * cos(beta) * cos(gama) -
                                                                                 sin(alpha) * sin(gama), -cos(alpha) *
                                                                                                         sin(beta),
                    sin(beta) * sin(gama), cos(gama) * sin(beta), cos(beta);

            for (auto pt: circle_set) {
                pt3_rot = Rot * pt;
                pt_random.x = pt3_rot(0) + x0 + 0.001;
                pt_random.y = pt3_rot(1) + y0 + 0.001;
                pt_random.z = pt3_rot(2) + z0 + 0.001;

                if (pt_random.z >= 0.0){
                    pt_random.z+=  cfg_.ground_height;
                    pcl_map_.points.push_back(pt_random);
                }

            }
        }
        print(fg(color::light_sea_green), " [Rings] Generate success.\n");
    }

    void GenerateCylinders() {
        pcl::PointXYZ pt_random;
        bool is_kdtree_empty = false;
        if (pcl_map_.points.size() > 0)
            kdtreeMap.setInputCloud(pcl_map_.makeShared());
        else
            is_kdtree_empty = true;

        if(cfg_.cylinder_density > 0){
            cfg_.cylinder_num = cfg_.map_size_x * cfg_.map_size_y / cfg_.cylinder_density;
        }
        print(" -- [INFO] Cylinder number: {}.\n", cfg_.cylinder_num);
        for (size_t i = 0; i < cfg_.cylinder_num; i++) {
            double x, y, r, h;
            x = rand_.unif_pm_1(rand_.eng) * cfg_.map_size_x / 2;
            y = rand_.unif_pm_1(rand_.eng) * cfg_.map_size_y / 2;
            bool conti{false};
            for(size_t i = 0 ; i< cfg_.clear_area.size();i++){
                Eigen::Vector3d cur(x,y,0);
                if((cur-cfg_.clear_area[i]).norm() < cfg_.clear_radius[i]){
                    conti = true;
                    break;
                }
            }
            if(conti){
                i--;
                continue;
            }

            r = rand_.unif_0_1(rand_.eng) * (cfg_.cylinder_radius[1] - cfg_.cylinder_radius[0]) +
                cfg_.cylinder_radius[0];
            h = rand_.unif_0_1(rand_.eng) * (cfg_.cylinder_height[1] - cfg_.cylinder_height[0]) +
                cfg_.cylinder_height[0];

            pcl::PointXYZ searchPoint(x, y, (cfg_.pillar_height[0] + cfg_.pillar_height[1]) / 2.0);
            pointIdxSearch.clear();
            pointSquaredDistance.clear();

            if (is_kdtree_empty == false) {
                if (kdtreeMap.nearestKSearch(searchPoint, 1, pointIdxSearch, pointSquaredDistance) > 0) {
                    if (sqrt(pointSquaredDistance[0]) < 1.0) {
                        i--;
                        continue;
                    }
                }
            }

            x = floor(x / cfg_.resolution) * cfg_.resolution + cfg_.resolution / 2.0;
            y = floor(y / cfg_.resolution) * cfg_.resolution + cfg_.resolution / 2.0;


            for (double theta = 0; theta < M_PI * 2; theta += cfg_.cylinder_angle_resolution) {
                for (double radius = r-2*cfg_.resolution; radius < r; radius += cfg_.resolution) {
                    for (double hei = 0.0; hei < h; hei += cfg_.resolution) {
                        pt_random.x = x + cos(theta) * radius;
                        pt_random.y = y + sin(theta) * radius;
                        pt_random.z = hei + cfg_.ground_height;
                        pcl_map_.points.push_back(pt_random);
                    }
                }
            }
        }

        for(size_t i = 0 ; i < cfg_.fix_cylinder.size(); i++){
            double x, y, r, h;
            x = cfg_.fix_cylinder[i].x();
            y = cfg_.fix_cylinder[i].y();

            r = rand_.unif_0_1(rand_.eng) * (cfg_.cylinder_radius[1] - cfg_.cylinder_radius[0]) +
                    cfg_.cylinder_radius[0];
            h = rand_.unif_0_1(rand_.eng) * (cfg_.cylinder_height[1] - cfg_.cylinder_height[0]) +
                    cfg_.cylinder_height[0];

            x = floor(x / cfg_.resolution) * cfg_.resolution + cfg_.resolution / 2.0;
            y = floor(y / cfg_.resolution) * cfg_.resolution + cfg_.resolution / 2.0;


            for (double theta = 0; theta < M_PI * 2; theta += cfg_.cylinder_angle_resolution) {
                for (double radius = 0.0; radius < r; radius += cfg_.resolution) {
                    for (double hei = 0.0; hei < h; hei += cfg_.resolution) {
                        pt_random.x = x + cos(theta) * radius;
                        pt_random.y = y + sin(theta) * radius;
                        pt_random.z = hei + cfg_.ground_height;
                        pcl_map_.points.push_back(pt_random);
                    }
                }
            }
        }

        print(fg(color::light_sea_green), " [Cylinder] Generate success.\n");
    }

    void GeneratePillars() {
        pcl::PointXYZ pt_random;
        bool is_kdtree_empty = false;
        if (pcl_map_.points.size() > 0)
            kdtreeMap.setInputCloud(pcl_map_.makeShared());
        else
            is_kdtree_empty = true;

        for (size_t i = 0; i < cfg_.pillar_num; i++) {
            double x, y, w, h, l;
            x = rand_.unif_pm_1(rand_.eng) * cfg_.map_size_x / 2;
            y = rand_.unif_pm_1(rand_.eng) * cfg_.map_size_y / 2;
            bool conti{false};
            for(size_t i = 0 ; i< cfg_.clear_area.size();i++){
                Eigen::Vector3d cur(x,y,0);
                if((cur-cfg_.clear_area[i]).norm() < cfg_.clear_radius[i]){
                    conti = true;
                    break;
                }
            }
            if(conti){
                i--;
                continue;
            }

            w = rand_.unif_0_1(rand_.eng) * (cfg_.pillar_width[1] - cfg_.pillar_width[0]) + cfg_.pillar_width[0];
            l = rand_.unif_0_1(rand_.eng) * (cfg_.pillar_length[1] - cfg_.pillar_length[0]) + cfg_.pillar_length[0];

            pcl::PointXYZ searchPoint(x, y, (cfg_.pillar_height[0] + cfg_.pillar_height[1]) / 2.0);
            pointIdxSearch.clear();
            pointSquaredDistance.clear();

            if (is_kdtree_empty == false) {
                if (kdtreeMap.nearestKSearch(searchPoint, 1, pointIdxSearch, pointSquaredDistance) > 0) {
                    if (sqrt(pointSquaredDistance[0]) < 1.0) {
                        i--;
                        continue;
                    }

                }
            }

            x = floor(x / cfg_.resolution) * cfg_.resolution + cfg_.resolution / 2.0;
            y = floor(y / cfg_.resolution) * cfg_.resolution + cfg_.resolution / 2.0;

            int widNum = ceil(w / cfg_.resolution);
            int lenNum = ceil(l / cfg_.resolution);

            for (int r = -widNum / 2.0; r < widNum / 2.0; r++) {
                for (int s = -lenNum / 2.0; s < lenNum / 2.0; s++) {
                    h = rand_.unif_0_1(rand_.eng) * (cfg_.pillar_height[1] - cfg_.pillar_height[0]) +
                        cfg_.cylinder_height[0];
                    int heiNum = 2.0 * ceil(h / cfg_.resolution);
                    for (int t = 0; t < heiNum; t++) {
                        pt_random.x = x + (r + 0.0) * cfg_.resolution + 0.001;
                        pt_random.y = y + (s + 0.0) * cfg_.resolution + 0.001;
                        pt_random.z = (t + 0.0) * cfg_.resolution * 0.5 + 0.001 + cfg_.ground_height;
                        pcl_map_.points.push_back(pt_random);
                    }
                }
            }
        }

        print(fg(color::light_sea_green), " [Pillar] Generate success.\n");
    }

    void GenerateGround(){
        if(!cfg_.use_ground){return;}

        pcl::PointXYZ pt_random;
        for (double i = -cfg_.ground_size_x / 2.0; i <cfg_.ground_size_x / 2.0; i += cfg_.ground_resolution)
            for (double j = -cfg_.ground_size_y / 2.0; j <cfg_.ground_size_y / 2.0; j += cfg_.ground_resolution) {
                pt_random.x = i;
                pt_random.y = j;
                pt_random.z = cfg_.ground_height;
                pcl_map_.points.push_back(pt_random);
                if(cfg_.ground_size_z == 0){
                    continue;
                }
                pt_random.z = cfg_.ground_size_z;
                pcl_map_.points.push_back(pt_random);
            }
        print(fg(color::light_sea_green), " [Groune] Generate success.\n");
    }

    void GenerateTree(pcl::PointCloud<pcl::PointXYZ> &input_map, bvp::Vec3 &center,
                      bvp::Vec3 &start_dir, bvp::Vec3 &end_dir, bvp::Vec3 &end_pos, int &b_num) {
        bvp_.reset(new bvp::ObvpSolver(1000));
        uniform_real_distribution<double> rand_b_xy = uniform_real_distribution<double>(-cfg_.tree_width, cfg_.tree_width);
        uniform_real_distribution<double> rand_bt = uniform_real_distribution<double>(0, 0.3);

        bvp::StatePVAJ start_state;
        bvp::StatePVAJ end_state;
        bvp::Vec3 zero3(0,0,0);
        start_state << zero3, start_dir * cfg_.bending_rate, zero3, zero3;
        end_state << end_pos, end_dir * cfg_.bending_rate, zero3, zero3;

        /* Generate Main Branch */
        bvp::Piece main_branch = bvp_->GenFixStateMinJerkOptT(start_state, end_state);
        double main_dur = main_branch.getDuration();
        double max_vel = main_branch.getMaxVelRate();
        double step = (cfg_.radius_max - cfg_.radius_min) / max_vel;
        for (double eval_t = 0.0; eval_t < main_dur; eval_t += cfg_.sample_t_reso) {
            double radius = step * (max_vel - main_branch.getVel(eval_t).norm()) + cfg_.radius_min;
            bvp::Vec3 pos = main_branch.getPos(eval_t);
            for (double theta = 0.0; theta < 2 * M_PI; theta += cfg_.sample_theta_reso) {
                pcl::PointXYZ p;
                p.x = center.x() + pos.x() + sin(theta) * radius;
                p.y = center.y() + pos.y() + cos(theta) * radius;
                p.z = pos.z() + cfg_.ground_height;
                input_map.points.push_back(p);
            }
        }
        bvp::StatePVAJ branch_start, branch_end;

        for (int i = 0; i < b_num; i++) {
            double split_t = rand_bt(rand_.eng) * main_dur;
            bvp::Vec3 endp(rand_b_xy(rand_.eng), rand_b_xy(rand_.eng), rand_b_xy(rand_.eng) / cfg_.tree_width + end_pos.z());
            branch_start << main_branch.getPos(split_t), main_branch.getVel(split_t), main_branch.getAcc(split_t), zero3;
            branch_end << endp, zero3, zero3, zero3;

            bvp::Piece branch = bvp_->GenFixStateMinJerkOptT(branch_start, branch_end, FIXP);
            double main_dur = branch.getDuration();
            //        double max_vel = branch.getMaxVelRate();
            double step = (cfg_.radius_max - cfg_.radius_min) / max_vel;
            for (double eval_t = 0.0; eval_t < main_dur; eval_t += cfg_.sample_t_reso) {
                double radius = step * max(0.0, max_vel - branch.getVel(eval_t).norm()) + cfg_.radius_min;
                bvp::Vec3 pos = branch.getPos(eval_t);
                for (double theta = 0.0; theta < 2 * M_PI; theta += cfg_.sample_theta_reso) {
                    pcl::PointXYZ p;
                    p.x = center.x() + pos.x() + sin(theta) * radius;
                    p.y = center.y() + pos.y() + cos(theta) * radius;
                    p.z = pos.z() + cfg_.ground_height;
                    input_map.points.push_back(p);
                }
            }
        }
    }

    void GenerateSpline(){
        pcl::PointXYZ pt_random;
        bool is_kdtree_empty = false;
        if (pcl_map_.points.size() > 0)
            kdtreeMap.setInputCloud(pcl_map_.makeShared());
        else
            is_kdtree_empty = true;
        uniform_real_distribution<double> rand_h = uniform_real_distribution<double>(cfg_.height_min, cfg_.height_max);
        uniform_int_distribution<int> rand_b_num = uniform_int_distribution<int>(0, cfg_.max_branch_num);
        for (size_t i = 0; i < cfg_.spline_num; i++) {
            double x, y, w, h, l;
            x = rand_.unif_pm_1(rand_.eng) * cfg_.map_size_x / 2;
            y = rand_.unif_pm_1(rand_.eng) * cfg_.map_size_y / 2;
            bool conti{false};
            for(size_t i = 0 ; i< cfg_.clear_area.size();i++){
                Eigen::Vector3d cur(x,y,0);
                if((cur-cfg_.clear_area[i]).norm() < cfg_.clear_radius[i]){
                    conti = true;
                    break;
                }
            }
            if(conti){
                i--;
                continue;
            }

            if (sqrt(pow(x - cfg_.init_x, 2) + pow(y - cfg_.init_y, 2)) < 0.8) {
                i--;
                continue;
            }


            pcl::PointXYZ searchPoint(x, y, (cfg_.pillar_height[0] + cfg_.pillar_height[1]) / 2.0);
            pointIdxSearch.clear();
            pointSquaredDistance.clear();

            if (is_kdtree_empty == false) {
                if (kdtreeMap.nearestKSearch(searchPoint, 1, pointIdxSearch, pointSquaredDistance) > 0) {
                    if (sqrt(pointSquaredDistance[0]) < 1.0) {
                        i--;
                        continue;
                    }

                }
            }

            bvp::Vec3 rand_pos(x,y,0);

            bvp::Vec3 start_dir(rand_.unif_pm_1(rand_.eng), rand_.unif_pm_1(rand_.eng), abs(rand_.unif_pm_1(rand_.eng)) + 1);
            bvp::Vec3 end_dir(rand_.unif_pm_1(rand_.eng), rand_.unif_pm_1(rand_.eng), abs(rand_.unif_pm_1(rand_.eng)) + 1);
            bvp::Vec3 end_pos(rand_.unif_pm_1(rand_.eng), rand_.unif_pm_1(rand_.eng), rand_h(rand_.eng));
            start_dir.normalized();
            end_dir.normalized();
            int b_num = rand_b_num(rand_.eng);
            GenerateTree(pcl_map_, rand_pos, start_dir, end_dir, end_pos, b_num);
        }

        print(fg(color::light_sea_green), " [Spline] Generate success.\n");
    }


};