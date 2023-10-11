//
// Created by yunfan on 2022/3/26.
//

#ifndef MISSION_PLANNER_CONFIG
#define MISSION_PLANNER_CONFIG

#include "ros/ros.h"
#include "vector"
#include "fmt/color.h"
#include "string"
#include "ros/package.h"
#include "Eigen/Core"

using namespace std;
using namespace fmt;

namespace racing_map_generator {
    class RacingMapConfig {
    public:
        // Bool Params

        vector<Eigen::Vector4d> pos_yaw_pair;
        double size,pub_freq,resolution;

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

        template<class T>
        bool LoadParam(string param_name, vector<T> &param_value, vector<T> default_value) {
            if (nh_.getParam(param_name, param_value)) {
                printf("\033[0;32m Load param %s succes: \033[0;0m", param_name.c_str());
                for (int i = 0; i < param_value.size(); i++) {
                    cout << param_value[i] << " ";
                }
                cout << endl;
                return true;
            } else {
                printf("\033[0;33m Load param %s failed, use default value: \033[0;0m", param_name.c_str());
                param_value = default_value;
                for (int i = 0; i < param_value.size(); i++) {
                    cout << param_value[i] << " ";
                }
                cout << endl;
                return false;
            }
        }

        ros::NodeHandle nh_;

        RacingMapConfig() {};

        RacingMapConfig(const ros::NodeHandle &nh_priv) {
            nh_ = nh_priv;
            vector<double> x, y, z, yaw;
            LoadParam("racing_map/waypoints/x", x, x);
            LoadParam("racing_map/waypoints/y", y, y);
            LoadParam("racing_map/waypoints/z", z, z);
            LoadParam("racing_map/waypoints/yaw", yaw, yaw);
            LoadParam("racing_map/size", size, size);
            LoadParam("racing_map/resolution", resolution, resolution);
            LoadParam("racing_map/pub_freq", pub_freq, pub_freq);

            for (int i = 0; i < x.size(); i++) {
                double yaw_rad = yaw[i] / 180.0 * M_PI;
                pos_yaw_pair.push_back(Eigen::Vector4d(x[i], y[i], z[i], yaw_rad));
            }

            print(fg(color::lime_green), " -- [MISSION] Load {} waypoints success.\n", pos_yaw_pair.size());
            for (int i = 0; i < pos_yaw_pair.size(); i++) {
                print(fg(color::light_sky_blue), "\t Waypoint {0} at ({1:.2f}, {2:.2f}, {3:.2f}, {4:.2f})\n", i,
                      pos_yaw_pair[i][0], pos_yaw_pair[i][1], pos_yaw_pair[i][2], pos_yaw_pair[i][3]);
            }

        }

    };
}
#endif //PLANNER_CONFIG_HPP
