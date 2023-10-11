#ifndef LOCAL_ASTAR_H
#define LOCAL_ASTAR_H

#include <ros/ros.h>
#include <Eigen/Eigen>

#include "astar.h"

#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

class LoaclAstarClass {
public:
    LoaclAstarClass(){}
    ~LoaclAstarClass(){}

    void InitMap(double resolution, Eigen::Vector3d gl_low, Eigen::Vector3d gl_upp) {
        center_  = Eigen::Vector3d::Zero();
        resolution_ = resolution;
        gl_low_ = gl_low;
        gl_upp_ = gl_upp;
        GL_SIZE_.x() = (gl_upp.x() - gl_low.x()) / resolution_;
        GL_SIZE_.y() = (gl_upp.y() - gl_low.y()) / resolution_;
        GL_SIZE_.z() = (gl_upp.z() - gl_low.z()) / resolution_;
        // std::cout << "GL_SIZE_: " << GL_SIZE_.transpose() << std::endl;
        
        grid_map_ = new uint8_t[GL_SIZE_(0)*GL_SIZE_(1)*GL_SIZE_(2)];
        memset(grid_map_, 0, GL_SIZE_(0)*GL_SIZE_(1)*GL_SIZE_(2)*sizeof(uint8_t));
    
        GridNodeMap_ = new GridNodePtr ** [GL_SIZE_(0)];
        for (int i = 0; i < GL_SIZE_(0); i++) {
            GridNodeMap_[i] = new GridNodePtr * [GL_SIZE_(1)];
            for (int j = 0; j < GL_SIZE_(1); j++) {
                GridNodeMap_[i][j] = new GridNodePtr [GL_SIZE_(2)];
                for (int k = 0; k < GL_SIZE_(2); k++) {
                    Eigen::Vector3i tmpIdx(i,j,k);
                    GridNodeMap_[i][j][k] = new GridNode(tmpIdx);
                }
            }
        }
    }
    void Reset(void) {
        for (int i = 0; i < GL_SIZE_(0) ; i++) {
            for (int j = 0; j < GL_SIZE_(1) ; j++) {
                for (int k = 0; k < GL_SIZE_(2) ; k++) {
                    GridNodeMap_[i][j][k]->id = 0;
                    GridNodeMap_[i][j][k]->cameFrom = NULL;
                    GridNodeMap_[i][j][k]->gScore = inf;
                    GridNodeMap_[i][j][k]->fScore = inf;
                }
            }
        }
    }
    void setObs(double coord_x, double coord_y, double coord_z) {
        coord_x -= center_.x();
        coord_y -= center_.y();
        coord_z -= center_.z();
        if( coord_x <  gl_low_(0) || coord_y <  gl_low_(1) || coord_z <  gl_low_(2) || 
            coord_x >= gl_upp_(0) || coord_y >= gl_upp_(1) || coord_z >= gl_upp_(2) ) {
            return;
        }

        int x = static_cast<int>((coord_x - gl_low_(0)) / resolution_);
        int y = static_cast<int>((coord_y - gl_low_(1)) / resolution_);
        int z = static_cast<int>((coord_z - gl_low_(2)) / resolution_);
        grid_map_[x*GL_SIZE_(1)*GL_SIZE_(2) + y*GL_SIZE_(2) + z] = 1;
    }
    void SetCenter(Eigen::Vector3d center) {
        center_ = center;
    }
    void setObsPcl(pcl::PointCloud<pcl::PointXYZ> &cloud, double radius = 0.2) {
        memset(grid_map_, 0, GL_SIZE_(0)*GL_SIZE_(1)*GL_SIZE_(2)*sizeof(uint8_t));
        
        for (int i = 0; i < cloud.points.size(); i++) {
            Eigen::Vector3d new_pt;
            new_pt << cloud.points[i].x, cloud.points[i].y, cloud.points[i].z;
            if (coorIsInMap(new_pt) == false) continue;
            new_pt << cloud.points[i].x, cloud.points[i].y, resolution_*2;
            bool flag = CheckPoint(new_pt);
            flag = false;
            for (double x = -radius; x <= radius; x += resolution_) {
                for (double y = -radius; y <= radius; y += resolution_) {
                    if (flag == false) {
                        for (double z = -radius; z <= radius; z += resolution_) {
                            setObs(cloud.points[i].x+x, cloud.points[i].y+y, cloud.points[i].z+z);
                        }
                    } else {
                        for (double z = 0; z <= cloud.points[i].z+radius; z += resolution_-0.02) {
                            setObs(cloud.points[i].x+x, cloud.points[i].y+y, z);
                        }
                    }
                }
            }
        }
    }
    void setObsVector(std::vector<Eigen::Vector3d> &cloud, double radius = 0.2) {
        memset(grid_map_, 0, GL_SIZE_(0)*GL_SIZE_(1)*GL_SIZE_(2)*sizeof(uint8_t));

        for (int i = 0; i < cloud.size(); i++) {
            Eigen::Vector3d new_pt;
            new_pt = cloud[i];
            if (coorIsInMap(new_pt) == false) continue;
            new_pt << cloud[i].x(), cloud[i].y(), resolution_*2;
            // bool flag = CheckPoint(new_pt);
            // flag = false;
            for (double x = -radius; x <= radius; x += resolution_) {
                for (double y = -radius; y <= radius; y += resolution_) {
                    // if (flag == false) {
                        for (double z = -radius; z <= radius; z += resolution_) {
                            setObs(cloud[i].x()+x, cloud[i].y()+y, cloud[i].z()+z);
                        }
                    // } else {
                        // for (double z = 0; z <= cloud[i].z()+radius; z += resolution_-0.02) {
                        //     setObs(cloud[i].x()+x, cloud[i].y()+y, z);
                        // }
                    // }
                }
            }
            // for (double x = -radius; x <= radius; x += radius) {
            //     for (double y = -radius; y <= radius; y += radius) {
            //         for (double z = -radius; z <= radius; z += radius) {
            //             setObs(cloud[i].x()+x, cloud[i].y()+y, cloud[i].z()+z);
            //         }
            //     }
            // }
        }
    }
    void GetOccupyPcl(pcl::PointCloud<pcl::PointXYZ> &cloud) {
        cloud.points.clear();
        for (int i = 0; i < GL_SIZE_.x(); i++) {
            for (int j = 0; j < GL_SIZE_.y(); j++) {
                for (int k = 0; k < GL_SIZE_.z(); k++) {
                    if (isOccupied(i, j, k)) {
                        pcl::PointXYZ pt;
                        pt.x = i * resolution_ + gl_low_.x() + center_.x();
                        pt.y = j * resolution_ + gl_low_.y() + center_.y();
                        pt.z = k * resolution_ + gl_low_.z() + center_.z();
                        cloud.points.push_back(pt);
                    }
                }
            }
        }
    }
    bool CheckPoint(const Eigen::Vector3d& pt) {
        if (coorIsInMap(pt) == false) return true;
        if (isOccupied(coord2gridIndex(pt))) return false;
        return true;
    }
    bool CheckStartEnd(const Eigen::Vector3d& pt) {
        Eigen::Vector3i idx = coord2gridIndex(pt);
        if (isOccupied(idx)) return false;
        else return true;
    }
    bool CheckLineObstacleFree(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) {
        Eigen::Vector3d vector = p2 - p1;
        int point_num = vector.norm() / resolution_;
        bool flag = true;
        for (int i = 1; i <= point_num; i++) {
            Eigen::Vector3d coor = p1 + vector * i / (point_num+1);
            if (coorIsInMap(coor) == false) continue;
            if (isOccupied(coord2gridIndex(coor))) {
                flag = false;
                break;
            }
        }
        return flag;
    }
    bool CheckPathFree(const std::vector<Eigen::Vector3d>& path) {
        bool flag = true;
        for (int i = 0; i < path.size(); i++) {
            if (coorIsInMap(path[i]) == false) continue;
            if (isOccupied(coord2gridIndex(path[i]))) {
                flag = false;
                break;
            }
        }
        return flag;
    }

    virtual bool SearchPath(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);
    void GetPath(std::vector<Eigen::Vector3d> &path);
    void SimplifyPath(const std::vector<Eigen::Vector3d>& astar_path, 
                            std::vector<Eigen::Vector3d>& waypoint);
    void FloydHandle(const std::vector<Eigen::Vector3d>& astar_path, 
                           std::vector<Eigen::Vector3d>& waypoint);
    
    double resolution_;
    Eigen::Vector3i GL_SIZE_;
    Eigen::Vector3d gl_low_, gl_upp_;
    Eigen::Vector3d center_;

private:
    inline bool coorIsInMap(const Eigen::Vector3d & pt) {
        Eigen::Vector3d error = pt - center_;
        if (error.x() > gl_low_.x() && error.x() < gl_upp_.x() && 
            error.y() > gl_low_.y() && error.y() < gl_upp_.y() && 
            error.z() > gl_low_.z() && error.z() < gl_upp_.z()) {
            return true;
        }
        return false;
    }
    inline Eigen::Vector3d gridIndex2coord(const Eigen::Vector3i & index) {
        Eigen::Vector3d pt;
        pt(0) = ((double)index(0) + 0.5) * resolution_ + gl_low_(0);
        pt(1) = ((double)index(1) + 0.5) * resolution_ + gl_low_(1);
        pt(2) = ((double)index(2) + 0.5) * resolution_ + gl_low_(2);
        pt += center_;
        return pt;
    }
    inline Eigen::Vector3i coord2gridIndex(const Eigen::Vector3d & pt) {
        Eigen::Vector3d c_pt = pt - center_;
        Eigen::Vector3i index;
        index(0) = std::min(std::max(int((c_pt(0) - gl_low_(0)) / resolution_), 0), GL_SIZE_(0)-1);
        index(1) = std::min(std::max(int((c_pt(1) - gl_low_(1)) / resolution_), 0), GL_SIZE_(1)-1);
        index(2) = std::min(std::max(int((c_pt(2) - gl_low_(2)) / resolution_), 0), GL_SIZE_(2)-1);
        return index;
    }
    inline bool isOccupied(int idx_x, int idx_y, int idx_z) {
        return  (idx_x >= 0 && idx_x < GL_SIZE_(0) && idx_y >= 0 && idx_y < GL_SIZE_(1) && idx_z >= 0 && idx_z < GL_SIZE_(2) && 
                (grid_map_[idx_x*GL_SIZE_(1)*GL_SIZE_(2) + idx_y*GL_SIZE_(2) + idx_z] == 1));
    }
    inline bool isOccupied(const Eigen::Vector3i index) {
        return isOccupied(index(0), index(1), index(2));
    }
    inline bool isFree(int idx_x, int idx_y, int idx_z) {
        return  (idx_x >= 0 && idx_x < GL_SIZE_(0) && idx_y >= 0 && idx_y < GL_SIZE_(1) && idx_z >= 0 && idx_z < GL_SIZE_(2) && 
                (grid_map_[idx_x*GL_SIZE_(1)*GL_SIZE_(2) + idx_y*GL_SIZE_(2) + idx_z] < 1));
    }
    inline bool isFree(const Eigen::Vector3i index) {
        return isFree(index(0), index(1), index(2));
    }

    double CalcHeu(GridNodePtr node1, GridNodePtr node2);
    void GetNeighbors(GridNodePtr cur, std::vector<GridNodePtr> &neighbors, std::vector<double> &costs);

    uint8_t *grid_map_;
    GridNodePtr ***GridNodeMap_;
    Eigen::Vector3d start_pt_, end_pt_;
    GridNodePtr terminatePtr_;
    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, Compare> open_lists_;
};

#endif
