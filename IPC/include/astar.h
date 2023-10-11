#ifndef ASTAR_H
#define ASTAR_H

#include <ros/ros.h>
#include <Eigen/Eigen>
#include <queue>

constexpr double inf = 1 >> 20;
struct GridNode;
typedef GridNode* GridNodePtr;

struct GridNode {
    GridNode(){
        GridNode(Eigen::Vector3i(0,0,0));
    }
    GridNode(Eigen::Vector3i _index){  
		id = 0;
		index = _index;
        dir = Eigen::Vector3i::Zero();
		gScore = inf;
		fScore = inf;
		cameFrom = NULL;
    }
    ~GridNode(){}

    int8_t id; // 1--> open set, -1 --> closed set
    Eigen::Vector3i index;
    Eigen::Vector3i dir;
    double gScore;
    double fScore;
    GridNodePtr cameFrom;
};

struct Compare {
    bool operator()(const GridNodePtr a, const GridNodePtr b){
        return a->fScore > b->fScore;
    }
};

class AstarClass {
public:
    AstarClass(){}
    ~AstarClass(){}

    void InitMap(double resolution, Eigen::Vector3d gl_low, Eigen::Vector3d gl_upp, Eigen::Vector3i GL_SIZE) {
        resolution_ = resolution;
        gl_low_ = gl_low;
        gl_upp_ = gl_upp;
        GL_SIZE_ = GL_SIZE;

        grid_map_ = new uint8_t[GL_SIZE_(0)*GL_SIZE_(1)*GL_SIZE_(2)];
        memset(grid_map_, 0, GL_SIZE_(0)*GL_SIZE_(1)*GL_SIZE_(2)*sizeof(uint8_t));

        GridNodeMap_ = new GridNodePtr ** [GL_SIZE_(0)];    //三级指针
        for (int i = 0; i < GL_SIZE_(0); i++) {
            GridNodeMap_[i] = new GridNodePtr * [GL_SIZE_(1)];
            for (int j = 0; j < GL_SIZE_(1); j++) {
                GridNodeMap_[i][j] = new GridNodePtr [GL_SIZE_(2)];
                for (int k = 0; k < GL_SIZE_(2); k++) {
                    Eigen::Vector3i tmpIdx(i,j,k);
                    GridNodeMap_[i][j][k] = new GridNode(tmpIdx); // 以栅格点坐标初始化一个栅格点
                }
            }
        }
    }
    void setObs(const double coord_x, const double coord_y, const double coord_z) {
        if( coord_x <  gl_low_(0) || coord_y <  gl_low_(1) || coord_z <  gl_low_(2) || 
            coord_x >= gl_upp_(0) || coord_y >= gl_upp_(1) || coord_z >= gl_upp_(2) ) {
            return;
        }
        
        int x = static_cast<int>((coord_x - gl_low_(0)) / resolution_);
        int y = static_cast<int>((coord_y - gl_low_(1)) / resolution_);
        int z = static_cast<int>((coord_z - gl_low_(2)) / resolution_);
        grid_map_[x*GL_SIZE_(1)*GL_SIZE_(2) + y*GL_SIZE_(2) + z] = 1;
    }
    Eigen::Vector3d coordRounding(Eigen::Vector3d &coord) {
        return gridIndex2coord(coord2gridIndex(coord));
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
    void GetVisitNodes(std::vector<Eigen::Vector3d> &nodes) {
        nodes.clear();
        for (int i = 0; i < GL_SIZE_(0) ; i++) {
            for (int j = 0; j < GL_SIZE_(1) ; j++) {
                for (int k = 0; k < GL_SIZE_(2) ; k++) {
                    if (GridNodeMap_[i][j][k]->id == -1) {
                        nodes.push_back(gridIndex2coord(GridNodeMap_[i][j][k]->index));
                    }
                }
            }
        }
    }
    bool CheckObstacle(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) {
        Eigen::Vector3d vector = p2 - p1;
        int point_num = vector.norm() / resolution_;
        bool flag = true;
        for (int i = 0; i <= point_num; i++) {
            Eigen::Vector3d coor = p1 + vector * i / (point_num+1);
            if (isOccupied(coord2gridIndex(coor))) {
                flag = false;
                break;
            }
        }
        return flag;
    }

    virtual bool SearchPath(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);
    double CalcHeu(GridNodePtr node1, GridNodePtr node2);
    void GetPath(std::vector<Eigen::Vector3d> &path);
    void SimplifyPath(const std::vector<Eigen::Vector3d>& astar_path, 
                            std::vector<Eigen::Vector3d>& waypoint);
    void FloydHandle(const std::vector<Eigen::Vector3d>& astar_path, 
                           std::vector<Eigen::Vector3d>& waypoint);
    void InsertPoints(std::vector<Eigen::Vector3d>& waypoint, double inter_dis);

private:
    void GetNeighbors(GridNodePtr cur, std::vector<GridNodePtr> &neighbors, std::vector<double> &costs);

protected:
    inline bool isOccupied(const int &idx_x, const int &idx_y, const int &idx_z) {
        return  (idx_x >= 0 && idx_x < GL_SIZE_(0) && idx_y >= 0 && idx_y < GL_SIZE_(1) && idx_z >= 0 && idx_z < GL_SIZE_(2) && 
                (grid_map_[idx_x*GL_SIZE_(1)*GL_SIZE_(2) + idx_y*GL_SIZE_(2) + idx_z] == 1));
    }
    inline bool isOccupied(const Eigen::Vector3i &index) {
        return isOccupied(index(0), index(1), index(2));
    }
    inline bool isFree(const int &idx_x, const int &idx_y, const int &idx_z) {
        return  (idx_x >= 0 && idx_x < GL_SIZE_(0) && idx_y >= 0 && idx_y < GL_SIZE_(1) && idx_z >= 0 && idx_z < GL_SIZE_(2) && 
                (grid_map_[idx_x*GL_SIZE_(1)*GL_SIZE_(2) + idx_y*GL_SIZE_(2) + idx_z] < 1));
    }
    inline bool isFree(const Eigen::Vector3i &index) {
        return isFree(index(0), index(1), index(2));
    }

    inline Eigen::Vector3d gridIndex2coord(const Eigen::Vector3i & index) {
        Eigen::Vector3d pt;
        pt(0) = ((double)index(0) + 0.5) * resolution_ + gl_low_(0);
        pt(1) = ((double)index(1) + 0.5) * resolution_ + gl_low_(1);
        pt(2) = ((double)index(2) + 0.5) * resolution_ + gl_low_(2);
        return pt;
    }
    inline Eigen::Vector3i coord2gridIndex(const Eigen::Vector3d & pt) {
        Eigen::Vector3i index;
        index(0) = std::min(std::max(int((pt(0) - gl_low_(0)) / resolution_), 0), GL_SIZE_(0)-1);
        index(1) = std::min(std::max(int((pt(1) - gl_low_(1)) / resolution_), 0), GL_SIZE_(1)-1);
        index(2) = std::min(std::max(int((pt(2) - gl_low_(2)) / resolution_), 0), GL_SIZE_(2)-1);
        return index;
    }

    uint8_t *grid_map_;
    GridNodePtr ***GridNodeMap_;

    double resolution_;
    Eigen::Vector3i GL_SIZE_;
    Eigen::Vector3d gl_low_, gl_upp_;
    Eigen::Vector3d start_pt_, end_pt_;

    GridNodePtr terminatePtr_;
    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, Compare> open_lists_;
};


#endif
