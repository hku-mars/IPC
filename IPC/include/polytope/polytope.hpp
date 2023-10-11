#ifndef _GEOMETYR_UTILS_POLYTOPE_
#define _GEOMETYR_UTILS_POLYTOPE_

#include "Eigen/Dense"
#include "vector"
#include "list"
#include "common_type_name.hpp"
#include "ellipsoid.hpp"
#include "line_segment.hpp"
// #include "geo_utils.hpp"

class Polytope {
private:
    bool undefined{true};
    bool is_known_free{false};
    MatD4f planes;
    bool have_seed_line{false};
public:
    double overlap_depth_with_last_one{0};
    Vec3f interior_pt_with_last_one;
    Ellipsoid ellipsoid_;
    std::pair<Vec3f, Vec3f> seed_line;
    double robot_r;

    Polytope() {}

    Polytope(MatD4f _planes) {
        planes = _planes;
        undefined = false;
    }

    bool HaveSeedLine() {
        return have_seed_line;
    }

    void SetSeedLine(const std::pair<Vec3f, Vec3f> &_seed_line, double r = 0) {
        robot_r = r;
        seed_line = _seed_line;
        have_seed_line = true;
    }

    int SurfNum() {
        if (undefined) {
            return 0;
        }
        return planes.rows();
    }

    Polytope CrossWith(Polytope &b) {
        MatD4f curIH;
        curIH.resize(this->SurfNum() + b.SurfNum(), 4);
        curIH << this->planes, b.GetPlanes();
        Polytope out;
        out.SetPlanes(curIH);
        return out;
    }

    Vec3f CrossCenter(Polytope &b) {
        MatD4f curIH;
        curIH.resize(this->SurfNum() + b.SurfNum(), 4);
        curIH << this->planes, b.GetPlanes();
        Mat3Df curIV; // 走廊的顶点
        if (!geo_utils::enumerateVs(curIH, curIV)) {
            print(fg(color::gold), " -- [processCorridor] Failed to get Overlap enumerateVs .\n");
            return Vec3f(-999, -999, -999);
        }
        double x = (curIV.row(0).maxCoeff() + curIV.row(0).minCoeff()) * 0.5;
        double y = (curIV.row(1).maxCoeff() + curIV.row(1).minCoeff()) * 0.5;
        double z = (curIV.row(2).maxCoeff() + curIV.row(2).minCoeff()) * 0.5;
        return Vec3f(x, y, z);
    }

    bool HaveOverlapWith(Polytope cmp, double eps = 1e-6) {
        return geo_utils::overlap(this->planes, cmp.GetPlanes(), eps);
    }

    MatD4f GetPlanes() {
        return planes;
    }

    void Reset() {
        undefined = true;
        is_known_free = false;
    }

    bool IsKnownFree() {
        if (undefined) {
            return false;
        }
        return is_known_free;
    }

    void SetKnownFree(bool is_free) {
        is_known_free = is_free;
    }

    void SetPlanes(MatD4f _planes) {
        planes = _planes;
        undefined = false;
    }

    void SetEllipsoid(const Ellipsoid &ellip) {
        ellipsoid_ = ellip;
    }

    bool PointIsInside(const Vec3f &pt, const double margin = 0.01) const {
        if (undefined) {
            return false;
        }
        Eigen::Vector4d pt_e;
        pt_e.head(3) = pt;
        pt_e(3) = 1;
        if ((planes * pt_e).maxCoeff() > margin) {
            return false;
        }
        return true;
    }

    /// Sort a set of points in a plane.
    static vec_E<Vec3f> SortPtsInClockWise(vec_E<Vec3f> &pts, Vec3f normal) {
        Vec3f center(0,0,0);
        for (auto pt: pts) {
            center += pt;
        }
        center = center / pts.size();
        Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(Vec3f(0, 0, 1), normal);
        center = q.matrix() * center;
        vector<std::pair<double, Vec3f>> pts_valued;
        pts_valued.resize(pts.size());
        Vec3f temp_p;
        for (int i = 0; i < pts.size(); i++) {
            temp_p = q.matrix() * pts[i];
            double theta = atan2(temp_p(1) - center(1), temp_p(0) - center(0));
            pts_valued[i] = make_pair(theta, pts[i]);
        }

        std::sort(
                pts_valued.begin(), pts_valued.end(),
                [](const std::pair<double, Vec3f> &i,
                   const std::pair<double, Vec3f> &j) { return i.first < j.first; });
        vec_E<Vec3f> pts_sorted(pts_valued.size());
        for (size_t i = 0; i < pts_valued.size(); i++)
            pts_sorted[i] = pts_valued[i].second;
        return pts_sorted;
    }

    void Visualize(const ros::Publisher &pub,
                   string ns = "mesh",
                   bool flag = false,
                   Color surf_color = Color::SteelBlue(),
                   Color edge_color = Color::Black(),
                   Color vertex_color = Color::Red(),
                   double alpha = 0.3) {
        // Due to the fact that H-representation cannot be directly visualized
        // We first conduct vertex enumeration of them, then apply quickhull
        // to obtain triangle meshs of polyhedra
        Eigen::Matrix3Xd mesh(3, 0), curTris(3, 0), oldTris(3, 0);

        if (!flag) ellipsoid_.Visualize(pub);
        LineSegment l(seed_line.first, seed_line.second);
        if (!flag) l.Visualize(pub);

        oldTris = mesh;
        Eigen::Matrix<double, 3, -1, Eigen::ColMajor> vPoly;
        MatD4f planes = GetPlanes();
        geo_utils::enumerateVs(planes, vPoly);

        quickhull::QuickHull<double> tinyQH;
        const auto polyHull = tinyQH.getConvexHull(vPoly.data(), vPoly.cols(), false, true);
        const auto &idxBuffer = polyHull.getIndexBuffer();
        int hNum = idxBuffer.size() / 3;

        curTris.resize(3, hNum * 3);
        for (int i = 0; i < hNum * 3; i++) {
            curTris.col(i) = vPoly.col(idxBuffer[i]);
        }
        mesh.resize(3, oldTris.cols() + curTris.cols());
        mesh.leftCols(oldTris.cols()) = oldTris;
        mesh.rightCols(curTris.cols()) = curTris;


        // RVIZ support tris for visualization
        visualization_msgs::MarkerArray mkr_arr;
        visualization_msgs::Marker meshMarker, edgeMarker, vertexMarker;
        static int mesh_id = 0;
        static int edge_id = 0;
        static int vertex_id = 0;
        if (flag) {
            mesh_id = 0;
            edge_id = 0;
            vertex_id = 0;
        }

        vertexMarker.id = vertex_id++;
        vertexMarker.header.stamp = ros::Time::now();
        vertexMarker.header.frame_id = "world";
        vertexMarker.pose.orientation.w = 1.00;
        if (flag) vertexMarker.action = visualization_msgs::Marker::DELETEALL;
        else {
            vertexMarker.action = visualization_msgs::Marker::ADD;
            vertexMarker.type = visualization_msgs::Marker::SPHERE_LIST;
            vertexMarker.ns = ns + " vertex";
            vertexMarker.color = vertex_color;
            vertexMarker.color.a = 1.0;
            vertexMarker.scale.x = 0.2;
            vertexMarker.scale.y = 0.2;
            vertexMarker.scale.z = 0.2;
        
            geometry_msgs::Point point;
            for(int i = 0 ; i < vPoly.cols() ; i++){
                point.x = vPoly(0,i);
                point.y = vPoly(1,i);
                point.z = vPoly(2,i);
                vertexMarker.points.push_back(point);
            }
        }

        meshMarker.id = mesh_id++;
        meshMarker.header.stamp = ros::Time::now();
        meshMarker.header.frame_id = "world";
        meshMarker.pose.orientation.w = 1.00;
        if (flag) meshMarker.action = visualization_msgs::Marker::DELETEALL;
        else {
            meshMarker.action = visualization_msgs::Marker::ADD;
            meshMarker.type = visualization_msgs::Marker::TRIANGLE_LIST;
            meshMarker.ns = ns + " mesh";
            meshMarker.color = surf_color;
            meshMarker.color.a = alpha;
            meshMarker.scale.x = 1.0;
            meshMarker.scale.y = 1.0;
            meshMarker.scale.z = 1.0;

            int ptnum = mesh.cols();
            geometry_msgs::Point point;
            for (int i = 0; i < ptnum; i++) {
                point.x = mesh(0, i);
                point.y = mesh(1, i);
                point.z = mesh(2, i);
                meshMarker.points.push_back(point);
            }
        }

        edgeMarker = meshMarker;
        edgeMarker.id = edge_id++;
        
        if (flag) edgeMarker.action = visualization_msgs::Marker::DELETEALL;
        else {
            edgeMarker.type = visualization_msgs::Marker::LINE_LIST;
            edgeMarker.ns = ns + " edge";
            edgeMarker.color = edge_color;
            edgeMarker.color.a = 1.00;
            edgeMarker.scale.x = 0.05;
        }

        if (!flag) {
            /// 接下来遍历选择定点属于哪些面
            vec_E<Vec3f> edges;
            // print(fg(color::orange), "PlaneNum = {}.\n", planes.rows());
            for (int i = 0; i < planes.rows(); i++) {
                edges.clear();
                Eigen::VectorXd temp = (planes.row(i).head(3) * vPoly);
                Eigen::VectorXd d(temp.size());
                d.setConstant(planes(i, 3));
                temp = temp + d;
                // cout<<temp.transpose()<<endl;
                for (int j = 0; j < temp.size(); j++) {
                    if (abs(temp(j)) < epsilon_) {
                        edges.push_back(vPoly.col(j));
                    }
                }
                if (edges.size() < 2) {
                    continue;
                }
                // print(fg(color::orange), "edge num = {}, vertex num = {}.\n",edges.size(), vPoly.cols());
                /// 对每个面的定点逆时针排序
                vec_E<Vec3f> pts_sorted = SortPtsInClockWise(edges, planes.row(i).head(3));
                int pts_num = pts_sorted.size();
                geometry_msgs::Point point;
                for (int k = 0; k < pts_num - 1; k++) {
                    point.x = pts_sorted[k].x();
                    point.y = pts_sorted[k].y();
                    point.z = pts_sorted[k].z();
                    edgeMarker.points.push_back(point);
                    point.x = pts_sorted[k + 1].x();
                    point.y = pts_sorted[k + 1].y();
                    point.z = pts_sorted[k + 1].z();
                    edgeMarker.points.push_back(point);
                }
                point.x = pts_sorted[0].x();
                point.y = pts_sorted[0].y();
                point.z = pts_sorted[0].z();
                edgeMarker.points.push_back(point);
                point.x = pts_sorted[pts_num - 1].x();
                point.y = pts_sorted[pts_num - 1].y();
                point.z = pts_sorted[pts_num - 1].z();
                edgeMarker.points.push_back(point);
            }
        }
        mkr_arr.markers.push_back(meshMarker);
        mkr_arr.markers.push_back(edgeMarker);
        mkr_arr.markers.push_back(vertexMarker);
        pub.publish(mkr_arr);
    }

};

typedef vector<Polytope> PolytopeVec;

#endif