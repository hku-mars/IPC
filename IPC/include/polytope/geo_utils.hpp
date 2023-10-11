/*
    MIT License

    Copyright (c) 2021 Zhepei Wang (wangzhepei@live.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef GEO_UTILS_HPP
#define GEO_UTILS_HPP

#include "quickhull.hpp"
#include "sdlp.hpp"

#include <Eigen/Eigen>

#include <cfloat>
#include <cstdint>
#include <set>
#include <chrono>

namespace geo_utils {

using namespace Eigen;

    inline Eigen::Matrix3d RotationFromVec3(const Eigen::Vector3d &v) {
        // zero roll
        Eigen::Vector3d rpy(0, std::atan2(-v(2), v.topRows<2>().norm()),
                            std::atan2(v(1), v(0)));
        Eigen::Quaterniond qx(cos(rpy(0) / 2), sin(rpy(0) / 2), 0, 0);
        Eigen::Quaterniond qy(cos(rpy(1) / 2), 0, sin(rpy(1) / 2), 0);
        Eigen::Quaterniond qz(cos(rpy(2) / 2), 0, 0, sin(rpy(2) / 2));
        return Eigen::Matrix3d(qz * qy * qx);
    }

    // 通过三点获得一个平面
    void FromPointsToPlane(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, const Eigen::Vector3d &p3,
                           Eigen::Vector4d &hPoly) {
        // Each row of hPoly is defined by h0, h1, h2, h3 as
        // h0*x + h1*y + h2*z + h3 <= 0
        hPoly(0) = ((p2.y() - p1.y()) * (p3.z() - p1.z()) - (p2.z() - p1.z()) * (p3.y() - p1.y()));
        hPoly(1) = ((p2.z() - p1.z()) * (p3.x() - p1.x()) - (p2.x() - p1.x()) * (p3.z() - p1.z()));
        hPoly(2) = ((p2.x() - p1.x()) * (p3.y() - p1.y()) - (p2.y() - p1.y()) * (p3.x() - p1.x()));
        hPoly(3) = (0 - (hPoly(0) * p1.x() + hPoly(1) * p1.y() + hPoly(2) * p1.z()));
    }

    void GetFovCheckPlane(const Eigen::Matrix3d R, const Eigen::Vector3d t, Eigen::MatrixX4d &fov_planes,
                          std::vector<Eigen::Matrix3d> &fov_pts) {
        // 只使用上下两个切面来约束。
        fov_planes.resize(2, 4);
        fov_pts.clear();
        static const double sqrt2 = sqrt(2);
        static Eigen::Vector3d inner_pt(1, 0, 0);
        static constexpr double rad60 = (35.0 / 180.0 * M_PI);
        static constexpr double radm5 = (-5.0 / 180.0 * M_PI);
        static double z60 = sin(rad60) * 5;
        static double r60 = cos(rad60) * 5;
        static double zm5 = sin(radm5) * 5;
        static double rm5 = cos(radm5) * 5;

        Eigen::Matrix3Xd four_pts(3, 4);
        four_pts.col(0) = Eigen::Vector3d(r60, -3, z60);
        four_pts.col(1) = Eigen::Vector3d(r60, 3, z60);
        four_pts.col(2) = Eigen::Vector3d(rm5, -3, zm5);
        four_pts.col(3) = Eigen::Vector3d(rm5, 3, zm5);
        four_pts = (R * four_pts).colwise() + t;
        Eigen::Vector3d fov_inner_pt = R * inner_pt + t;
        Eigen::Vector4d temp;
        FromPointsToPlane(four_pts.col(0), four_pts.col(1), t, temp);
        if (temp.head(3).dot(fov_inner_pt) + temp(3) > 0) {
            temp = -temp;
        }

        Eigen::Matrix3d temp_p;
        temp_p << four_pts.col(0), four_pts.col(1), t;
        fov_pts.push_back(temp_p);
        temp_p << four_pts.col(2), four_pts.col(3), t;
        fov_pts.push_back(temp_p);
        fov_planes.row(0) = temp;
        FromPointsToPlane(four_pts.col(2), four_pts.col(3), t, temp);
        if (temp.head(3).dot(fov_inner_pt) + temp(3) > 0) {
            temp = -temp;
        }

        fov_planes.row(1) = temp;
    }

    void GetFovPlanes(const Eigen::Matrix3d R, const Eigen::Vector3d t, Eigen::MatrixX4d &fov_planes,
                      std::vector<Eigen::Matrix3d> &fov_pts) {
        fov_planes.resize(10, 4);
        fov_pts.clear();
        // 1) 获得前三个面
        static const double sqrt2 = sqrt(2);
        static Eigen::Vector3d center(0, 0, 0);
        static Eigen::Vector3d inner_pt(0, 0, -1);
        static Eigen::Vector3d inner_pt2(0, 0, 1);
        static constexpr double rad60 = (35.0 / 180.0 * M_PI);
        static double z = sin(rad60) * 5;
        static double r = cos(rad60) * 5;
        Eigen::Matrix3Xd ten_pts(3, 10);
        ten_pts.col(0) = (Eigen::Vector3d(0, r, z));
        ten_pts.col(1) = (Eigen::Vector3d(r / sqrt2, r / sqrt2, z));
        ten_pts.col(2) = (Eigen::Vector3d(r, 0, z));
        ten_pts.col(3) = (Eigen::Vector3d(r / sqrt2, -r / sqrt2, z));
        ten_pts.col(4) = (Eigen::Vector3d(0, -r, z));


        // 2) 获得下三个面
        constexpr double rad85 = (55.0 / 180.0 * M_PI);
        static double z85 = -cos(rad85) * 5;
        static double r85 = sin(rad85) * 5;
        ten_pts.col(5) = (Eigen::Vector3d(0, r85, z85));
        ten_pts.col(6) = (Eigen::Vector3d(r85 / sqrt2, r85 / sqrt2, z85));
        ten_pts.col(7) = (Eigen::Vector3d(r85, 0, z85));
        ten_pts.col(8) = (Eigen::Vector3d(r85 / sqrt2, -r85 / sqrt2, z85));
        ten_pts.col(9) = (Eigen::Vector3d(0, -r85, z85));

        // 3) 旋转所有点
        ten_pts = (R * ten_pts).colwise() + t;
        Eigen::Vector3d fov_inner_pt = R * inner_pt + t;
        Eigen::Vector3d fov_inner_pt2 = R * inner_pt2 + t;
        Eigen::Matrix3d fov_pt;
        for (int i = 0; i < 4; i++) {
            Eigen::Vector4d temp;
            FromPointsToPlane(ten_pts.col(i), ten_pts.col(i + 1), t, temp);
            fov_pt << ten_pts.col(i), ten_pts.col(i + 1), t;
            fov_pts.push_back(fov_pt);
            if (temp.head(3).dot(fov_inner_pt) + temp(3) > 0) {
                temp = -temp;
            }
            fov_planes.row(i) = temp;
        }
        for (int i = 5; i < 9; i++) {
            Eigen::Vector4d temp;
            FromPointsToPlane(ten_pts.col(i), ten_pts.col(i + 1), t, temp);
            fov_pt << ten_pts.col(i), ten_pts.col(i + 1), t;
            fov_pts.push_back(fov_pt);
            if (temp.head(3).dot(fov_inner_pt2) + temp(3) > 0) {
                temp = -temp;
            }
            fov_planes.row(i) = temp;
        }

        // 4）转为平面方程输出
    }


    double findInteriorDist(const Eigen::MatrixX4d &hPoly,
                            Eigen::Vector3d &interior) {
        const int m = hPoly.rows();

        Eigen::MatrixX4d A(m, 4);
        Eigen::VectorXd b(m);
        Eigen::Vector4d c, x;
        const Eigen::ArrayXd hNorm = hPoly.leftCols<3>().rowwise().norm();
        A.leftCols<3>() = hPoly.leftCols<3>().array().colwise() / hNorm;
        A.rightCols<1>().setConstant(1.0);
        b = -hPoly.rightCols<1>().array() / hNorm;
        c.setZero();
        c(3) = -1.0;

        const double minmaxsd = sdlp::linprog<4>(c, A, b, x);
        interior = x.head<3>();
        return -minmaxsd;
    }

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    inline bool findInterior(const Eigen::MatrixX4d &hPoly,
                             Eigen::Vector3d &interior) {
        const int m = hPoly.rows();

        Eigen::MatrixX4d A(m, 4);
        Eigen::VectorXd b(m);
        Eigen::Vector4d c, x;
        const Eigen::ArrayXd hNorm = hPoly.leftCols<3>().rowwise().norm();
        A.leftCols<3>() = hPoly.leftCols<3>().array().colwise() / hNorm;
        A.rightCols<1>().setConstant(1.0);
        b = -hPoly.rightCols<1>().array() / hNorm;
        c.setZero();
        c(3) = -1.0;

        const double minmaxsd = sdlp::linprog<4>(c, A, b, x);
        interior = x.head<3>();

        return minmaxsd < 0.0 && !std::isinf(minmaxsd);
    }

    inline bool overlap(const Eigen::MatrixX4d &hPoly0,
                        const Eigen::MatrixX4d &hPoly1,
                        const double eps = 1.0e-6) {
        const int m = hPoly0.rows();
        const int n = hPoly1.rows();
        Eigen::MatrixX4d A(m + n, 4);
        Eigen::Vector4d c, x;
        Eigen::VectorXd b(m + n);
        A.leftCols<3>().topRows(m) = hPoly0.leftCols<3>();
        A.leftCols<3>().bottomRows(n) = hPoly1.leftCols<3>();
        A.rightCols<1>().setConstant(1.0);
        b.topRows(m) = -hPoly0.rightCols<1>();
        b.bottomRows(n) = -hPoly1.rightCols<1>();
        c.setZero();
        c(3) = -1.0;

        const double minmaxsd = sdlp::linprog<4>(c, A, b, x);

        return minmaxsd < -eps && !std::isinf(minmaxsd);
    }

    struct filterLess {
        inline bool operator()(const Eigen::Vector3d &l,
                               const Eigen::Vector3d &r) {
            return l(0) < r(0) ||
                   (l(0) == r(0) &&
                    (l(1) < r(1) ||
                     (l(1) == r(1) &&
                      l(2) < r(2))));
        }
    };

    inline void filterVs(const Eigen::Matrix3Xd &rV,
                         const double &epsilon,
                         Eigen::Matrix3Xd &fV) {
        const double mag = std::max(fabs(rV.maxCoeff()), fabs(rV.minCoeff()));
        const double res = mag * std::max(fabs(epsilon) / mag, DBL_EPSILON);
        std::set<Eigen::Vector3d, filterLess> filter;
        fV = rV;
        int offset = 0;
        Eigen::Vector3d quanti;
        for (int i = 0; i < rV.cols(); i++) {
            quanti = (rV.col(i) / res).array().round();
            if (filter.find(quanti) == filter.end()) {
                filter.insert(quanti);
                fV.col(offset) = rV.col(i);
                offset++;
            }
        }
        fV = fV.leftCols(offset).eval();
        return;
    }

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    // proposed epsilon is 1.0e-6
    inline void enumerateVs(const Eigen::MatrixX4d &hPoly,
                            const Eigen::Vector3d &inner,
                            Eigen::Matrix3Xd &vPoly,
                            const double epsilon = 1.0e-6) {
        const Eigen::VectorXd b = -hPoly.rightCols<1>() - hPoly.leftCols<3>() * inner;
        const Eigen::Matrix<double, 3, -1, Eigen::ColMajor> A =
                (hPoly.leftCols<3>().array().colwise() / b.array()).transpose();

        quickhull::QuickHull<double> qh;
        const double qhullEps = std::min(epsilon, quickhull::defaultEps<double>());
        // CCW is false because the normal in quickhull towards interior
        const auto cvxHull = qh.getConvexHull(A.data(), A.cols(), false, true, qhullEps);
        const auto &idBuffer = cvxHull.getIndexBuffer();
        const int hNum = idBuffer.size() / 3;
        Eigen::Matrix3Xd rV(3, hNum);
        Eigen::Vector3d normal, point, edge0, edge1;
        for (int i = 0; i < hNum; i++) {
            point = A.col(idBuffer[3 * i + 1]);
            edge0 = point - A.col(idBuffer[3 * i]);
            edge1 = A.col(idBuffer[3 * i + 2]) - point;
            normal = edge0.cross(edge1); //cross in CW gives an outter normal
            rV.col(i) = normal / normal.dot(point);
        }
        filterVs(rV, epsilon, vPoly);
        vPoly = (vPoly.array().colwise() + inner.array()).eval();
        return;
    }

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    // proposed epsilon is 1.0e-6
    inline bool enumerateVs(const Eigen::MatrixX4d &hPoly,
                            Eigen::Matrix3Xd &vPoly,
                            const double epsilon = 1.0e-6) {
        Eigen::Vector3d inner;
        if (findInterior(hPoly, inner)) {
            enumerateVs(hPoly, inner, vPoly, epsilon);
            return true;
        } else {
            return false;
        }
    }

} // namespace geo_utils

namespace uav_utils {

    template<typename Scalar_t>
    Scalar_t toRad(const Scalar_t &x) {
        return x / 180.0 * M_PI;
    }

    template<typename Scalar_t>
    Scalar_t toDeg(const Scalar_t &x) {
        return x * 180.0 / M_PI;
    }

    template<typename Scalar_t>
    Eigen::Matrix<Scalar_t, 3, 3> rotx(Scalar_t t) {
        Eigen::Matrix<Scalar_t, 3, 3> R;
        R(0, 0) = 1.0;
        R(0, 1) = 0.0;
        R(0, 2) = 0.0;
        R(1, 0) = 0.0;
        R(1, 1) = std::cos(t);
        R(1, 2) = -std::sin(t);
        R(2, 0) = 0.0;
        R(2, 1) = std::sin(t);
        R(2, 2) = std::cos(t);

        return R;
    }

    template<typename Scalar_t>
    Eigen::Matrix<Scalar_t, 3, 3> roty(Scalar_t t) {
        Eigen::Matrix<Scalar_t, 3, 3> R;
        R(0, 0) = std::cos(t);
        R(0, 1) = 0.0;
        R(0, 2) = std::sin(t);
        R(1, 0) = 0.0;
        R(1, 1) = 1.0;
        R(1, 2) = 0;
        R(2, 0) = -std::sin(t);
        R(2, 1) = 0.0;
        R(2, 2) = std::cos(t);

        return R;
    }

    template<typename Scalar_t>
    Eigen::Matrix<Scalar_t, 3, 3> rotz(Scalar_t t) {
        Eigen::Matrix<Scalar_t, 3, 3> R;
        R(0, 0) = std::cos(t);
        R(0, 1) = -std::sin(t);
        R(0, 2) = 0.0;
        R(1, 0) = std::sin(t);
        R(1, 1) = std::cos(t);
        R(1, 2) = 0.0;
        R(2, 0) = 0.0;
        R(2, 1) = 0.0;
        R(2, 2) = 1.0;

        return R;
    }

    template<typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 3, 3> ypr_to_R(const Eigen::DenseBase<Derived> &ypr) {
        EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
        EIGEN_STATIC_ASSERT(Derived::RowsAtCompileTime == 3, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);
        EIGEN_STATIC_ASSERT(Derived::ColsAtCompileTime == 1, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);

        typename Derived::Scalar c, s;

        Eigen::Matrix<typename Derived::Scalar, 3, 3> Rz = Eigen::Matrix<typename Derived::Scalar, 3, 3>::Zero();
        typename Derived::Scalar y = ypr(0);
        c = cos(y);
        s = sin(y);
        Rz(0, 0) = c;
        Rz(1, 0) = s;
        Rz(0, 1) = -s;
        Rz(1, 1) = c;
        Rz(2, 2) = 1;

        Eigen::Matrix<typename Derived::Scalar, 3, 3> Ry = Eigen::Matrix<typename Derived::Scalar, 3, 3>::Zero();
        typename Derived::Scalar p = ypr(1);
        c = cos(p);
        s = sin(p);
        Ry(0, 0) = c;
        Ry(2, 0) = -s;
        Ry(0, 2) = s;
        Ry(2, 2) = c;
        Ry(1, 1) = 1;

        Eigen::Matrix<typename Derived::Scalar, 3, 3> Rx = Eigen::Matrix<typename Derived::Scalar, 3, 3>::Zero();
        typename Derived::Scalar r = ypr(2);
        c = cos(r);
        s = sin(r);
        Rx(1, 1) = c;
        Rx(2, 1) = s;
        Rx(1, 2) = -s;
        Rx(2, 2) = c;
        Rx(0, 0) = 1;

        Eigen::Matrix<typename Derived::Scalar, 3, 3> R = Rz * Ry * Rx;
        return R;
    }

    template<typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 3, 1> R_to_ypr(const Eigen::DenseBase<Derived> &R) {
        EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
        EIGEN_STATIC_ASSERT(Derived::RowsAtCompileTime == 3, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);
        EIGEN_STATIC_ASSERT(Derived::ColsAtCompileTime == 3, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);

        Eigen::Matrix<typename Derived::Scalar, 3, 1> n = R.col(0);
        Eigen::Matrix<typename Derived::Scalar, 3, 1> o = R.col(1);
        Eigen::Matrix<typename Derived::Scalar, 3, 1> a = R.col(2);

        Eigen::Matrix<typename Derived::Scalar, 3, 1> ypr(3);
        typename Derived::Scalar y = atan2(n(1), n(0));
        typename Derived::Scalar p = atan2(-n(2), n(0) * cos(y) + n(1) * sin(y));
        typename Derived::Scalar r =
                atan2(a(0) * sin(y) - a(1) * cos(y), -o(0) * sin(y) + o(1) * cos(y));
        ypr(0) = y;
        ypr(1) = p;
        ypr(2) = r;

        return ypr;
    }

    template<typename Derived>
    Eigen::Quaternion<typename Derived::Scalar> ypr_to_quaternion(const Eigen::DenseBase<Derived> &ypr) {
        EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
        EIGEN_STATIC_ASSERT(Derived::RowsAtCompileTime == 3, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);
        EIGEN_STATIC_ASSERT(Derived::ColsAtCompileTime == 1, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);

        const typename Derived::Scalar cy = cos(ypr(0) / 2.0);
        const typename Derived::Scalar sy = sin(ypr(0) / 2.0);
        const typename Derived::Scalar cp = cos(ypr(1) / 2.0);
        const typename Derived::Scalar sp = sin(ypr(1) / 2.0);
        const typename Derived::Scalar cr = cos(ypr(2) / 2.0);
        const typename Derived::Scalar sr = sin(ypr(2) / 2.0);

        Eigen::Quaternion<typename Derived::Scalar> q;

        q.w() = cr * cp * cy + sr * sp * sy;
        q.x() = sr * cp * cy - cr * sp * sy;
        q.y() = cr * sp * cy + sr * cp * sy;
        q.z() = cr * cp * sy - sr * sp * cy;

        return q;
    }

    template<typename Scalar_t>
    Eigen::Matrix<Scalar_t, 3, 1> quaternion_to_ypr(const Eigen::Quaternion<Scalar_t> &q_) {
        Eigen::Quaternion<Scalar_t> q = q_.normalized();

        Eigen::Matrix<Scalar_t, 3, 1> ypr;
        ypr(2) = atan2(2 * (q.w() * q.x() + q.y() * q.z()), 1 - 2 * (q.x() * q.x() + q.y() * q.y()));
        ypr(1) = asin(2 * (q.w() * q.y() - q.z() * q.x()));
        ypr(0) = atan2(2 * (q.w() * q.z() + q.x() * q.y()), 1 - 2 * (q.y() * q.y() + q.z() * q.z()));

        return ypr;
    }

    template<typename Scalar_t>
    Scalar_t get_yaw_from_quaternion(const Eigen::Quaternion<Scalar_t> &q) {
        return quaternion_to_ypr(q)(0);
    }

    template<typename Scalar_t>
    Eigen::Quaternion<Scalar_t> yaw_to_quaternion(Scalar_t yaw) {
        return Eigen::Quaternion<Scalar_t>(rotz(yaw));
    }

    template<typename Scalar_t>
    Scalar_t normalize_angle(Scalar_t a) {
        int cnt = 0;
        while (true) {
            cnt++;

            if (a < -M_PI) {
                a += M_PI * 2.0;
            } else if (a > M_PI) {
                a -= M_PI * 2.0;
            }

            if (-M_PI <= a && a <= M_PI) {
                break;
            };

            assert(cnt < 10 && "[uav_utils/geometry_msgs] INVALID INPUT ANGLE");
        }

        return a;
    }

    template<typename Scalar_t>
    Scalar_t angle_add(Scalar_t a, Scalar_t b) {
        Scalar_t c = a + b;
        c = normalize_angle(c);
        assert(-M_PI <= c && c <= M_PI);
        return c;
    }

    template<typename Scalar_t>
    Scalar_t yaw_add(Scalar_t a, Scalar_t b) {
        return angle_add(a, b);
    }

    template<typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 3, 3> get_skew_symmetric(const Eigen::DenseBase<Derived> &v) {
        EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
        EIGEN_STATIC_ASSERT(Derived::RowsAtCompileTime == 3, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);
        EIGEN_STATIC_ASSERT(Derived::ColsAtCompileTime == 1, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);

        Eigen::Matrix<typename Derived::Scalar, 3, 3> M;
        M.setZero();
        M(0, 1) = -v(2);
        M(0, 2) = v(1);
        M(1, 0) = v(2);
        M(1, 2) = -v(0);
        M(2, 0) = -v(1);
        M(2, 1) = v(0);
        return M;
    }

    template<typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 3, 1> from_skew_symmetric(const Eigen::DenseBase<Derived> &M) {
        EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
        EIGEN_STATIC_ASSERT(Derived::RowsAtCompileTime == 3, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);
        EIGEN_STATIC_ASSERT(Derived::ColsAtCompileTime == 3, THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);

        Eigen::Matrix<typename Derived::Scalar, 3, 1> v;
        v(0) = M(2, 1);
        v(1) = -M(2, 0);
        v(2) = M(1, 0);

        assert(v.isApprox(Eigen::Matrix<typename Derived::Scalar, 3, 1>(-M(1, 2), M(0, 2), -M(0, 1))));

        return v;
    }


}  // end of namespace uav_utils

#endif
