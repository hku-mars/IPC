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

/* This is an old version of FIRI for temporary usage here. */

#ifndef EMVP_HPP
#define EMVP_HPP

#include <Eigen/Eigen>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <vector>
#include "queue"

#include "lbfgs.hpp"
#include "sdlp.hpp"
#include "geo_utils.hpp"
#include "polytope.hpp"
#include "ellipsoid.hpp"

// #include "utils/scope_timer.hpp"
// #include "fmt/color.h"

// #include <Eigen/SVD>
// #include <Eigen/Dense>

namespace emvp {

/**
 * @brief chol3d
 *	An implementation of 3D cholesky decomposition
 *   https://en.wikipedia.org/wiki/Cholesky_decomposition
 */
    inline void chol3d(const Eigen::Matrix3d &A,
                       Eigen::Matrix3d &L) {
        L(0, 0) = sqrt(A(0, 0));
        L(0, 1) = 0.0;
        L(0, 2) = 0.0;
        L(1, 0) = 0.5 * (A(0, 1) + A(1, 0)) / L(0, 0);
        L(1, 1) = sqrt(A(1, 1) - L(1, 0) * L(1, 0));
        L(1, 2) = 0.0;
        L(2, 0) = 0.5 * (A(0, 2) + A(2, 0)) / L(0, 0);
        L(2, 1) = (0.5 * (A(1, 2) + A(2, 1)) - L(2, 0) * L(1, 0)) / L(1, 1);
        L(2, 2) = sqrt(A(2, 2) - L(2, 0) * L(2, 0) - L(2, 1) * L(2, 1));
        return;
    }

/**
 * @brief smoothedL1
 *	A smooth L1 loss for high-accuracy constrain.
 *   https://mohitjainweb.files.wordpress.com/2018/03/smoothl1loss.pdf
 */
    inline bool smoothedL1(const double &mu,
                           const double &x,
                           double &f,
                           double &df) {
        if (x < 0.0) {
            return false;
        } else if (x > mu) {
            f = x - 0.5 * mu;
            df = 1.0;
            return true;
        } else {
            const double xdmu = x / mu;
            const double sqrxdmu = xdmu * xdmu;
            const double mumxd2 = mu - 0.5 * x;
            f = mumxd2 * sqrxdmu * xdmu;
            df = sqrxdmu * ((-0.5) * xdmu + 3.0 * mumxd2 / mu);
            return true;
        }
    }

/**
 * @brief costMVIE
 *
 *   The cost function of Maximum Volume Inscribe Ellipsoid problem
 *
 */
    inline double costMVIE(void *data,
                           const Eigen::VectorXd &x,
                           Eigen::VectorXd &grad) {
        const int *pM = (int *) data;
        const double *pSmoothEps = (double *) (pM + 1);
        const double *pPenaltyWt = pSmoothEps + 1;
        const double *pA = pPenaltyWt + 1;
        const int M = *pM;
        const double smoothEps = *pSmoothEps;
        const double penaltyWt = *pPenaltyWt;
        bool *fix_p = (bool *) (pSmoothEps + (2 + 3 * M));
        // 先从 pA中取出所有多面体
        Eigen::Map<const Eigen::MatrixX3d> A(pA, M, 3);
        // 取出seed
        Eigen::Map<const Eigen::Vector3d> p(x.data());
        // 取出  L 矩阵中的非零元素 前三个是对角元素
        Eigen::Map<const Eigen::Vector3d> rtd(x.data() + 3);
        // 后三个是左下角元素
        Eigen::Map<const Eigen::Vector3d> cde(x.data() + 6);
        // 取出seed和半径的梯度
        Eigen::Map<Eigen::Vector3d> gdp(grad.data());
        Eigen::Map<Eigen::Vector3d> gdrtd(grad.data() + 3);
        Eigen::Map<Eigen::Vector3d> gdcde(grad.data() + 6);

        double cost = 0;
        gdp.setZero();
        gdrtd.setZero();
        gdcde.setZero();

        Eigen::Matrix3d L;
        L(0, 0) = rtd(0) * rtd(0) + DBL_EPSILON;
        L(0, 1) = 0.0;
        L(0, 2) = 0.0;
        L(1, 0) = cde(0);
        L(1, 1) = rtd(1) * rtd(1) + DBL_EPSILON;
        L(1, 2) = 0.0;
        L(2, 0) = cde(2);
        L(2, 1) = cde(1);
        L(2, 2) = rtd(2) * rtd(2) + DBL_EPSILON;

        const Eigen::MatrixX3d AL = A * L;
        const Eigen::VectorXd normAL = AL.rowwise().norm();
        const Eigen::Matrix3Xd adjNormAL = (AL.array().colwise() / normAL.array()).transpose();
        const Eigen::VectorXd consViola = (normAL + A * p).array() - 1.0;

        double c, dc;
        Eigen::Vector3d vec;
        for (int i = 0; i < M; ++i) {
            if (smoothedL1(smoothEps, consViola(i), c, dc)) {
                cost += c;
                vec = dc * A.row(i).transpose();
                gdp += vec;
                gdrtd += adjNormAL.col(i).cwiseProduct(vec);
                gdcde(0) += adjNormAL(0, i) * vec(1);
                gdcde(1) += adjNormAL(1, i) * vec(2);
                gdcde(2) += adjNormAL(0, i) * vec(2);
            }
        }
        cost *= penaltyWt;
        gdp *= penaltyWt;
        gdrtd *= penaltyWt;
        gdcde *= penaltyWt;

        cost -= log(L(0, 0)) + log(L(1, 1)) + log(L(2, 2));
        gdrtd(0) -= 1.0 / L(0, 0);
        gdrtd(1) -= 1.0 / L(1, 1);
        gdrtd(2) -= 1.0 / L(2, 2);

        gdrtd(0) *= 2.0 * rtd(0);
        gdrtd(1) *= 2.0 * rtd(1);
        gdrtd(2) *= 2.0 * rtd(2);
        if (*fix_p == true) {
            gdp.setZero();
        }

        return cost;
    }

/**
 * @brief maxVolInsEllipsoid
 * @param hPoly X rows of hyperplane
 * @param offset_x offset added to the long semi-axis, default is 0
 */
    inline bool maxVolInsEllipsoid(const Eigen::MatrixX4d &hPoly,
                                   Ellipsoid &E,
//							   const Eigen::Vector3d &interior,
                                   const bool &_fix_p = false) {
        // Find the deepest interior point [ Anylitical center]
        const int M = hPoly.rows();
        Eigen::MatrixX4d Alp(M, 4);
        Eigen::VectorXd blp(M);
        Eigen::Vector4d clp, xlp;
        const Eigen::ArrayXd hNorm = hPoly.leftCols<3>().rowwise().norm();
        Alp.leftCols<3>() = hPoly.leftCols<3>().array().colwise() / hNorm;
        Alp.rightCols<1>().setConstant(1.0);
        blp = -hPoly.rightCols<1>().array() / hNorm;
        clp.setZero();
        clp(3) = -1.0;
        const double maxdepth = -sdlp::linprog<4>(clp, Alp, blp, xlp);
        if (!(maxdepth > 0.0) || std::isinf(maxdepth)) {
            return false;
        }
        const Eigen::Vector3d interior = xlp.head<3>();
//  const Eigen::Vector3d interior = E.d();
// Prepare the data for MVIE optimization
// Maximum Volume Inscribed Ellipsoid
        uint8_t *optData = new uint8_t[sizeof(int) + (2 + 3 * M) * sizeof(double) + sizeof(bool)];
        int *pM = (int *) optData;
        double *pSmoothEps = (double *) (pM + 1);
        double *pPenaltyWt = pSmoothEps + 1;
        double *pA = pPenaltyWt + 1;
        bool *fix_p = (bool *) (pSmoothEps + (2 + 3 * M));

        *fix_p = _fix_p;
        *pM = M;
        Eigen::Map<Eigen::MatrixX3d> A(pA, M, 3);
        A = Alp.leftCols<3>().array().colwise() /
            (blp - Alp.leftCols<3>() * interior).array();

        Eigen::VectorXd x(9);
        const Eigen::Matrix3d Q = E.R() * (E.r().cwiseProduct(E.r())).asDiagonal() * E.R().transpose();
        Eigen::Matrix3d L;
        chol3d(Q, L);

// seed到多边形中心的距离
        x.head<3>() = E.d() - interior;
// L 矩阵中的非零元素
        x(3) = sqrt(L(0, 0));
        x(4) = sqrt(L(1, 1));
        x(5) = sqrt(L(2, 2));
        x(6) = L(1, 0);
        x(7) = L(2, 1);
        x(8) = L(2, 0);

        double minCost;
        lbfgs::lbfgs_parameter_t paramsMVIE;
        paramsMVIE.mem_size = 18;
        paramsMVIE.g_epsilon = 0.0;
        paramsMVIE.min_step = 1.0e-32;
        paramsMVIE.past = 3;
        paramsMVIE.delta = 1.0e-7;
        *pSmoothEps = 1.0e-2;
        *pPenaltyWt = 1.0e+3;

        int ret = lbfgs::lbfgs_optimize(x,
                                        minCost,
                                        &costMVIE,
                                        nullptr,
                                        nullptr,
                                        optData,
                                        paramsMVIE);

        if (ret < 0) {
            printf("FIRI WARNING: %s\n", lbfgs::lbfgs_strerror(ret));
        }
        Eigen::Vector3d d, r;
        Eigen::Matrix3d C;
        d = x.head<3>() + interior;
        L(0, 0) = x(3) * x(3);
        L(0, 1) = 0.0;
        L(0, 2) = 0.0;
        L(1, 0) = x(6);
        L(1, 1) = x(4) * x(4);
        L(1, 2) = 0.0;
        L(2, 0) = x(8);
        L(2, 1) = x(7);
        L(2, 2) = x(5) * x(5);
        Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::FullPivHouseholderQRPreconditioner> svd2(L, Eigen::ComputeFullU);
        const Eigen::Matrix3d U = svd2.matrixU();
        const Eigen::Vector3d S = svd2.singularValues();
        if (U.determinant() < 0.0) {
            C.col(0) = U.col(1);
            C.col(1) = U.col(0);
            C.col(2) = U.col(2);
            r(0) = S(1);
            r(1) = S(0);
            r(2) = S(2);
        } else {
            C = U;
            r = S;
        }
        E = Ellipsoid(C, r, d);
        delete[] optData;
        return ret >= 0;
    }

/**
 * @brief findEllipsoid: find maximum ellipsoid with RILS
 * @param pc the obstacle points
 * @param a the start point of the line segment seed
 * @param b the end point of the line segment seed
 * @param out_ell the output ellipsoid
 * @param r_robot the robot_size, decide if the polytope need to be shrink
 * @param _fix_p decide if the ellipsoid center need to be optimized
 * @param iterations number of the alternating optimization
 */

    inline void findEllipsoid(
            const Eigen::Matrix3Xd &pc,
            const Eigen::Vector3d &a,
            const Eigen::Vector3d &b,
            Ellipsoid &out_ell,
            const double &robot_r = 0.0) {

        double f = (a - b).norm() / 2;
        Mat3f C = f * Mat3f::Identity();
        Vec3f r = Vec3f::Constant(f);
        Vec3f center = (a + b) / 2;
        C(0, 0) += robot_r;
        r(0) += robot_r;
        if (r(0) > 0) {
            double ratio = r(1) / r(0);
            r *= ratio;
            C *= ratio;
        }

        Mat3f Ri = Eigen::Quaterniond::FromTwoVectors(Vec3f::UnitX(), (b - a)).toRotationMatrix();
        Ellipsoid E(Ri, r, center);
        Mat3f Rf = Ri;
        Mat3Df obs;
        int min_dis_id;
        Vec3f pw;
        if (E.pointsInside(pc, obs, min_dis_id)) {
            pw = obs.col(min_dis_id);
        } else {
            out_ell = E;
            return;
        }
        Mat3Df obs_inside = obs;
        int max_iter = 100;
        while (max_iter--) {
            Vec3f p_e = Ri.transpose() * (pw - E.d());
            const double roll = atan2(p_e(2), p_e(1));
            Rf = Ri * Eigen::Quaterniond(cos(roll / 2), sin(roll / 2), 0, 0);
            p_e = Rf.transpose() * (pw - E.d());
            if (p_e(0) < r(0)) {
                r(1) = std::abs(p_e(1)) / std::sqrt(1 - std::pow(p_e(0) / r(0), 2));
            }
            E = Ellipsoid(Rf, r, center);
            if (E.pointsInside(obs_inside, obs_inside, min_dis_id)) {
                pw = obs_inside.col(min_dis_id);
            } else {
                break;
            }
        }
        if (max_iter == 0) {
            print(fg(color::gold), " -- [EMVP] Find Ellipsoid reach max iteration, may cause error.\n");
        }
        max_iter = 100;


        if (E.pointsInside(obs, obs_inside, min_dis_id)) {
            pw = obs_inside.col(min_dis_id);
        } else {
            out_ell = E;
            return;
        }

        while (max_iter--) {
            Vec3f p = Rf.transpose() * (pw - E.d());
            double dd = 1 - std::pow(p(0) / r(0), 2) -
                        std::pow(p(1) / r(1), 2);
            if (dd > epsilon_) {
                r(2) = std::abs(p(2)) / std::sqrt(dd);
            }
            E = Ellipsoid(Rf, r, center);
            if (E.pointsInside(obs_inside, obs_inside, min_dis_id)) {
                pw = obs_inside.col(min_dis_id);
            } else {
                out_ell = E;
                break;
            }
        }

        if (max_iter == 0) {
            print(fg(color::gold), " -- [EMVP] Find Ellipsoid reach max iteration, may cause error.\n");
        }
        E = Ellipsoid(Rf, r, center);
        out_ell = E;
    }

    inline void findTangentPlaneOfSphere(const Eigen::Vector3d &center, const double &r,
                                         const Eigen::Vector3d &pass_point,
                                         const Eigen::Vector3d &seed_p,
                                         Eigen::Vector4d &outter_plane) {
        Eigen::Vector3d P = pass_point - center;
        Eigen::Vector3d norm_ = (pass_point - center).cross(seed_p - center).normalized();
        Eigen::Matrix3d R = Eigen::Quaterniond::FromTwoVectors(norm_, Vec3f(0, 0, 1)).matrix();
        P = R * P;
        Eigen::Vector3d C = R * (seed_p - center);
        Eigen::Vector3d Q;
        double r2 = r * r;
        double p1p2n = P.head(2).squaredNorm();
        double d = sqrt(p1p2n - r2);
        double rp1p2n = r / p1p2n;
        double q11 = rp1p2n * (P(0) * r - P(1) * d);
        double q21 = rp1p2n * (P(1) * r + P(0) * d);

        double q12 = rp1p2n * (P(0) * r + P(1) * d);
        double q22 = rp1p2n * (P(1) * r - P(0) * d);
        if (q11 * C(0) + q21 * C(1) > 0) {
            Q(0) = q12;
            Q(1) = q22;
        } else {
            Q(0) = q11;
            Q(1) = q21;
        }
        Q(2) = 0;
        // point(Q) + normal (AQ)
        outter_plane.head(3) = R.transpose() * Q;
        Q = outter_plane.head(3) + center;
        outter_plane(3) = -Q.dot(outter_plane.head(3));
        if (outter_plane.head(3).dot(center) + outter_plane(3) > epsilon_) {
            outter_plane = -outter_plane;
        }
    }

/**
 * @brief embed maximum volume polytope
 * @param bd bounding box with 6 faces
 * @param pc the obstacle points
 * @param a the start point of the line segment seed
 * @param b the end point of the line segment seed
 * @param hPoly the output polytope
 * @param r_robot the robot_size, decide if the polytope need to be shrink
 * @param _fix_p decide if the ellipsoid center need to be optimized
 * @param iterations number of the alternating optimization
 */
    bool emvp(const Eigen::MatrixX4d &bd,
              const Eigen::Matrix3Xd &pc,
              const Eigen::Vector3d &a,
              const Eigen::Vector3d &b,
              Polytope &out_poly,
              const double r_robot = 0.1,
              const bool _fix_p = false,
              const int iterations = 3) {
        const Eigen::Vector4d ah(a(0), a(1), a(2), 1.0);
        const Eigen::Vector4d bh(b(0), b(1), b(2), 1.0);

        /// planes
        MatD4f hPoly;

        /// force return if the seed is not inside the boundary
        if ((bd * ah).maxCoeff() > 0.0 ||
            (bd * bh).maxCoeff() > 0.0) {
            print(fg(color::gold), " -- [WARN] ah, bh not in BD, forced return.\n");
            return false;
        }

        /// Maximum M boundary constraints and N point constraints
        const int M = bd.rows();
        const int N = pc.cols();

        Ellipsoid E(Mat3f::Identity(), (a + b) / 2);
        if ((a - b).norm() > 0.1) {
            /// use line seed
            findEllipsoid(pc, a, b, E, r_robot);
        }

        vector<Eigen::Vector4d> planes;
        for (int loop = 0; loop < iterations; ++loop) {
            // Initialize the boundary in ellipsoid frame
            const Eigen::MatrixX4d bd_e = E.toEllipsoidFrame(bd);
            // Initialize the seed points
            const Eigen::Vector3d fwd_a = E.toEllipsoidFrame(a);
            const Eigen::Vector3d fwd_b = E.toEllipsoidFrame(b);
            const Eigen::VectorXd distDs = bd_e.rightCols<1>().cwiseAbs().cwiseQuotient(
                    bd_e.leftCols<3>().rowwise().norm());
            const Eigen::Matrix3Xd pc_e = E.toEllipsoidFrame(pc);
            Eigen::VectorXd distRs = pc_e.colwise().norm();
            Eigen::Matrix<uint8_t, -1, 1> bdFlags = Eigen::Matrix<uint8_t, -1, 1>::Constant(M, 1);
            Eigen::Matrix<uint8_t, -1, 1> pcFlags = Eigen::Matrix<uint8_t, -1, 1>::Constant(N, 1);

            bool completed = false;
            int bdMinId, pcMinId;
            double minSqrD = distDs.minCoeff(&bdMinId);
            double minSqrR = distRs.minCoeff(&pcMinId);

            Eigen::Vector3d pt_w, pt_e;
            Eigen::Vector4d temp_tangent, temp_tange_W;

            planes.clear();
            planes.reserve(30);
            for (int i = 0; !completed && i < (M + N); ++i) {
                // Get the min dis point of this round.
                pt_w = pc.col(pcMinId);
                pt_e = pc_e.col(pcMinId);
                if (minSqrD < minSqrR) {
                    // enable the boundary constrain.
                    temp_tangent = bd_e.row(bdMinId);
                    bdFlags(bdMinId) = 0;
                } else {
                    // enable the obstacle point constarin.
                    // First generate a plane in E frame
                    temp_tangent(3) = -distRs(pcMinId);
                    temp_tangent.head(3) = pt_e.transpose() / distRs(pcMinId);
                    if (r_robot < epsilon_) {
                        if (temp_tangent.head(3).dot(fwd_a) + temp_tangent(3) > epsilon_) {
                            const Eigen::Vector3d delta = pc_e.col(pcMinId) - fwd_a;
                            temp_tangent.head(3) = fwd_a - (delta.dot(fwd_a) / delta.squaredNorm()) * delta;
                            distRs(pcMinId) = temp_tangent.head(3).norm();
                            temp_tangent(3) = -distRs(pcMinId);
                            temp_tangent.head(3) /= distRs(pcMinId);
                        }
                        if (temp_tangent.head(3).dot(fwd_b) + temp_tangent(3) > epsilon_) {
                            const Eigen::Vector3d delta = pc_e.col(pcMinId) - fwd_b;
                            temp_tangent.head(3) = fwd_b - (delta.dot(fwd_b) / delta.squaredNorm()) * delta;
                            distRs(pcMinId) = temp_tangent.head(3).norm();
                            temp_tangent(3) = -distRs(pcMinId);
                            temp_tangent.head(3) /= distRs(pcMinId);
                        }
                        if (temp_tangent.head(3).dot(fwd_b) + temp_tangent(3) > epsilon_) {
                            const Eigen::Vector3d delta = pc_e.col(pcMinId) - fwd_b;
                            temp_tangent.head(3) = fwd_b - (delta.dot(fwd_b) / delta.squaredNorm()) * delta;
                            distRs(pcMinId) = temp_tangent.head(3).norm();
                            temp_tangent(3) = -distRs(pcMinId);
                            temp_tangent.head(3) /= distRs(pcMinId);
                        }
                    } else {
                        // Then convert the plane to world frame
                        temp_tange_W = E.toWorldFrame(temp_tangent);
                        // the check in the w frame
                        if (temp_tange_W.head(3).dot(a) + temp_tange_W(3) + r_robot > epsilon_) {
                            findTangentPlaneOfSphere(a, r_robot, pt_w, E.d(), temp_tange_W);
                        }
                        if (temp_tange_W.head(3).dot(b) + temp_tange_W(3) + r_robot > epsilon_) {
                            findTangentPlaneOfSphere(b, r_robot, pt_w, E.d(), temp_tange_W);
                        }
                        temp_tangent = E.toEllipsoidFrame(temp_tange_W);
                    }
//		distRs(pcMinId) = -temp_tangent(3);
                    pcFlags(pcMinId) = 0;
                }
                // update pcMinId and bdMinId
                completed = true;
                minSqrD = INFINITY;
                for (int j = 0; j < M; ++j) {
                    if (bdFlags(j)) {
                        completed = false;
                        if (minSqrD > distDs(j)) {
                            bdMinId = j;
                            minSqrD = distDs(j);
                        }
                    }
                }
                minSqrR = INFINITY;
                for (int j = 0; j < N; ++j) {
                    if (pcFlags(j)) {
                        if (temp_tangent.head(3).dot(pc_e.col(j)) + temp_tangent(3) > -epsilon_) {
                            pcFlags(j) = 0;
                        } else {
                            completed = false;
                            if (minSqrR > distRs(j)) {
                                pcMinId = j;
                                minSqrR = distRs(j);
                            }
                        }
                    }
                }
                planes.push_back(temp_tangent);
            }
            hPoly.resize(planes.size(), 4);
            for (int i = 0; i < planes.size(); ++i) {
                hPoly.row(i) = E.toWorldFrame(planes[i]);
            }
            if (loop == iterations - 1) {
                break;
            }
            maxVolInsEllipsoid(hPoly, E, _fix_p);
        }

        /// shrink the polytope with robot_r
        for (int i = 0; i < hPoly.rows(); i++) {
            // A B C D / n + r )* n
            double n = hPoly.row(i).head(3).norm();
            hPoly(i, 3) += r_robot * n;
        }
        out_poly.SetPlanes(hPoly);
        out_poly.SetSeedLine(make_pair(a, b));
        out_poly.SetEllipsoid(E);
        return true;
    }

}

#endif