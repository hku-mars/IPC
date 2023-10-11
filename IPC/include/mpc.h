#ifndef MPC_H_
#define MPC_H_

#include <ros/ros.h>

#include <Eigen/Eigen>
#include "OsqpEigen/OsqpEigen.h"

/* QP formulation:
    min 1/2* x^T H x + f^T x   subject to
    b <= Ax <= b (Ax = b),  d <= Ax <= f,  l <= x <= u
*/
struct mpc_osqp_t {
    Eigen::MatrixXd Ax;
    Eigen::MatrixXd Bx;
    Eigen::MatrixXd M;
    Eigen::MatrixXd C;
    Eigen::MatrixXd Q_bar;
    Eigen::MatrixXd R_bar, R_con_bar;

    Eigen::VectorXd u_low, u_upp, a_low, a_upp, v_low, v_upp;
    Eigen::VectorXd B_a, B_v, B_p;
    Eigen::MatrixXd A_a, A_v, A_p;
    Eigen::MatrixXd M_a, M_v, M_p;

    Eigen::MatrixXd T;
    Eigen::VectorXd D_T;

    Eigen::MatrixXd A_sys;
    Eigen::VectorXd A_sys_low, A_sys_upp;
    Eigen::MatrixXd A_sfc;
    Eigen::VectorXd A_sfc_low, A_sfc_upp;

    Eigen::MatrixXd H;
    Eigen::VectorXd f;
    Eigen::MatrixXd A;
    Eigen::VectorXd Alow, Aupp;
    Eigen::SparseMatrix<double> H_sparse;
    Eigen::SparseMatrix<double> A_sparse;

    Eigen::VectorXd u_optimal;
};

class MPCPlannerClass {
public:
    MPCPlannerClass(){}
    MPCPlannerClass(ros::NodeHandle& nh) {
        Drag_.setZero();

        nh.param("/ipc_node/mpc/horizon", MPC_HORIZON, 5);
        nh.param("/ipc_node/mpc/step",    MPC_STEP, 0.1);
        nh.param("/ipc_node/mpc/ctrl_delay", ctrl_delay_, 0);
        u_last_.resize(ctrl_delay_);
        for (int i = 0; i < u_last_.size(); i++) {
            u_last_[i] = Eigen::Vector3d::Zero();
        }

        nh.param("/ipc_node/mpc/R_p",  R_p_, 100.0);
        nh.param("/ipc_node/mpc/R_v",  R_v_, 0.0);
        nh.param("/ipc_node/mpc/R_a",  R_a_, 0.0);
        nh.param("/ipc_node/mpc/R_u",  R_u_, 10.0);
        nh.param("/ipc_node/mpc/R_u_con",  R_u_con_, 1.0);
        nh.param("/ipc_node/mpc/R_pN", R_pN_, 0.0);
        nh.param("/ipc_node/mpc/R_vN", R_vN_, 0.0);
        nh.param("/ipc_node/mpc/R_aN", R_aN_, 0.0);

        nh.param("/ipc_node/mpc/D_x", Drag_(0, 0), 0.0);
        nh.param("/ipc_node/mpc/D_y", Drag_(1, 1), 0.0);
        nh.param("/ipc_node/mpc/D_z", Drag_(2, 2), 0.0);

        nh.param("/ipc_node/mpc/vx_min", v_min_.x(), -1.0);
        nh.param("/ipc_node/mpc/vy_min", v_min_.y(), -1.0);
        nh.param("/ipc_node/mpc/vz_min", v_min_.z(), -1.0);
        nh.param("/ipc_node/mpc/vx_max", v_max_.x(),  1.0);
        nh.param("/ipc_node/mpc/vy_max", v_max_.y(),  1.0);
        nh.param("/ipc_node/mpc/vz_max", v_max_.z(),  1.0);

        nh.param("/ipc_node/mpc/ax_min", a_min_.x(), -1.0);
        nh.param("/ipc_node/mpc/ay_min", a_min_.y(), -1.0);
        nh.param("/ipc_node/mpc/az_min", a_min_.z(), -1.0);
        nh.param("/ipc_node/mpc/ax_max", a_max_.x(),  1.0);
        nh.param("/ipc_node/mpc/ay_max", a_max_.y(),  1.0);
        nh.param("/ipc_node/mpc/az_max", a_max_.z(),  1.0);

        nh.param("/ipc_node/mpc/ux_min", u_min_.x(), -1.0);
        nh.param("/ipc_node/mpc/uy_min", u_min_.y(), -1.0);
        nh.param("/ipc_node/mpc/uz_min", u_min_.z(), -1.0);
        nh.param("/ipc_node/mpc/ux_max", u_max_.x(),  1.0);
        nh.param("/ipc_node/mpc/uy_max", u_max_.y(),  1.0);
        nh.param("/ipc_node/mpc/uz_max", u_max_.z(),  1.0);

        ProblemFormation();
        X_0_.resize(mpc_.M.cols(), 1);
        X_r_.resize(mpc_.M.rows(), 1);
        planes_.resize(MPC_HORIZON);
    }
    ~MPCPlannerClass() {}

    void GetOptimCmd(Eigen::Vector3d& u, int segment) {
        if (segment >= MPC_HORIZON) {
            segment = MPC_HORIZON - 1;
        }
        u.x() = mpc_.u_optimal(segment * 3 + 0, 0);
        u.y() = mpc_.u_optimal(segment * 3 + 1, 0);
        u.z() = mpc_.u_optimal(segment * 3 + 2, 0);
    }
    void SetGoal(Eigen::Vector3d pr, Eigen::Vector3d vr, Eigen::Vector3d ar, int step) {
        if (step > MPC_HORIZON || step < 0) {
            ROS_WARN("[MPC]: Check goal index! Error index: %d", step);
            return ;
        }
        Eigen::VectorXd x_0(X_0_.rows(), 1);
        x_0.block(0, 0, pr.rows(), 1) = pr;
        x_0.block(pr.rows(), 0, vr.rows(), 1) = vr;
        x_0.block(pr.rows()+vr.rows(), 0, ar.rows(), 1) = ar;
        X_r_.block(x_0.rows()*step, 0, x_0.rows(), 1) = x_0;
    }
    void SetStatus(Eigen::Vector3d p0, Eigen::Vector3d v0, Eigen::Vector3d a0) {
        StatusSaturation(v0, a0);
        X_0_.block(0, 0, p0.rows(), 1) = p0;
        X_0_.block(p0.rows(), 0, v0.rows(), 1) = v0;
        X_0_.block(p0.rows()+v0.rows(), 0, a0.rows(), 1) = a0;
        // for (int i = 0; i < ctrl_delay_; i++) {
        //     X_0_ = mpc_.Ax * X_0_ + mpc_.Bx * u_last_[i];
        // }
    }
    void StatusSaturation(Eigen::Vector3d& v0, Eigen::Vector3d& a0) {
        for (int i = 0; i < 3; i++) {
            if (v0(i, 0) > v_max_(i, 0)) v0(i, 0) = v_max_(i, 0);
            if (v0(i, 0) < v_min_(i, 0)) v0(i, 0) = v_min_(i, 0);
            if (a0(i, 0) > a_max_(i, 0)) a0(i, 0) = a_max_(i, 0);
            if (a0(i, 0) < a_min_(i, 0)) a0(i, 0) = a_min_(i, 0);
        }
    }
    void UpdateOutputHistory(Eigen::Vector3d u) {
        if (ctrl_delay_ == 0) return ;
        u_last_.erase(u_last_.begin());
        u_last_.push_back(u);
    }

    bool Run(void);
    void SystemModel(Eigen::MatrixXd& A, Eigen::MatrixXd& B, double t);
    void SetFSC(Eigen::Matrix<double, Eigen::Dynamic, 4>& planes, int step);
    bool IsInFSC(Eigen::Vector3d pos, Eigen::Matrix<double, Eigen::Dynamic, 4>& planes){
        if (planes.rows() == 0) return false;
        for (int i = 0; i < planes.rows(); i++) {
            if (pos.x()*planes(i,0) + pos.y()*planes(i,1) + 
                pos.z()*planes(i,2) + planes(i,3) > 0) {
                // std::cout << "plane:" << i << " " << planes_.block<1, 4>(i, 0) << std::endl;
                return false;
            }
        }
        return true;
    }

    Eigen::VectorXd X_0_, X_r_;

    // mpc param
    int MPC_HORIZON;
    double MPC_STEP;
    int ctrl_delay_;
    std::vector<Eigen::Vector3d> u_last_;

private:
    void ProblemFormation(void);
    void MPCModel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, 
                  Eigen::MatrixXd& M, Eigen::MatrixXd& C);
    void QuadraticTerm(mpc_osqp_t& mpc, const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R, 
                        const Eigen::MatrixXd& R_con, const Eigen::MatrixXd& F);
    void LinearTerm(mpc_osqp_t& mpc, const Eigen::VectorXd& x_0, const Eigen::VectorXd& x_r);
    void ALLConstraint(mpc_osqp_t& mpc);
    void UpdateBound(mpc_osqp_t& mpc, const Eigen::VectorXd& x_0);

    mpc_osqp_t mpc_;
    ros::Time print_time_;
    int fps_;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 4>> planes_;

    // mpc param
    double R_p_, R_v_, R_a_, R_u_, R_u_con_, R_pN_, R_vN_, R_aN_;
    Eigen::Matrix3d Drag_;
    Eigen::Vector3d v_min_, v_max_, a_min_, a_max_, u_min_, u_max_;
};

#endif
