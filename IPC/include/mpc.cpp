#include "mpc.h"

bool MPCPlannerClass::Run(void)
{
    ros::Time time_0 = ros::Time::now();

    // update reference
    LinearTerm(mpc_, X_0_, X_r_);

    // update system status and their constraints
    UpdateBound(mpc_, X_0_);
    
    // update safe flight corridor constraints
    bool add_sfc_flag = true;
    if (mpc_.T.rows() == 0) add_sfc_flag = false;
    if (add_sfc_flag) {
        mpc_.A_sfc.resize(mpc_.T.rows(), mpc_.A_p.cols());
        mpc_.A_sfc_upp.resize(mpc_.D_T.rows(), 1);
        mpc_.A_sfc = mpc_.T * mpc_.A_p;
        mpc_.A_sfc_upp = -mpc_.D_T - mpc_.T * mpc_.B_p;
    }

    // generate all inequality constraints
    if (add_sfc_flag) {
        mpc_.A.resize(mpc_.A_sys.rows() + mpc_.A_sfc.rows(), mpc_.A_sys.cols());
        mpc_.A.block(0, 0, mpc_.A_sys.rows(), mpc_.A_sys.cols()) = mpc_.A_sys;
        mpc_.A.block(mpc_.A_sys.rows(), 0, mpc_.A_sfc.rows(), mpc_.A_sfc.cols()) = mpc_.A_sfc;
        mpc_.Alow.resize(mpc_.A_sys_low.rows() + mpc_.A_sfc_low.rows(), 1);
        mpc_.Alow.block(0, 0, mpc_.A_sys_low.rows(), 1) = mpc_.A_sys_low;
        mpc_.Alow.block(mpc_.A_sys_low.rows(), 0, mpc_.A_sfc_low.rows(), 1) = mpc_.A_sfc_low;
        mpc_.Aupp.resize(mpc_.A_sys_upp.rows() + mpc_.A_sfc_upp.rows(), 1);
        mpc_.Aupp.block(0, 0, mpc_.A_sys_upp.rows(), 1) = mpc_.A_sys_upp;
        mpc_.Aupp.block(mpc_.A_sys_upp.rows(), 0, mpc_.A_sfc_upp.rows(), 1) = mpc_.A_sfc_upp;
    } else {
        mpc_.A.resize(mpc_.A_sys.rows(), mpc_.A_sys.cols());
        mpc_.A.block(0, 0, mpc_.A_sys.rows(), mpc_.A_sys.cols()) = mpc_.A_sys;
        mpc_.Alow.resize(mpc_.A_sys_low.rows(), 1);
        mpc_.Alow.block(0, 0, mpc_.A_sys_low.rows(), 1) = mpc_.A_sys_low;
        mpc_.Aupp.resize(mpc_.A_sys_upp.rows(), 1);
        mpc_.Aupp.block(0, 0, mpc_.A_sys_upp.rows(), 1) = mpc_.A_sys_upp;
    }

    mpc_.H_sparse = mpc_.H.sparseView();
    mpc_.A_sparse = mpc_.A.sparseView();

    // data reset
    mpc_.T.resize(0, 0);
    mpc_.A_sfc.resize(0, 0);
    mpc_.A_sfc_low.resize(0, 1);
    mpc_.A_sfc_upp.resize(0, 1);
    mpc_.D_T.resize(0, 1);

    // use osqp-eigen to solve MPC problem
    OsqpEigen::Solver solver;
    // solver.settings()->setTimeLimit(0.008);
    solver.settings()->setVerbosity(0); // osqp stop print
    solver.settings()->setWarmStart(true);
    // solver.setWarmStart()
    solver.data()->setNumberOfConstraints(mpc_.A_sparse.rows());
    solver.data()->setNumberOfVariables(mpc_.f.rows());
    solver.data()->setHessianMatrix(mpc_.H_sparse);
    solver.data()->setGradient(mpc_.f);
    solver.data()->setLinearConstraintsMatrix(mpc_.A_sparse);
    solver.data()->setLowerBound(mpc_.Alow);
    solver.data()->setUpperBound(mpc_.Aupp);

    bool init_flag = solver.initSolver();
    bool solve_flag = true;
    if (init_flag) {
        solve_flag = solver.solve();
    } else {
        ROS_ERROR("[MPC]: Can't set mpc problem!");
        solve_flag = false;
    }
    fps_++;
    if (solve_flag == true && init_flag == true) {
        mpc_.u_optimal = solver.getSolution();
        static double vel_max = 0.0;
        if (vel_max < X_0_.block(3,0,3,1).norm()) vel_max = X_0_.block(3,0,3,1).norm();
        if ((ros::Time::now()-print_time_).toSec() > 2.0) {
            print_time_ = ros::Time::now();
            std::cout << "mpc fps: " << fps_/2 << ", this time is: " << (ros::Time::now()-time_0).toSec()*1000 << " ms. " 
                "Velocity now is: " << vel_max << "m/s. " << std::endl;
            fps_ = 0;
        }
        return true;
    } else {
        bool flag = IsInFSC(X_0_.block(0,0,3,1), planes_[0]);
        // std::cout << "planes: " << std::endl << planes_ << std::endl;
        std::cout << "pos: " << X_0_.block(0,0,3,1).transpose() << "  fsc: " << std::boolalpha << flag << std::endl;
        if (init_flag) {
            OsqpEigen::Status status = solver.getStatus();
            if (status == OsqpEigen::Status::DualInfeasibleInaccurate) {
                ROS_ERROR("[MPC]: Error status: Dual Infeasible Inaccurate");
            }
            if (status == OsqpEigen::Status::PrimalInfeasibleInaccurate) {
                ROS_ERROR("[MPC]: Error status: Primal Infeasible Inaccurate");
            }
            if (status == OsqpEigen::Status::SolvedInaccurate) {
                ROS_ERROR("[MPC]: Error status: Solved Inaccurate");
            }
            if (status == OsqpEigen::Status::Sigint) {
                ROS_ERROR("[MPC]: Error status: Sigint");
            }
            if (status == OsqpEigen::Status::MaxIterReached) {
                ROS_ERROR("[MPC]: Error status: Max Iter Reached");
            }
            if (status == OsqpEigen::Status::PrimalInfeasible) {
                ROS_ERROR("[MPC]: Error status: Primal Infeasible");
                // std::cout << "init state: " << X_0_.transpose() << std::endl;
            }
            if (status == OsqpEigen::Status::DualInfeasible) {
                ROS_ERROR("[MPC]: Error status: Dual Infeasible");
            }
            if (status == OsqpEigen::Status::NonCvx) {
                ROS_ERROR("[MPC]: Error status: NonCvx");
            }
        }
        return false;
    }
}

void MPCPlannerClass::SetFSC(Eigen::Matrix<double, Eigen::Dynamic, 4>& planes, int step)
{
    if (step >= MPC_HORIZON || step < 0) {
        ROS_WARN("[MPC]: Check sfc index! Error index: %d", step);
    }
    planes_[step] = planes;
    mpc_.T.conservativeResize(mpc_.T.rows() + planes.rows(), mpc_.A_p.rows());
    mpc_.A_sfc_low.conservativeResize(mpc_.A_sfc_low.rows() + planes.rows(), 1);
    mpc_.D_T.conservativeResize(mpc_.D_T.rows() + planes.rows(), 1);
    for (int i = 0; i < planes.rows(); i++) {
        mpc_.T.row(mpc_.T.rows()-1 - i).setZero();
        mpc_.T.block(mpc_.T.rows()-1 - i, step*3, 1, 3) = Eigen::Vector3d(planes(i, 0), planes(i, 1), planes(i, 2)).transpose(); 
        mpc_.A_sfc_low(mpc_.A_sfc_low.rows()-1 - i, 0) = -OSQP_INFTY;
        mpc_.D_T(mpc_.D_T.rows()-1 - i, 0) = planes(i, 3);
    }
}

void MPCPlannerClass::ProblemFormation(void)
{
    /* system status: {p1, v1, a1, p2, v2, a2, ... , pN, vN, aN} 
       input: {u0, u1, u2, ... , u(N-1)}
    */

    // system model
    SystemModel(mpc_.Ax, mpc_.Bx, MPC_STEP);
    MPCModel(mpc_.Ax, mpc_.Bx, mpc_.M, mpc_.C);

    // cost function (Quadratic term and Linear term)
    Eigen::MatrixXd Q = Eigen::MatrixXd(9, 9).setZero();
    Q.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity() * R_p_;
    Q.block(3, 3, 3, 3) = Eigen::Matrix3d::Identity() * R_v_;
    Q.block(6, 6, 3, 3) = Eigen::Matrix3d::Identity() * R_a_;
    Eigen::MatrixXd R = Eigen::MatrixXd(3, 3).setIdentity() * R_u_;
    Eigen::MatrixXd R_con = Eigen::MatrixXd(3, 3).setIdentity() * R_u_con_;
    Eigen::MatrixXd F = Eigen::MatrixXd(9, 9).setZero();
    F.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity() * R_pN_;
    F.block(3, 3, 3, 3) = Eigen::Matrix3d::Identity() * R_vN_;
    F.block(6, 6, 3, 3) = Eigen::Matrix3d::Identity() * R_aN_;
    QuadraticTerm(mpc_, Q, R, R_con, F);
    LinearTerm(mpc_, Eigen::VectorXd(9, 1).setZero(), Eigen::VectorXd(9*MPC_HORIZON, 1).setZero());

    // system status and input constrains
    ALLConstraint(mpc_);

    // Inequality constrains (none)

    mpc_.u_optimal.resize(mpc_.f.rows(), 1);

    ROS_INFO("\033[1;32mMPC problem formation done!\033[0m");
}

void MPCPlannerClass::SystemModel(Eigen::MatrixXd& A, Eigen::MatrixXd& B, double t)
{
    A.resize(9, 9);
    A.setZero();
    A.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();
    A.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity() * t;
    A.block(0, 6, 3, 3) = Eigen::Matrix3d::Identity() * t * t * 0.5;
    A.block(3, 3, 3, 3) = Eigen::Matrix3d::Identity() * (Eigen::Matrix3d::Identity() - Drag_ * t);
    A.block(3, 6, 3, 3) = Eigen::Matrix3d::Identity() * t;
    A.block(6, 6, 3, 3) = Eigen::Matrix3d::Identity();

    B.resize(9, 3);
    B.setZero();
    B.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity() * 1.0/6.0 * std::pow(t, 3);
    B.block(3, 0, 3, 3) = Eigen::Matrix3d::Identity() * 1.0/2.0 * std::pow(t, 2);
    B.block(6, 0, 3, 3) = Eigen::Matrix3d::Identity() * t;
}

void MPCPlannerClass::MPCModel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, 
                               Eigen::MatrixXd& M, Eigen::MatrixXd& C)
{
    M.resize(MPC_HORIZON * A.rows(), A.cols());
    M.setZero();
    C.resize(MPC_HORIZON * B.rows(), MPC_HORIZON * B.cols());
    C.setZero();

    Eigen::MatrixXd temp = Eigen::MatrixXd(A.rows(), A.cols()).setIdentity();
    for (int i = 0; i < MPC_HORIZON; i++) {
        if (i == 0) {
            C.block(0, 0, B.rows(), B.cols()) = B;
        } else {
            Eigen::MatrixXd temp_c = Eigen::MatrixXd(B.rows(), C.cols());
            temp_c << temp * B, C.block((i-1) * B.rows(), 0, B.rows(), B.cols() * (MPC_HORIZON-1));
            C.block(B.rows() * i, 0, B.rows(), C.cols()) = temp_c;
        }

        temp = temp * A;
        M.block(A.rows() * i, 0, A.rows(), A.cols()) = temp;
    }
}

void MPCPlannerClass::QuadraticTerm(mpc_osqp_t& mpc, const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R, 
                                    const Eigen::MatrixXd& R_con, const Eigen::MatrixXd& F)
{
    mpc.Q_bar.resize(Q.rows() * MPC_HORIZON, Q.cols() * MPC_HORIZON);
    mpc.Q_bar.setZero();
    mpc.R_bar.resize(R.rows() * MPC_HORIZON, R.cols() * MPC_HORIZON);
    mpc.R_bar.setZero();
    mpc.R_con_bar.resize(R_con.rows() * MPC_HORIZON, R_con.cols() * MPC_HORIZON);
    mpc.R_con_bar.setZero();
    for (int i = 0; i < MPC_HORIZON; i++) {
        mpc.Q_bar.block(i * Q.rows(), i * Q.cols(), Q.rows(), Q.cols()) = Q;
        mpc.R_bar.block(i * R.rows(), i * R.cols(), R.rows(), R.cols()) = R;
        if (i == 0) mpc.R_con_bar.block(0, 0, R_con.rows(), R_con.cols()) = R_con;
        else if (i == MPC_HORIZON - 1) {
            mpc.R_con_bar.block(i * R_con.rows(), i * R_con.cols(), R_con.rows(), R_con.cols()) = R_con;
            mpc.R_con_bar.block(i * R_con.rows(), (i-1) * R_con.cols(), R_con.rows(), R_con.cols()) = -2 * R_con;
        } else {
            mpc.R_con_bar.block(i * R_con.rows(), i * R_con.cols(), R_con.rows(), R_con.cols()) = 2 * R_con;
            mpc.R_con_bar.block(i * R_con.rows(), (i-1) * R_con.cols(), R_con.rows(), R_con.cols()) = -2 * R_con;
        }
    }
    mpc.Q_bar.block((MPC_HORIZON-1) * Q.rows(), (MPC_HORIZON-1) * Q.cols(), F.rows(), F.cols()) = F;

    mpc.H.resize(mpc.C.rows(), mpc.C.cols());
    mpc.H = mpc.C.transpose() * mpc.Q_bar * mpc.C + mpc.R_bar + mpc.R_con_bar;
}

void MPCPlannerClass::LinearTerm(mpc_osqp_t& mpc, const Eigen::VectorXd& x_0, const Eigen::VectorXd& x_r)
{
    if (x_r.rows() != mpc.M.rows()) {
        ROS_ERROR("[MPC]: MPC linear term set goal error!");
        return;
    }
    mpc.f.resize(mpc.C.rows(), 1);
    mpc.f = ((x_0.transpose()*mpc.M.transpose() - x_r.transpose()) * mpc.Q_bar * mpc.C).transpose();
}

void MPCPlannerClass::ALLConstraint(mpc_osqp_t& mpc)
{
    mpc.A_sys.resize(3 * MPC_HORIZON * 3, 3 * MPC_HORIZON);
    mpc.A_sys_low.resize(3 * MPC_HORIZON * 3, 1);
    mpc.A_sys_upp.resize(3 * MPC_HORIZON * 3, 1);

    /* --system constraint: input, acceleration, velocity limit-- */
    mpc.u_low.resize(3 * MPC_HORIZON, 1);
    mpc.u_upp.resize(3 * MPC_HORIZON, 1);
    mpc.a_low.resize(3 * MPC_HORIZON, 1);
    mpc.a_upp.resize(3 * MPC_HORIZON, 1);
    mpc.v_low.resize(3 * MPC_HORIZON, 1);
    mpc.v_upp.resize(3 * MPC_HORIZON, 1);
    for (int i = 0; i < MPC_HORIZON; i++) {
        mpc.u_low.block(i * u_min_.rows(), 0, u_min_.rows(), 1) = u_min_;
        mpc.u_upp.block(i * u_max_.rows(), 0, u_max_.rows(), 1) = u_max_;
        mpc.a_low.block(i * a_min_.rows(), 0, a_min_.rows(), 1) = a_min_;
        mpc.a_upp.block(i * a_max_.rows(), 0, a_max_.rows(), 1) = a_max_;
        mpc.v_low.block(i * v_min_.rows(), 0, v_min_.rows(), 1) = v_min_;
        mpc.v_upp.block(i * v_max_.rows(), 0, v_max_.rows(), 1) = v_max_;
    }
    
    // input constraint
    Eigen::MatrixXd A_u;
    A_u.resize(3 * MPC_HORIZON, 3 * MPC_HORIZON);
    A_u.setIdentity();

    // calculate A_p, A_v, A_a: system state position p, v, a transform to input u
    mpc.A_p.resize(3 * MPC_HORIZON, mpc.C.cols());
    mpc.A_v.resize(3 * MPC_HORIZON, mpc.C.cols());
    mpc.A_a.resize(3 * MPC_HORIZON, mpc.C.cols());
    mpc.M_p.resize(3 * MPC_HORIZON, mpc.M.cols());
    mpc.M_v.resize(3 * MPC_HORIZON, mpc.M.cols());
    mpc.M_a.resize(3 * MPC_HORIZON, mpc.M.cols());
    mpc.B_a.resize(mpc.M_p.rows(), 1);
    mpc.B_v.resize(mpc.M_v.rows(), 1);
    mpc.B_p.resize(mpc.M_a.rows(), 1);

    for (int i = 0; i < MPC_HORIZON; i++) {
        mpc.A_p.block(3 * i, 0, 3, mpc.A_p.cols()) = mpc.C.block(9 * i + 0, 0, 3, mpc.C.cols());
        mpc.A_v.block(3 * i, 0, 3, mpc.A_v.cols()) = mpc.C.block(9 * i + 3, 0, 3, mpc.C.cols());
        mpc.A_a.block(3 * i, 0, 3, mpc.A_a.cols()) = mpc.C.block(9 * i + 6, 0, 3, mpc.C.cols());
        mpc.M_p.block(3 * i, 0, 3, mpc.M_p.cols()) = mpc.M.block(9 * i + 0, 0, 3, mpc.M.cols());
        mpc.M_v.block(3 * i, 0, 3, mpc.M_v.cols()) = mpc.M.block(9 * i + 3, 0, 3, mpc.M.cols());
        mpc.M_a.block(3 * i, 0, 3, mpc.M_a.cols()) = mpc.M.block(9 * i + 6, 0, 3, mpc.M.cols());
    }

    mpc.A_sys.block(0, 0, A_u.rows(), A_u.cols()) = A_u;
    mpc.A_sys.block(A_u.rows(), 0, mpc.A_a.rows(), mpc.A_a.cols()) = mpc.A_a;
    mpc.A_sys.block(A_u.rows() + mpc.A_a.rows(), 0, mpc.A_v.rows(), mpc.A_v.cols()) = mpc.A_v;

    // std::cout << "mpc.A_a: " << std::endl << mpc.A_a << std::endl;
    // std::cout << "mpc.A_v: " << std::endl << mpc.A_v << std::endl;
    // std::cout << "mpc.A_p: " << std::endl << mpc.A_p << std::endl;

    Eigen::VectorXd x_0;
    x_0.resize(9, 1);
    x_0.setZero();
    UpdateBound(mpc, x_0);
    // x_0 << 1.0,0.5,0, 2,1,0.5, 0.1,0,0;
    // UpdateBound(mpc, x_0);
}

void MPCPlannerClass::UpdateBound(mpc_osqp_t& mpc, const Eigen::VectorXd& x_0)
{
    if (x_0.rows() != mpc.M_p.cols() || x_0.rows() != mpc.M_v.cols() || x_0.rows() != mpc.M_a.cols()) {
        ROS_ERROR("[MPC]: Update bound error!");
        return;
    }
    mpc.B_p = mpc.M_p * x_0;
    mpc.B_v = mpc.M_v * x_0;
    mpc.B_a = mpc.M_a * x_0;

    mpc.A_sys_low.block(0, 0, mpc.u_low.rows(), 1) = mpc.u_low;
    mpc.A_sys_upp.block(0, 0, mpc.u_upp.rows(), 1) = mpc.u_upp;
    mpc.A_sys_low.block(mpc.u_low.rows(), 0, mpc.a_low.rows(), 1) = mpc.a_low - mpc.B_a;
    mpc.A_sys_upp.block(mpc.u_upp.rows(), 0, mpc.a_upp.rows(), 1) = mpc.a_upp - mpc.B_a;
    mpc.A_sys_low.block(mpc.u_low.rows() + mpc.a_low.rows(), 0, mpc.v_low.rows(), 1) = mpc.v_low - mpc.B_v;
    mpc.A_sys_upp.block(mpc.u_upp.rows() + mpc.a_upp.rows(), 0, mpc.v_upp.rows(), 1) = mpc.v_upp - mpc.B_v;

    // std::cout << "mpc.B_a: " << mpc.B_a.transpose() << std::endl;
    // std::cout << "mpc.B_v: " << mpc.B_v.transpose() << std::endl;
    // std::cout << "mpc.B_p: " << mpc.B_p.transpose() << std::endl;
}

