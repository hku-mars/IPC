#include "quadrotor_msgs/PositionCommand.h"
#include <ros/ros.h>
#include "geometry_msgs/PoseStamped.h"

class Interface
{
public:
    Interface();

private:
    ros::NodeHandle nh;
    ros::Publisher pub,new_pub;
    ros::Subscriber sub;

    quadrotor_msgs::PositionCommand cmd;
    geometry_msgs::PoseStamped pos_cmd;
    int _n_seq;

    void messageCallback(const geometry_msgs::PoseStampedConstPtr &msg);
};

Interface::Interface() {
    pub = nh.advertise<quadrotor_msgs::PositionCommand>
            ("/planning/pos_cmd", 10);
    sub = nh.subscribe<geometry_msgs::PoseStamped>
            ("/goal", 1, &Interface::messageCallback, this);
            ///move_base_simple/goal.
    new_pub = nh.advertise<geometry_msgs::PoseStamped>
            ("/planning/pos_new_cmd", 10);

    _n_seq = 0;

    /* kP */
    double pos_gain[3] = { 5.7, 5.7, 6.2 };
    double vel_gain[3] = { 3.4, 3.4, 4.0 };

    /* control parameter */
    cmd.kx[0] = pos_gain[0];
    cmd.kx[1] = pos_gain[1];
    cmd.kx[2] = pos_gain[2];

    cmd.kv[0] = vel_gain[0];
    cmd.kv[1] = vel_gain[1];
    cmd.kv[2] = vel_gain[2];

    ros::Time start_t = ros::Time::now();
    while(ros::ok()){
        ros::Time cur_t = ros::Time::now();
        if((cur_t - start_t).toSec()>1){
            break;
        }
        cmd.header.stamp = cur_t;
        cmd.header.frame_id = "world";

        cmd.trajectory_id = 0;
        cmd.trajectory_flag = 1;
        cmd.position.x = 0;
        cmd.position.y = 0;
        cmd.position.z = 2;
        cmd.velocity.x = 0;
        cmd.velocity.y = 0;
        cmd.velocity.z = 0;
        cmd.acceleration.x = 0;
        cmd.acceleration.y = 0;
        cmd.acceleration.z = 0;

        pub.publish(cmd);
        ros::Duration(0.1).sleep();
    }

    ros::spin();
}

void Interface::messageCallback(const geometry_msgs::PoseStampedConstPtr &msg) {
    // header
    cmd.header.stamp = msg->header.stamp;
    cmd.header.frame_id = "world";

    cmd.trajectory_id = 0;
    cmd.trajectory_flag = 1;
    cmd.position.x = msg->pose.position.x;
    cmd.position.y = msg->pose.position.y;
    cmd.position.z = msg->pose.position.z;
    cmd.velocity.x = 0;
    cmd.velocity.y = 0;
    cmd.velocity.z = 0;
    cmd.acceleration.x = 0;
    cmd.acceleration.y = 0;
    cmd.acceleration.z = 0;

    pub.publish(cmd);

    pos_cmd.header.stamp = msg->header.stamp;
    pos_cmd.header.frame_id = "world";
    pos_cmd.pose.position.x = msg->pose.position.x;
    pos_cmd.pose.position.y = msg->pose.position.y;
    pos_cmd.pose.position.z = msg->pose.position.z;
    // pos_cmd.pose.orientation.z = msg->pose.orientation.z;

    new_pub.publish(pos_cmd);

    ROS_INFO("Goal set %f, %f, %f", msg->pose.position.x, msg->pose.position.y, msg->pose.position.z);

}

int main(int argc, char** argv)
{
    ROS_WARN("*****START*****");
    ros::init(argc, argv, "test_interface");
    Interface Int;

    return 0;
}
