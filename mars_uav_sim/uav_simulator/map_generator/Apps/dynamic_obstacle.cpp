#include "random_map_generator/random_map_generator.hpp"

int main(int argc, char **argv) {
  ros::init(argc, argv, "random_complex_scene");
  ros::NodeHandle n("~");

  ros::Publisher dyn_cloud_pub = n.advertise<sensor_msgs::PointCloud2>("/dyn_cloud", 10000);
  pcl::PointCloud<pcl::PointXYZ> pc_to_pub;
  Eigen::Matrix3Xd points(3, 1000), temp_pt(3, 1000);
  for (int x = 0; x < 10; x++) {
	for (int y = 0; y < 10; y++) {
	  for (int z = 0; z < 10; z++) {
		points.col(x * 100 + y * 10 + z) = Eigen::Vector3d(x, y, z)*0.1;
	  }
	}
  }
  Eigen::Vector3d center_pos;
  ros::Time start_t = ros::Time::now();
  sensor_msgs::PointCloud2 pc2;
  while (ros::ok()) {
	pc_to_pub.clear();
	double t = (start_t - ros::Time::now()).toSec();
	Eigen::Vector3d offset(3 * sin(t), 3 * cos(t), 1 * sin(t) + 1);
	temp_pt = points.colwise() + offset;
	for (int i = 0; i < temp_pt.cols(); i++) {
	  pc_to_pub.push_back(pcl::PointXYZ(temp_pt.col(i).x(), temp_pt.col(i).y(), temp_pt.col(i).z()));
	}

	pcl::toROSMsg(pc_to_pub, pc2);
	pc2.header.frame_id = "world";
	pc2.header.stamp = ros::Time::now();
	dyn_cloud_pub.publish(pc2);
	ros::Duration(0.1).sleep();
  }

  ros::AsyncSpinner spinner(0);
  spinner.start();

  ros::waitForShutdown();
  return 0;
}
