
#include "depth_gene/depth_grid.hpp"

int main(int argc, char **argv) {
    ros::init(argc, argv, "local_depth_node");
    ros::NodeHandle n("~");

    DepthGrid dg(n);

    ros::AsyncSpinner spinner(0);
    spinner.start();

    ros::waitForShutdown();
}
