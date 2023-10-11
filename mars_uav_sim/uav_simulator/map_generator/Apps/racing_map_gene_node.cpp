
#include "racing_map_generator/racing_map_generator.hpp"

int main (int argc, char** argv)
{
    ros::init (argc, argv, "racing_map_gene_node");
    ros::NodeHandle n( "~" );

    racing_map_generator::RacingMapGenerator rdm(n);

    ros::AsyncSpinner spinner(0);
    spinner.start();

    ros::waitForShutdown();
    return 0;
}
