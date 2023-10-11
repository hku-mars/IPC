#include "random_map_generator/random_map_generator.hpp"

int main (int argc, char** argv)
{
    ros::init (argc, argv, "random_complex_scene");
    ros::NodeHandle n( "~" );

    RandomMap rdm(n);

    ros::AsyncSpinner spinner(0);
    spinner.start();

    ros::waitForShutdown();
    return 0;
}
