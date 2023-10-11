#include "perfect_drone_sim/perfect_drone_model.hpp"

#define BACKWARD_HAS_DW 1
#include "backward.hpp"
namespace backward {
backward::SignalHandling sh;
}


int main(int argc, char **argv) {
  ros::init(argc, argv, "perfect_tracking");
  ros::NodeHandle n("~");


  PerfectDrone dp(n);


  ros::AsyncSpinner spinner(0);
  spinner.start();

  ros::waitForShutdown();
  return 0;
}

