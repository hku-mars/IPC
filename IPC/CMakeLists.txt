cmake_minimum_required(VERSION 3.0.2)
project(ipc)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")
# set(CMAKE_BUILD_TYPE "Debug")

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)

find_package(OsqpEigen REQUIRED)

find_package(fmt REQUIRED)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  nav_msgs
  std_msgs
  quadrotor_msgs
  eigen_conversions
  visualization_msgs
  sensor_msgs
  pcl_ros
  mavros
)

catkin_package(

)

include_directories(
  include 
  ${catkin_INCLUDE_DIRS}
  ${Eigen3_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
)

add_executable(ipc_node
  src/main.cpp
  src/planner.cpp
  include/astar.cpp
  include/mpc.cpp
  include/local_astar.cpp
)
add_dependencies(ipc_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(ipc_node
  ${catkin_LIBRARIES}
  ${EIGEN3_LIBRARIES}
  ${PCL_LIBRARIES}
  -ldw
  tbb
  OsqpEigen::OsqpEigen
  fmt::fmt
)

add_executable(fast_avoid
  test/fast_avoid.cpp
)
add_dependencies(fast_avoid ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(fast_avoid
  ${catkin_LIBRARIES}
  ${EIGEN3_LIBRARIES}
)


