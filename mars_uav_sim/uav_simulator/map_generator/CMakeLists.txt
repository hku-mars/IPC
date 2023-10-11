cmake_minimum_required(VERSION 2.8.3)
project(map_generator)

add_compile_options(-std=c++14)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g -fPIC")

find_package(catkin REQUIRED COMPONENTS
        roscpp
        std_msgs
        pcl_ros
        geometry_msgs
        nav_msgs
        message_generation
        )
find_package(fmt)
find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
find_package(OpenCV REQUIRED)
catkin_package(
        INCLUDE_DIRS include
        LIBRARIES random_map_generator
        CATKIN_DEPENDS roscpp rospy std_msgs
        #  DEPENDS system_lib
)

include_directories(
        SYSTEM
        include
        ../../../../devel/include
        ${catkin_INCLUDE_DIRS}
        ${Eigen3_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
)
link_directories(${PCL_LIBRARY_DIRS})

set(THIRD_PARTY
        ${catkin_LIBRARIES}
        ${PCL_LIBRARIES}
        fmt::fmt
        )

add_executable(random_map_node
        Apps/random_map_gene_node.cpp
        )
target_link_libraries(random_map_node
        ${THIRD_PARTY}
        )

add_executable(racing_map_gene_node
        Apps/racing_map_gene_node.cpp
        )
target_link_libraries(racing_map_gene_node
        ${THIRD_PARTY}
        )

add_executable(color_map_for_pointcloud
        Apps/color_map_for_pointcloud.cpp
        )
target_link_libraries(color_map_for_pointcloud
        ${THIRD_PARTY}
        )

add_executable(dynamic_obstacle
        Apps/dynamic_obstacle.cpp
        )
target_link_libraries(dynamic_obstacle
        ${THIRD_PARTY}
        )

add_executable(map_recorder src/map_recorder.cpp)
target_link_libraries(map_recorder
        ${catkin_LIBRARIES}
        ${PCL_LIBRARIES})

add_executable(map_pub src/map_publisher.cpp)
target_link_libraries(map_pub
        ${catkin_LIBRARIES}
        ${PCL_LIBRARIES})

