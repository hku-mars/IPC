cmake_minimum_required(VERSION 2.8.3)
project(perfect_drone_sim)

add_compile_options(-std=c++14)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g -fPIC")
add_compile_options(-DQT_NO_VERSION_TAGGING)
add_definitions(-DROOT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\")
find_package(catkin REQUIRED COMPONENTS
        roscpp
        std_msgs
        pcl_ros
        geometry_msgs
        nav_msgs
        message_generation
        quadrotor_msgs
        visualization_msgs
        message_filters
        )
find_package(PythonLibs REQUIRED)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
find_package(OpenCV REQUIRED)
find_package(fmt)
catkin_package(
        INCLUDE_DIRS include
        LIBRARIES
        CATKIN_DEPENDS roscpp rospy std_msgs
        #  DEPENDS system_lib
)

include_directories(
        SYSTEM
        include
        ../../../devel/include
        /usr/include/x86_64-linux-gnu/qt5
        /usr/include/x86_64-linux-gnu/qt5/QtGui
        /usr/include/x86_64-linux-gnu/qt5/QtCore
        ${catkin_INCLUDE_DIRS}
        ${Eigen3_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
        ${OpenCV_INCLUDE_DIRS}
        ${PYTHON_INCLUDE_DIRS}
)

set(THIRD_PARTY
        ${catkin_LIBRARIES}
        ${PCL_LIBRARIES}
        ${OpenCV_LIBRARIES}
        tbb
        fmt::fmt
        -ldw
        )


add_executable(perfect_tracker
        src/perfect_drone.cpp
        )
target_link_libraries(perfect_tracker
        ${THIRD_PARTY}
        )