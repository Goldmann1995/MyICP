cmake_minimum_required(VERSION 2.8.3)
project(calibration_rot_axis)

# add_compile_options(-std=c++11)
set(CMAKE_CXX_STANDARD 14)

find_package(catkin REQUIRED COMPONENTS
pcl_ros
roscpp
sensor_msgs
livox_ros_driver
)

find_package(PCL 1.8 REQUIRED)
find_package(Eigen3  REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS roscpp sensor_msgs pcl_ros 
)

include_directories(
 include
 ${catkin_INCLUDE_DIRS}
 ${EIGEN3_INCLUDE_DIR}
 ${PCL_INCLUDE_DIRS}
)



link_directories(${PCL_LIBRARY_DIRS})

add_executable(${PROJECT_NAME}_node src/main.cpp src/preprocess.cpp)


target_link_libraries(${PROJECT_NAME}_node
  ${catkin_LIBRARIES}
  ${PCL_LIBRARIES}
  -lgomp
)
