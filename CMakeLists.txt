cmake_minimum_required(VERSION 3.0.2)
project(parkerstage)

# Data-only package: URDF, meshes and launch files.  Works with catkin_make
# and with `colcon build` (catkin build type) in ROS 1 workspaces.
find_package(catkin REQUIRED)

catkin_package()

# Runtime assets.  The launch file resolves the URDF via $(find parkerstage)
# and the URDF references meshes via package://parkerstage/meshes/..., so
# both the urdf and meshes directories must end up in the package share dir.
install(DIRECTORY launch urdf meshes
        DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
