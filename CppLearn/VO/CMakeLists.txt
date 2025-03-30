cmake_minimum_required(VERSION 3.16)

project(vo)
set(CMAKE_CXX_STANDARD 17)


find_package(Ceres REQUIRED)
find_package(OpenCV 4 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(g2o REQUIRED)
find_package(Sophus REQUIRED)

include_directories(
    &{Ceres_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
    ${Sophus_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${G2O_INCLUDE_DIR}
    )
    

add_executable(gaussian_newton gaussian_newton.cpp)
target_link_libraries(gaussian_newton
    ${OpenCV_LIBS}
    )

add_executable(curvefitting_ceres curvefitting_ceres.cpp)
target_link_libraries(curvefitting_ceres
    ${OpenCV_LIBS}
    ${CERES_LIBRARIES}
    )

add_executable(curvefitting_g2o curvefitting_g2o.cpp)
target_link_libraries(curvefitting_g2o
    ${OpenCV_LIBS}
    g2o_core g2o_stuff
    )
