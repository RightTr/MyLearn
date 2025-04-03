cmake_minimum_required(VERSION 3.16)

project(backend)
set(CMAKE_CXX_FLAGS "-O3")


find_package(Ceres REQUIRED)
find_package(OpenCV 4 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(g2o REQUIRED)

include_directories(
    &{Ceres_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${G2O_INCLUDE_DIR}
    )

add_library(bal BAL.cpp)

add_executable(BA_ceres BA_ceres.cpp)
target_link_libraries(BA_ceres
    ${CERES_LIBRARIES}
    bal)