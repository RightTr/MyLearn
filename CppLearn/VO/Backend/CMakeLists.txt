cmake_minimum_required(VERSION 3.16)

project(backend)
set(CMAKE_CXX_FLAGS "-O3")


find_package(Ceres REQUIRED)
find_package(OpenCV 4 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(g2o REQUIRED)
find_package(fmt REQUIRED)

include_directories(
    &{Ceres_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${G2O_INCLUDE_DIR}
    )

SET(G2O_LIBS g2o_csparse_extension g2o_stuff g2o_core g2o_solver_csparse)

add_library(bal BAL.cpp)

add_executable(BA_ceres BA_ceres.cpp)
target_link_libraries(BA_ceres
    ${CERES_LIBRARIES}
    bal)

add_executable(BA_g2o BA_g2o.cpp)
target_link_libraries(BA_g2o
    ${CERES_LIBRARIES}
    ${G2O_LIBS}
    fmt::fmt
    bal)