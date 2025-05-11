cmake_minimum_required(VERSION 3.16)

project(tcp)

add_executable(tcp_server tcp_srv.cpp)
add_executable(tcp_client tcp_cli.cpp)
add_executable(tcp_process tcp_process.cpp)
