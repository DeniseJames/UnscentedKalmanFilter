cmake_minimum_required(VERSION 3.5)
project(Unscented_Kalman_Filter)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp
         ukf.cpp ukf.h measurement_package.h ground_truth_package.h tools.cpp tools.h)

add_executable(Unscented_Kalman_Filter ${SOURCE_FILES} ground_truth_package.h tools.cpp tools.h)

file(COPY ${CMAKE_SOURCE_DIR}/sample-laser-radar-measurement-data-1.txt DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/sample-laser-radar-measurement-data-2.txt DESTINATION ${CMAKE_BINARY_DIR})