cmake_minimum_required(VERSION 3.0)

project(a1)
set(CMAKE_CXX_STANDARD 11)

find_package(GLEW REQUIRED)
find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS})

add_library(col781 src/hw.cpp src/camera.cpp)
target_link_libraries(col781 GLEW::GLEW glm::glm OpenGL::GL ${SDL2_LIBRARIES})

add_executable(example src/example.cpp)
target_link_libraries(example col781)
