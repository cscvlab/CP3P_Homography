cmake_minimum_required(VERSION 3.22)
project(CP3P_Homography_IJCV)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE " -O3 -march=native -mtune=native -DNDEBUG")

#Change to your Eigen package address
include_directories("./eigen-3.4.0")
include_directories("./colmap_util")


set(epnp-srcs "")
list(APPEND epnp-srcs EPNP/EPNP.h)
list(APPEND epnp-srcs EPNP/EPNP.cpp)

set(colmap_util-src "")
list(APPEND colmap_util-src colmap_util/eigen_alignment.h)
list(APPEND colmap_util-src colmap_util/models.h)
list(APPEND colmap_util-src colmap_util/models.cc)
list(APPEND colmap_util-src colmap_util/types.h)

set(util-src "")
list(APPEND util-src utils/quarticSolver.h)


set(cost_function-src "")
list(APPEND cost_function-src cost_function/cost_function.h)

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

add_executable(CP3P_Homography_IJCV
        P3P_HD.h
        demo_P3P.cpp
        demo_Ransac.cpp
        main.cpp
        demo_Ransac.h
        demo_P3P.h
        ${epnp-srcs}
        ${colmap_util-src}
        ${cost_function-src}
        ${util-src}
)
target_link_libraries(CP3P_Homography_IJCV ${CERES_LIBRARIES})
