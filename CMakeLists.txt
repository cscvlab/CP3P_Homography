cmake_minimum_required(VERSION 3.22)
project(CP3P_Homography)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE " -O3 -march=native -mtune=native -DNDEBUG")

#Change to your Eigen package address
include_directories("./eigen-3.4.0")
include_directories("utils/colmap_util")

set(util-src "")
list(APPEND util-src utils/functionSolver.h)

list(APPEND util-src utils/colmap_util/eigen_alignment.h)
list(APPEND util-src utils/colmap_util/models.h)
list(APPEND util-src utils/colmap_util/models.cc)
list(APPEND util-src utils/colmap_util/types.h)

list(APPEND util-src utils/cost_function/cost_function.h)

list(APPEND util-src utils/EPNP/EPNP.h)
list(APPEND util-src utils/EPNP/EPNP.cpp)

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

add_executable(CP3P_Homography
        P3P_HD.h
        demo_P3P.cpp
        demo_Ransac.cpp
        main.cpp
        demo_Ransac.h
        demo_P3P.h
        ${util-src}
)
target_link_libraries(CP3P_Homography ${CERES_LIBRARIES})
