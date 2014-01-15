cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(TSDFOctree)

file(GLOB srcs
   "cpp/*.cpp"
)

file(GLOB hdrs
	"cpp/*.h" "cpp/*.hpp"
)

set(OpenCV_DIR "F:/Libraries/opencv/build" CACHE PATH "description")
find_package(OpenCV REQUIRED)

set(OctomapIncludeDirectory "F:/Libraries/octomap-1.6.1/octomap/include" CACHE PATH "description")
set(OctomapLibDirectory "F:/Libraries/octomap-1.6.1/lib/Release" CACHE PATH "description")
set(EIGEN_INCLUDE_DIRS "F:/Libraries/Eigen" CACHE PATH "description")

set(OctomapLibraries ${OctomapLibDirectory}/octomap.lib ${OctomapLibDirectory}/octomath.lib)

include_directories(${OctomapIncludeDirectory})
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})

link_directories(${OpenCV_LIBRARY_DIRS})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "./../bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "./../bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE "./../bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG "./../bin")
set(CMAKE_DEBUG_POSTFIX "d")

set(target TSDFOctree)
add_executable(${target} ${srcs}  ${hdrs})
target_link_libraries(${target} ${OctomapLibraries} ${OpenCV_LIBRARIES})