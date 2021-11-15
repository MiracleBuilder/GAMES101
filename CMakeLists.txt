cmake_minimum_required (VERSION 3.1.11)
project (Transformation)

INCLUDE_DIRECTORIES(
    D:/ProgramFiles/vcpkg/installed/x86-windows/include/
)

LINK_DIRECTORIES(
    D:/ProgramFiles/vcpkg/installed/x86-windows/lib/
)

# each target need
add_executable (Transformation src/assignment/assignment0/Assignment0.cpp)
add_executable(Rasterizer src/assignment/assignment1/Assignment1.cpp 
                    src/assignment/assignment1/rasterizer.hpp 
                    src/assignment/assignment1/Triangle.hpp 
                    src/assignment/assignment1/rasterizer.cpp 
                    src/assignment/assignment1/Triangle.cpp)


find_package(Eigen3 REQUIRED)

# import opencv
set(OpenCV_INCLUDE_DIRS D:/ProgramFiles/vcpkg/installed/x86-windows/include)
set(OpenCV_LIBS D:/ProgramFiles/vcpkg/installed/x86-windows/lib)
find_package(OpenCV REQUIRED)
target_include_directories(Rasterizer PRIVATE ${OpenCV_INCLUDE_DIRS})
target_link_libraries(Rasterizer PRIVATE  ${OpenCV_LIBS})


