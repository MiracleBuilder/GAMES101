cmake_minimum_required (VERSION 3.1.11)
project (Transformation)

INCLUDE_DIRECTORIES(
    D:/ProgramFiles/vcpkg/installed/x86-windows/include/
)

LINK_DIRECTORIES(
    D:/ProgramFiles/vcpkg/installed/x86-windows/lib/
)

# each target need
add_executable (Assignment0 src/assignment/assignment0/Assignment0.cpp)
add_executable(Assignment1 src/assignment/assignment1/Assignment1.cpp 
                    src/assignment/assignment1/rasterizer.hpp 
                    src/assignment/assignment1/Triangle.hpp 
                    src/assignment/assignment1/rasterizer.cpp 
                    src/assignment/assignment1/Triangle.cpp)
add_executable(Assignment2 src/assignment/assignment2/Assignment2.cpp 
                    src/assignment/assignment2/rasterizer.hpp 
                    src/assignment/assignment2/rasterizer.cpp 
                    src/assignment/assignment2/global.hpp 
                    src/assignment/assignment2/Triangle.hpp 
                    src/assignment/assignment2/Triangle.cpp)


find_package(Eigen3 REQUIRED)

# import opencv
set(OpenCV_INCLUDE_DIRS D:/ProgramFiles/vcpkg/installed/x86-windows/include)
set(OpenCV_LIBS D:/ProgramFiles/vcpkg/installed/x86-windows/lib)
find_package(OpenCV REQUIRED)
target_include_directories(Assignment1 PRIVATE ${OpenCV_INCLUDE_DIRS})
target_link_libraries(Assignment1 PRIVATE  ${OpenCV_LIBS})

target_include_directories(Assignment2 PRIVATE ${OpenCV_INCLUDE_DIRS})
target_link_libraries(Assignment2 PRIVATE  ${OpenCV_LIBS})


