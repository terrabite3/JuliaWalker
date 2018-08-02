cmake_minimum_required(VERSION 3.7)

project(JuliaWalker)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

add_executable(${PROJECT_NAME} 
    test.cpp
    )

target_link_libraries(${PROJECT_NAME}
    CONAN_PKG::libpng
)