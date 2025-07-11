cmake_minimum_required(VERSION 3.10)
project(HaploGI VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Add executable
add_executable(HaploGI src/HaploGI.cpp)

# Add this to enable 'make install'
install(TARGETS HaploGI
    RUNTIME DESTINATION bin
)

# Optional: Add compiler warnings (recommended)
if(MSVC)
    target_compile_options(HaploGI PRIVATE /W4 /permissive-)
else()
    target_compile_options(HaploGI PRIVATE -Wall -Wextra -pedantic)
endif()

