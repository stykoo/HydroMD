cmake_minimum_required(VERSION 2.8)

PROJECT(HydroMD)

set(EXECUTABLE_NAME "HydroMD")

# Build in debug mode by default
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Debug)
endif()
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "-W -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-Ofast -march=native -g -pg") # Profiling
set(CMAKE_CXX_STANDARD 14)

# Packages needed
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
# find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Threads REQUIRED)
#include_directories(${HDF5_INCLUDE_DIRS})

find_package(MKL)
include_directories(${MKL_INCLUDE_DIRS})

file(
	GLOB
	source_files
	src/*
)

add_executable(
	${EXECUTABLE_NAME}
	${source_files}
)

if(MKL_FOUND)
	target_compile_definitions(${EXECUTABLE_NAME} PRIVATE USE_MKL)
	target_link_libraries(${EXECUTABLE_NAME} -Wl,--start-group ${MKL_LIBRARIES} -Wl,--end-group pthread dl)
endif()


target_link_libraries(
	${EXECUTABLE_NAME}
	#${HDF5_LIBRARIES}
	${Boost_LIBRARIES}
)
