project(DFVS)

cmake_minimum_required(VERSION 3.5.1)

# set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
add_compile_options(
	-Wall
	-Wextra
	-Werror
)
if("cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
	add_compile_options(-std=c++20)
else()
	add_compile_options(-std=c++17)
endif()
add_compile_definitions(
	$<$<CONFIG:DEBUG>:PROFILING>
	$<$<CONFIG:RELWITHDEBINFO>:PROFILING>
)
file(GLOB SOURCES "src/*.cpp")


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
add_executable(${PROJECT_NAME} ${SOURCES})
if(${CMAKE_BUILD_TYPE} STREQUAL Release)
	message("Static compilation")
	target_link_libraries(${PROJECT_NAME} -static)
else()
	message("Non static compilation")
endif()
