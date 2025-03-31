cmake_minimum_required(VERSION 3.5.1)

project(DFVS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU OR CMAKE_CXX_COMPILER_ID MATCHES Clang)
	set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -DPROFILING")
	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DPROFILING")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
	add_compile_options(
		-Wall
		-Wextra
		-Werror
		-fdiagnostics-color
	)
endif()

file(GLOB SOURCES "src/*.cpp")


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
add_executable(${PROJECT_NAME} ${SOURCES})
if(${CMAKE_BUILD_TYPE} STREQUAL Release)
	message("Static compilation")
	target_link_libraries(${PROJECT_NAME} -static)
else()
	message("Non static compilation")
endif()
