CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
PROJECT(lambdatwist-pnp CXX) # defines ${CMAKE_PROJECT_NAME}(top level) PROJECT_NAME(most recent)
add_subdirectory(pybind11)



macro(print_list name list)
    message("\n Listing: ${name}\n")
    foreach(item IN LISTS ${list})
        message("     ${name}: ${item}")
    endforeach()
    message("\n Listing: ${name} - done \n")
endmacro()

if(NOT CMAKE_BUILD_TYPE)
    MESSAGE("-- No build type specified; defaulting to CMAKE_BUILD_TYPE=Release.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: Debug Release"  FORCE)
else()
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        MESSAGE("\n${line}")
        MESSAGE("\n-- Build type: Debug. Performance will be terrible!")
        MESSAGE("-- Add -DCMAKE_BUILD_TYPE=Release to the CMake command line to get an optimized build.")
        MESSAGE("\n${line}")
    endif()
endif()
# add warnings
#set(WARN "-Wall -Wextra -Wpedantic -Werror")
set(CMAKE_CXX_FLAGS_DEBUG " -fno-omit-frame-pointer -g -pg  -rdynamic ${WARN} ") # dynamic is for the improved asserts
set(CMAKE_CXX_FLAGS_RELEASE " -O3 -march=native -DNDEBUG ${WARN} ")






option(RANDOM_SEED_FROM_TIME "generate a seed from time, default" OFF)
if(RANDOM_SEED_FROM_TIME)
    add_definitions(-DRANDOM_SEED_FROM_TIME)
else()
    add_definitions(-DRANDOM_SEED_VALUE=${RANDOM_SEED})
endif()

add_definitions("-std=c++11")
include_directories(".")
include_directories("lambdatwist-p3p/utils")

find_package(Ceres REQUIRED )
INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

#including headers here make them show up in older ides

set(lambdatwist-srcs "")
list(APPEND lambdatwist-srcs lambdatwist/p3p_timers.h)
list(APPEND lambdatwist-srcs lambdatwist/lambdatwist.p3p.h)
list(APPEND lambdatwist-srcs lambdatwist/solve_cubic.h)
list(APPEND lambdatwist-srcs lambdatwist/solve_eig0.h)
list(APPEND lambdatwist-srcs lambdatwist/refine_lambda.h)


set(util-srcs "")

list(APPEND util-srcs utils/mlibtime.cpp)

list(APPEND util-srcs utils/mlibtime.h )
list(APPEND util-srcs utils/matlab_helpers.h)
list(APPEND util-srcs utils/random.h)
list(APPEND util-srcs utils/random_vectors.h)

list(APPEND util-srcs utils/string_helpers.h)
list(APPEND util-srcs utils/vector_utils.h)
list(APPEND util-srcs utils/cvl/matrix.h)
list(APPEND util-srcs utils/cvl/pose.h)
list(APPEND util-srcs utils/cvl/rotation_helpers.h)

list(APPEND util-srcs pnp_python_binding.h)

#list(APPEND util-srcs pnp.h)


print_list("utils" util-srcs)
#file(GLOB_RECURSE utils2 )


#list(APPEND util-srcs data.h)
#list(APPEND util-srcs solver.h)
#list(APPEND util-srcs mnumerics.h)
#list(APPEND util-srcs p3p_generator.h)

add_library(libpnp SHARED ${util-srcs} ${lambdatwist-srcs} pnp_ransac.cpp p4p.cpp)
target_link_libraries(libpnp ${CERES_LIBRARIES})

# headers are included for ide benefit...
add_executable(main test_pnp.cpp)
target_link_libraries(main libpnp)



# Build pybind11 bindinis
pybind11_add_module(pnp_python_binding pnp_python_binding.cpp)
target_link_libraries(pnp_python_binding PUBLIC libpnp ${CERES_LIBRARIES})
# OBS! Pybind11 require that you specify python exectuable to build for, e.g. cmake .. -DPYTHON_EXECUTABLE=python3


message("================================================================================")
