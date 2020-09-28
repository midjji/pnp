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
