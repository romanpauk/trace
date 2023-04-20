cmake_minimum_required(VERSION 3.22)
project(trace)
include(FetchContent)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS NO)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
enable_testing()

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache) 
endif(CCACHE_FOUND)

option(TRACE_ENABLE_TESTING "enable testing using googletest" ON)
option(TRACE_ENABLE_BENCHMARK "enable benchmark using googlebenchmark" ON)

if(TRACE_ENABLE_TESTING)
    FetchContent_Declare(googletest GIT_REPOSITORY https://github.com/google/googletest.git GIT_TAG main)
    if(WIN32)
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    endif()
    FetchContent_MakeAvailable(googletest)
endif()

if(TRACE_ENABLE_BENCHMARK)
    FetchContent_Declare(googlebenchmark GIT_REPOSITORY https://github.com/google/benchmark.git GIT_TAG main) 
    set(BENCHMARK_ENABLE_TESTING off)
    FetchContent_MakeAvailable(googlebenchmark)
endif()

add_library(trace INTERFACE)
target_include_directories(trace INTERFACE include)

target_sources(trace INTERFACE
    include/trace/trace.h
)

if(TRACE_ENABLE_TESTING)
    add_executable(trace_test
        test/trace.cpp
    )

    add_test(trace_test COMMAND trace_test)
    target_link_libraries(trace_test trace gtest_main)
    target_include_directories(trace_test PRIVATE test)
endif()

if(TRACE_ENABLE_BENCHMARK)
    add_executable(trace_benchmark
        benchmark/trace.cpp
    )

    target_link_libraries(trace_benchmark trace benchmark::benchmark)
    target_include_directories(trace_benchmark PRIVATE benchmark)
endif()
