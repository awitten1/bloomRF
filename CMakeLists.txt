cmake_minimum_required(VERSION 3.22)
project("BloomRF filters")
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(cmake/sanitize.cmake)

add_subdirectory(src)
