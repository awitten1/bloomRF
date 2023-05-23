option (SANITIZE "Enable one of the code sanitizers" "")

set (SAN_FLAGS "${SAN_FLAGS} -g -fno-omit-frame-pointer")

if (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
    message(STATUS "Detected Apple Clang compiler.")
    set (COMPILER_CLANG 1)
endif()


if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Detected Clang compiler.")
    set (COMPILER_CLANG 1)
endif()

if (SANITIZE AND COMPILER_CLANG)
    if (SANITIZE STREQUAL "undefined")
        set (UBSAN_FLAGS "-fsanitize=undefined")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SAN_FLAGS} ${UBSAN_FLAGS}")
    elseif(SANITIZE STREQUAL "memory")
        set (MSAN_FLAGS "-fsanitize=memory")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SAN_FLAGS} ${MSAN_FLAGS}")
    else ()
        message (FATAL_ERROR "Unknown sanitizer type: ${SANITIZE}")
    endif ()
endif()