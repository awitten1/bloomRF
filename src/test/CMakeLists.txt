include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  test_bloomrf
  test_bloomrf.cpp
  test_bloomrf_signed.cpp
  test_bloomrf_floats.cpp
)
target_link_libraries(
  test_bloomrf
  bloomRF
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_bloomrf)

add_executable(
  experiments
  experiments.cpp
)

target_link_libraries(
  experiments
  bloomRF
)
