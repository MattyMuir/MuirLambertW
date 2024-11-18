# EXECUTABLE PROJECT - tests

# === Includes ===
cmake_policy(SET CMP0069 NEW)
include(../macros.cmake)

# === Create Executable ===
add_executable(tests "tests.cpp" "types.h" "ulp.cpp" "ulp.h" "wrappers.h"  "ReciprocalDistribution.h"  "ReciprocalDistributionEx.h" "../bench/others/BarryLambertW.cpp" "../bench/others/FukushimaLambertW.cc" "histogram.h" "histogram.cpp" "constants.h")

# === Libraries ===
target_link_libraries(tests PUBLIC MuirLambertW)
target_include_directories(tests PUBLIC "../include/")

find_package(PkgConfig)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
target_link_libraries(tests PRIVATE PkgConfig::mpfr)

find_package(sleef REQUIRED)
target_link_libraries(tests PUBLIC sleef::sleef)
target_include_directories(tests PRIVATE ${SLEEF_INCLUDE_DIR})

find_package(ReferenceLambertW)
target_link_libraries(tests PRIVATE ReferenceLambertW::ReferenceLambertW)

find_package(boost_math REQUIRED CONFIG)

target_include_directories(tests PRIVATE Boost::math)

# === Feature Enables ===
if (MUIRW_MSVC_STATIC_RUNTIME)
    set_property(TARGET tests PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_compile_features(tests PUBLIC cxx_std_20)
enable_ipo(tests)
set_arch(tests)