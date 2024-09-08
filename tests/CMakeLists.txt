# EXECUTABLE PROJECT - tests

# === Includes ===
cmake_policy(SET CMP0069 NEW)
include(../macros.cmake)

# === Create Executable ===
add_executable(tests "tests.cpp" "reference/ReferenceLambertW.cpp" "reference/ReferenceLambertW.h" "types.h" "ulp.cpp" "ulp.h" "wrappers.h" "ReciprocalDistribution.cpp" "ReciprocalDistribution.h" "ReciprocalDistributionEx.cpp" "ReciprocalDistributionEx.h" "../bench/others/BarryLambertW.cpp" "../bench/others/FukushimaLambertW.cc")

# === Libraries ===
target_link_libraries(tests PUBLIC MuirLambertW)
target_include_directories(tests PUBLIC "../include/")

find_package(PkgConfig)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
target_link_libraries(tests PRIVATE PkgConfig::mpfr)

# === Feature Enables ===
target_compile_features(tests PUBLIC cxx_std_20)
enable_ipo(tests)
set_arch(tests)