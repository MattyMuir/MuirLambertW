# EXECUTABLE PROJECT - bench

# === Includes ===
cmake_policy(SET CMP0069 NEW)
include(../macros.cmake)

# === Create Executable ===
add_executable(benchw0 "benchw0.cpp" "Timer.h" "others/BarryLambertW.cpp" "others/BarryLambertW.h" "others/FukushimaLambertW.cc" "others/FukushimaLambertW.h")
add_executable(benchwm1 "benchwm1.cpp" "Timer.h")

# === Libraries ===
target_link_libraries(benchw0 PUBLIC MuirLambertW)
target_include_directories(benchw0 PUBLIC "../include/")

target_link_libraries(benchwm1 PUBLIC MuirLambertW)
target_include_directories(benchwm1 PUBLIC "../include/")

# === Feature Enables ===
target_compile_features(benchw0 PUBLIC cxx_std_20)
enable_ipo(benchw0)
set_arch(benchw0)

target_compile_features(benchwm1 PUBLIC cxx_std_20)
enable_ipo(benchwm1)
set_arch(benchwm1)