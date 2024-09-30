# EXECUTABLE PROJECT - bench

# === Includes ===
cmake_policy(SET CMP0069 NEW)
include(../macros.cmake)

# === Create Executable ===
add_executable(benchw0 "benchw0.cpp" "Timer.h" "others/BarryLambertW.cpp" "others/BarryLambertW.h" "others/FukushimaLambertW.cc" "others/FukushimaLambertW.h" "../tests/ReciprocalDistributionEx.cpp" "../tests/ReciprocalDistribution.cpp")
add_executable(benchwm1 "benchwm1.cpp" "Timer.h" "others/BarryLambertW.cpp" "others/BarryLambertW.h" "others/FukushimaLambertW.cc" "others/FukushimaLambertW.h" "../tests/ReciprocalDistributionEx.cpp" "../tests/ReciprocalDistribution.cpp")
add_executable(arraybench "arraybench.cpp" "others/BarryLambertW.cpp" "others/BarryLambertW.h" "others/FukushimaLambertW.cc" "others/FukushimaLambertW.h")

# === Libraries ===
target_link_libraries(benchw0 PUBLIC MuirLambertW)
target_include_directories(benchw0 PUBLIC "../include/")
target_include_directories(benchw0 PRIVATE "C:/sdks/boost_1_86_0")

target_link_libraries(benchwm1 PUBLIC MuirLambertW)
target_include_directories(benchwm1 PUBLIC "../include/")
target_include_directories(benchwm1 PRIVATE "C:/sdks/boost_1_86_0")

target_link_libraries(arraybench PUBLIC MuirLambertW)
target_include_directories(arraybench PUBLIC "../include/")
target_include_directories(arraybench PRIVATE "C:/sdks/boost_1_86_0")

# === Feature Enables ===
target_compile_features(benchw0 PUBLIC cxx_std_20)
enable_ipo(benchw0)
set_arch(benchw0)

target_compile_features(benchwm1 PUBLIC cxx_std_20)
enable_ipo(benchwm1)
set_arch(benchwm1)

target_compile_features(arraybench PUBLIC cxx_std_20)
enable_ipo(arraybench)
set_arch(arraybench)