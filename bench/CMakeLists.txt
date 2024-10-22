# EXECUTABLE PROJECT - bench

# === Includes ===
cmake_policy(SET CMP0069 NEW)
include(../macros.cmake)

# === Create Executable ===
add_executable(benchw0 "benchw0.cpp" "Timer.h" "others/BarryLambertW.cpp" "others/BarryLambertW.h" "others/FukushimaLambertW.cc" "others/FukushimaLambertW.h" "../tests/ReciprocalDistributionEx.cpp" "../tests/ReciprocalDistribution.cpp")
add_executable(benchwm1 "benchwm1.cpp" "Timer.h" "others/BarryLambertW.cpp" "others/BarryLambertW.h" "others/FukushimaLambertW.cc" "others/FukushimaLambertW.h" "../tests/ReciprocalDistributionEx.cpp" "../tests/ReciprocalDistribution.cpp")
add_executable(arraybench "arraybench.cpp" "others/BarryLambertW.cpp" "others/BarryLambertW.h" "others/FukushimaLambertW.cc" "others/FukushimaLambertW.h")

# === Libraries ===
find_package(boost_math REQUIRED CONFIG)

target_link_libraries(benchw0 PUBLIC MuirLambertW)
target_include_directories(benchw0 PUBLIC "../include/")
target_link_libraries(benchw0 PRIVATE Boost::math)

target_link_libraries(benchwm1 PUBLIC MuirLambertW)
target_include_directories(benchwm1 PUBLIC "../include/")
target_link_libraries(benchwm1 PRIVATE Boost::math)

target_link_libraries(arraybench PUBLIC MuirLambertW)
target_include_directories(arraybench PUBLIC "../include/")
target_link_libraries(arraybench PRIVATE Boost::math)

# === Feature Enables ===
if (MUIRW_MSVC_STATIC_RUNTIME)
    set_property(TARGET benchw0 PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_compile_features(benchw0 PUBLIC cxx_std_20)
enable_ipo(benchw0)
set_arch(benchw0)

if (MUIRW_MSVC_STATIC_RUNTIME)
    set_property(TARGET benchwm1 PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_compile_features(benchwm1 PUBLIC cxx_std_20)
enable_ipo(benchwm1)
set_arch(benchwm1)

if (MUIRW_MSVC_STATIC_RUNTIME)
    set_property(TARGET arraybench PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_compile_features(arraybench PUBLIC cxx_std_20)
enable_ipo(arraybench)
set_arch(arraybench)