# EXECUTABLE PROJECT - example

# === Includes ===
cmake_policy(SET CMP0069 NEW)
include(../macros.cmake)

# === Create Executable ===
add_executable(example "example.cpp")

# === Libraries ===
target_link_libraries(example PUBLIC MuirLambertW)
target_include_directories(example PUBLIC "../include/")

# === Feature Enables ===
target_compile_features(example PUBLIC cxx_std_20)
enable_ipo(example)
set_arch(example)