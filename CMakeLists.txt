# GLOBAL - MuirLambertW

# === Includes ===
include(macros.cmake)

# === Configuration ===
cmake_minimum_required(VERSION 3.13)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Enable multi-processor build
enable_mp()

# Options
option(WARNINGS "Enables warnings" ON)
if (WARNINGS)
	add_compile_options("-Wall")
    add_compile_options("-Wextra")
    add_compile_options("-Wpedantic")
    enable_clangtidy()
endif()

# === Projects ===
project("MuirLambertW")
add_subdirectory("src")

add_subdirectory("examples")
add_subdirectory("bench")

# === Testing ===
add_subdirectory("tests")