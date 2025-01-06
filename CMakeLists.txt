# GLOBAL - MuirLambertW

# === Includes ===
cmake_policy(SET CMP0069 NEW)
include(cmake/macros.cmake)

# === Configuration ===
cmake_minimum_required(VERSION 3.21)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Enable multi-processor build
enable_mp()

# Options
set(WARNINGS OFF CACHE BOOL "Enables clang-tidy and warnings on clang-like compilers")
if (WARNINGS)
	add_compile_options("-Wall")
    add_compile_options("-Wextra")
    add_compile_options("-Wpedantic")
    enable_clangtidy()
endif()

set(MUIRW_NO_ITERATOR_DEBUGGING OFF CACHE BOOL "Disable iterator debugging")
if (MUIRW_NO_ITERATOR_DEBUGGING)
    add_compile_definitions(_HAS_ITERATOR_DEBUGGING=0)
endif()

set(MUIRW_MSVC_STATIC_RUNTIME, OFF CACHE BOOL "Use static C runtime library when taregting MSVC ABI")
if (MUIRW_MSVC_STATIC_RUNTIME)
    use_static_msvc_crt()
endif()

# === Projects ===
project("MuirLambertW")
add_subdirectory("src")

add_subdirectory("examples")
add_subdirectory("bench")

# === Tests ===
enable_testing()
set(CTEST_OUTPUT_ON_FAILURE ON)

add_subdirectory("tests")

# Exhaustive tests
add_test(NAME W0fserialExhaustive COMMAND tests 0)
add_test(NAME Wm1fserialExhaustive COMMAND tests 1)
add_test(NAME W0fsimdExhaustive COMMAND tests 2)
add_test(NAME Wm1fsimdExhaustive COMMAND tests 3)

# Unit tests
add_test(NAME W0fserialUnit COMMAND tests 4)
add_test(NAME Wm1fserialUnit COMMAND tests 5)
add_test(NAME W0fsimdUnit COMMAND tests 6)
add_test(NAME Wm1fsimdUnit COMMAND tests 7)
add_test(NAME W0serialUnit COMMAND tests 8)
add_test(NAME Wm1serialUnit COMMAND tests 9)
add_test(NAME W0simdUnit COMMAND tests 10)
add_test(NAME Wm1simdUnit COMMAND tests 11)

# Sample tests
add_test(NAME W0serialSample COMMAND tests 12)
add_test(NAME Wm1serialSample COMMAND tests 13)
add_test(NAME W0simdSample COMMAND tests 14)
add_test(NAME Wm1simdSample COMMAND tests 15)