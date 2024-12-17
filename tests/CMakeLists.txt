# EXECUTABLE PROJECT - tests

# === Create Executable ===
add_executable(tests "tests.cpp"         "../bench/others/BarryLambertW.cpp" "../bench/others/FukushimaLambertW.cc"   "constants.h")

# === Libraries ===
target_link_libraries(tests PUBLIC MuirLambertW)
target_include_directories(tests PUBLIC "../include/")

find_package(ReferenceLambertW)
target_link_libraries(tests PRIVATE ReferenceLambertW::ReferenceLambertW)

find_package(flttestlib)
target_link_libraries(tests PRIVATE flttestlib::flttestlib)

find_package(boost_math REQUIRED CONFIG)
target_include_directories(tests PRIVATE Boost::math)

# === Feature Enables ===
if (MUIRW_MSVC_STATIC_RUNTIME)
    set_property(TARGET tests PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_compile_features(tests PUBLIC cxx_std_20)
enable_ipo(tests)
set_arch(tests)