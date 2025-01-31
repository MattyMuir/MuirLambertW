# EXECUTABLE PROJECT - example

# === Create Executable ===
add_executable(example "example.cpp" "../bench/others/FukushimaLambertW.cc" "../bench/others/VebericLambertW.cc" "../bench/others/VebericLambertWOld.cc" "../bench/others/BarryLambertW.cpp" "../bench/others/MuirFukushima.cpp")

# === Libraries ===
target_link_libraries(example PUBLIC MuirLambertW)
target_include_directories(example PUBLIC "../include/")

find_package(ReferenceLambertW)
target_link_libraries(example PRIVATE ReferenceLambertW::ReferenceLambertW)

find_package(flttestlib)
target_link_libraries(example PRIVATE flttestlib::flttestlib)

find_package(boost_math REQUIRED CONFIG)
target_include_directories(example PRIVATE Boost::math)

# === Feature Enables ===
if (MUIRW_MSVC_STATIC_RUNTIME)
    set_property(TARGET example PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_compile_features(example PUBLIC cxx_std_20)
enable_ipo(example)
set_arch(example)