# EXECUTABLE PROJECT - example

# === Create Executable ===
add_executable(example "example.cpp")

# === Libraries ===
target_link_libraries(example PUBLIC MuirLambertW)
target_include_directories(example PUBLIC "../include/")

# === Feature Enables ===
if (MUIRW_MSVC_STATIC_RUNTIME)
    set_property(TARGET example PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_compile_features(example PUBLIC cxx_std_20)
enable_ipo(example)
set_arch(example)