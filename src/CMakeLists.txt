# LIBRARY PROJECT - MuirLambertW

# === Includes ===
include(GNUInstallDirs)

# === Create Library ===
add_library(MuirLambertW  "w0.cpp"   "wm1.cpp"   "wm1serial.cpp" "w0serial.cpp"  "wm1fserial.cpp"  "wm1f.cpp" "w0f.cpp" "w0fserial.cpp" "w0fver2.cpp")

find_package(sleef REQUIRED)
target_link_libraries(MuirLambertW PUBLIC sleef::sleef)

# === Feature Enables ===
if (MUIRW_MSVC_STATIC_RUNTIME)
    set_property(TARGET MuirLambertW PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

set_property(TARGET MuirLambertW PROPERTY CXX_STANDARD 20)
enable_ipo(MuirLambertW)
set_arch(MuirLambertW)

# === Install ===
# paths for binaries and headers
install(TARGETS ${PROJECT_NAME}
    EXPORT "${PROJECT_NAME}Config"
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME} # include/SomeLibrary
)

# Export includes
install(FILES "../include/MuirLambertW.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY "../src"
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        FILES_MATCHING
        PATTERN "*.h")

# CMake config name, namespace and path
install(
    EXPORT "${PROJECT_NAME}Config"
    FILE "${PROJECT_NAME}Config.cmake"
    NAMESPACE MuirLambertW::
    DESTINATION cmake
)