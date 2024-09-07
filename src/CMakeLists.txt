# LIBRARY PROJECT - MuirLambertW

# === Includes ===
cmake_policy(SET CMP0069 NEW)
include(../macros.cmake)
include(GNUInstallDirs)

# === Create Executable ===
add_library(MuirLambertW  "MuirLambertWSimd.cpp" "MuirLambertWSimd.h" "mysleef.h")

# === Link Libraries ===
find_package(sleef REQUIRED)
target_link_libraries(MuirLambertW PRIVATE sleef::sleef)
target_include_directories(MuirLambertW PRIVATE ${SLEEF_INCLUDE_DIR})

# === Feature Enables ===
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