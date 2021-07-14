cmake_minimum_required(VERSION 3.4)
project(SPH)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules/" ${CMAKE_MODULE_PATH})

include_directories("${CMAKE_SOURCE_DIR}/include")

add_subdirectory(extern/glfw/ EXCLUDE_FROM_ALL)

# Add Corrade as a subproject
add_subdirectory(extern/corrade EXCLUDE_FROM_ALL)

# Add Magnum as a subproject, enable Sdl2Application
#set(WITH_SDL2APPLICATION ON CACHE BOOL "" FORCE)
set(WITH_CGLCONTEXT ON CACHE BOOL "" FORCE)
set(WITH_GLFWAPPLICATION ON CACHE BOOL "" FORCE)
add_subdirectory(extern/magnum EXCLUDE_FROM_ALL)

set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/imgui)
set(WITH_IMGUI ON CACHE BOOL "" FORCE)
add_subdirectory(extern/magnum-integration EXCLUDE_FROM_ALL)

find_package(Corrade REQUIRED Main)

find_package(Magnum REQUIRED
    DebugTools
    GL
    MeshTools
    Primitives
    SceneGraph
    Shaders
    GlfwApplication)

find_package(MagnumIntegration REQUIRED ImGui)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_subdirectory(lib)
add_subdirectory(src)
