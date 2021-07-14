corrade_add_resource(FluidSimulation_RESOURCES ../resources.conf)

add_executable(sph main.cpp ${FluidSimulation_RESOURCES})

target_link_libraries(sph PRIVATE
    sphlib
    Magnum::Application
    Magnum::Primitives
    Magnum::SceneGraph
    MagnumIntegration::ImGui)

# Make the executable a default target to build & run in Visual Studio
set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT sph)
