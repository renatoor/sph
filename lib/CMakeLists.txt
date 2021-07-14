add_library(sphlib
    Objects/ArcBall.cpp
    Objects/ParticleGroup.cpp
    Shaders/ParticleSphereShader.cpp
    SPH.cpp)

target_link_libraries(sphlib PRIVATE Magnum::Magnum)
