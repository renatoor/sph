#pragma once

#include <iostream>

#include <vector>
#include <Magnum/Magnum.h>

#include "SPHKernels.h"

namespace Magnum {

class SPH {
public:

    //SPH(Float particleRadius, std::vector<Vector3> &positions);
    SPH(Float particleRadius, size_t n);

    const std::vector<Vector3>& particlePositions() { return _positions; }

    void step(Float timestep);

    void addParticle(Vector3 position, Vector3 initialVelocity);

private:

    void computeUncorrectedDensities();

    void computeDensities();

    void computeForces();

    Float computePressure(Float density);

    void integrate(Float timestep);

    void checkCollisions(size_t idx);

    std::vector<Vector3> _positions; // _x
    std::vector<Vector3> _accelerations; // _a
    std::vector<Vector3> _velocities; // _v
    std::vector<float> _densities; // _rho
    std::vector<float> _densityCorrections; // _rho
    std::vector<Vector3> _f;

    Vector3 _gravity = Vector3(0.0f, -9.82f, 0.0f);

    Float _particleRadius; // _h
    Float _particleRadiusSqr; // _h2
    Float _deltaParticleRadius;
    Float _virtualDistance;
    Float _virtualDistanceSqr;
    Float _particleMass = 0.02f; // _m
    //Float _particleViscosity = 3.5f; // _mu
    Float _particleViscosity = 0.01f; // _mu
    Float _gasConstant = 3.0f; // _k
    Float _restDensity = 998.29f / 1000.0f; // _rho0
    //Float _surfaceTension = 0.0728f;

    SPHKernels _kernels;

    size_t _maxParticles;
    size_t _numParticles = 0;
    Float _delta = 0.4f;
    Vector3 virtualParticle = Vector3(0.01f, _delta, 0.01f);

};

}

