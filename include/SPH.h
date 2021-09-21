#pragma once

#include <iostream>

#include <vector>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>

#include "SPHKernels.h"

#include "Objects/ParticleVertex.h"

namespace Magnum {

class SPH {
public:

    //SPH(Float particleRadius, std::vector<Vector3> &positions);
    SPH(Float particleRadius, size_t n);

    //const std::vector<Vector3>& particlePositions() { return _positions; }

    const std::vector<ParticleVertex>& particleVertex() { return _particleVertex; }

    void step(Float timestep);

    void addParticle(Vector3 position, Vector3 initialVelocity, float temperature = 25.0f);

    void showTemperature();
    void showDensity();
    void showVelocity();
    void showAcceleration();
    void showPosition();

private:

    void computeUncorrectedDensities();

    void computeDensities();

    void computeForces();

    Float computePressure(Float density);

    void integrate(Float timestep);

    void checkCollisions(size_t idx);

    Color3 floatToColor(float temp, float min, float max);

    Color3 positionToColor(Vector3 vec);

    //std::vector<Vector3> _positions; // _x
    std::vector<ParticleVertex> _particleVertex; // _x
    std::vector<Vector3> _virtualParticle;
    std::vector<Vector3> _accelerations; // _a
    std::vector<Vector3> _velocities; // _v
    std::vector<float> _densities; // _rho
    std::vector<float> _densityCorrections; // _rho
    std::vector<float> _temperatures; // _rho
    std::vector<float> _dTemp; // _rho
    std::vector<Vector3> _f;

    Vector3 _gravity = Vector3(0.0f, -9.82f, 0.0f);

    Float _particleRadius; // _h
    Float _particleRadiusSqr; // _h2
    //Float _deltaParticleRadius;
    //Float _virtualDistance;
    //Float _virtualDistanceSqr;
    Float _particleMass = 0.02f; // _m
    //Float _particleViscosity = 3.5f; // _mu
    Float _particleViscosity = 0.01f; // _mu
    Float _gasConstant = 3.0f; // _k
    Float _restDensity = 998.29f / 1000.0f; // _rho0
    //Float _surfaceTension = 0.0728f;

    Float _delta = 0.4f * _particleRadius;

    SPHKernels _kernels;

    size_t _maxParticles;
    size_t _numParticles = 0;
    //Vector3 virtualParticle = Vector3(0.0f, _delta, 0.0f);

    Float _dampingCoefficient = 100.0f;
    Float _dampingThreshold = 1.0f;

    //Float _buoyancyCoefficient = 0.392f; // aceleração bouyancy compensar a gravidade e particula ficar parada
    //Float _buoyancyCoefficient = 1.0f;
    Float _buoyancyCoefficient = 0.8f;
    //Float _temperature = 25.0f;
    //Float _thermalConductivity = 0.003;//0.598f;
    Float _thermalConductivity = 0.598f;
    Float _radiationHalfLife = 1000.0f;
    Float _smallPositive = 0.0001;

    Vector3 _buoyancyDirection = Vector3(0.0f, 1.0f, 0.0f);

    float _minTemp = 0.0f;
    float _maxTemp = 0.0f;

    float _minDensity = 0.0f;
    float _maxDensity = 0.0f;

    float _minVel = 0.0f;
    float _maxVel = 0.0f;

    float _minAcc = 0.0f;
    float _maxAcc = 0.0f;

    bool _showTemperature = false;
    bool _showDensity = false;
    bool _showVelocity = false;
    bool _showAcceleration = false;
    bool _showPosition = false;

};

}

