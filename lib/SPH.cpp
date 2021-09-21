#include "SPH.h"

#include <algorithm>

namespace Magnum {

float maxTemp = 0.0f;

Color3 SPH::floatToColor(float temp, float min, float max) {
    float val = max - min;
    float temp0 = temp - min;
    float deg = (temp0 / val) * 360.0f;

    return Color3::fromHsv(ColorHsv(Deg(deg), 0.7f, 0.5f));
}

Color3 SPH::positionToColor(Vector3 vec) {
    return Color3 { vec[0], vec[1], vec[2] };
}

SPH::SPH(Float particleRadius, size_t n) :
    _accelerations(n, {0.0f, 0.0f, 0.0f}),
    _velocities(n, {0.0f, 0.0f, 0.0f}),
    _densities(n, 0.0f),
    _densityCorrections(n, 0.0f),
    _temperatures(n, 0.0f),
    _dTemp(n, 0.0f),
    _f(n, {0.0f, 0.0f, 0.0f}),
    _particleRadius(particleRadius),
    _particleRadiusSqr(particleRadius * particleRadius),
    _kernels(_particleRadius),
    _maxParticles(n) {}

void SPH::addParticle(Vector3 position, Vector3 initialVelocity, float temperature)
{
    if (_numParticles < _maxParticles) {
        _particleVertex.push_back({position, 0x00ff00_rgbf});
        _virtualParticle.push_back({0.0f, _delta, 0.0f});
        _velocities[_numParticles] = initialVelocity;
        _temperatures[_numParticles] = temperature;
        _numParticles++;
    }
}

void SPH::computeUncorrectedDensities()
{
    for (size_t i = 0; i < _particleVertex.size(); i++) {
        _densities[i] = _kernels.W0();

        for (size_t j = 0; j < _particleVertex.size(); j++) {
            if (i == j) continue;

            Vector3 diff = _particleVertex[i].position - _particleVertex[j].position;

            Float r = diff.length();

            if (r > 0.0f && r <= _particleRadius) {
                _densities[i] += _particleMass * _kernels.W(diff);
            }
        }
    }
}

void SPH::computeDensities()
{
    for (size_t i = 0; i < _particleVertex.size(); i++) {
        Vector3 aux (0.0f, 0.0f, 0.0f);

        for (size_t j = 0; j < _particleVertex.size(); j++) {
            if (i == j) continue;

            Vector3 diff = _particleVertex[i].position - _particleVertex[j].position;

            Float r = diff.length();

            if (r > 0.0f && r <= _particleRadius) {
                aux -= _particleMass / _densities[j] * _kernels.poly6GradW(diff);
            }

        }

        Float v0 = aux.length() / _kernels.poly6GradW(_virtualParticle[i]).length();
        _densityCorrections[i] = _densities[i] * (1.0f + v0 * _kernels.W(_virtualParticle[i]));
    }
}

void SPH::computeForces()
{
    for (size_t i = 0; i < _particleVertex.size(); i++) {
        Vector3 damping {0.0f, 0.0f, 0.0f};
        Vector3 buoyancy {0.0f, 0.0f, 0.0f};

        Vector3 fAtmosfericPressure {0.0f, 0.0f, 0.0f};
        Vector3 fPressure {0.0f, 0.0f, 0.0f};
        Vector3 fViscosity {0.0f, 0.0f, 0.0f};

        Float dTemperature = 0.0f;

        for (size_t j = 0; j < _particleVertex.size(); j++) {
            if (i == j) continue;

            Float pi = computePressure(_densities[i]);
            Float pj = computePressure(_densities[j]);
            Float pk = computePressure(_densityCorrections[i]);

            Vector3 diff = _particleVertex[i].position - _particleVertex[j].position;

            Float r = diff.dot();

            if (r > 0.0f && r <= _particleRadiusSqr) {
                fAtmosfericPressure += (_particleMass / pj) * _kernels.spikyGradW(diff);

                fPressure -= (_particleMass / _densities[j]) * ((pi + pj) / 2.0f) * _kernels.spikyGradW(diff)
                    + (_particleMass / _densityCorrections[i]) * ((pi + pk) / 2.0f) * _kernels.spikyGradW(_virtualParticle[i]);

                fViscosity += _particleMass * (_velocities[j] - _velocities[i]) / _densities[j] * _kernels.laplW(diff);

                dTemperature += (_particleMass / (pi * pj)) * _thermalConductivity * (_temperatures[i] - _temperatures[j])
                    * (dot(diff, _kernels.spikyGradW(diff)) / (diff.dot() + _smallPositive));
            }
        }

        fViscosity *= _particleViscosity;

        if (fAtmosfericPressure.length() > _dampingThreshold) {
          dTemperature -= _temperatures[i] / _radiationHalfLife;
          damping = -_dampingCoefficient * _velocities[i];
        }

        _dTemp[i] = dTemperature;

        buoyancy = _buoyancyCoefficient * _temperatures[i] * _buoyancyDirection;

        _f[i] = (fPressure + 1.0f * fAtmosfericPressure) + fViscosity + _densities[i] * (_gravity +  buoyancy + damping);
    }
}

Float SPH::computePressure(Float density)
{
    return _gasConstant * (density - _restDensity);
}

void SPH::integrate(Float timestep)
{
    timestep = 0.001;

    for (size_t i = 0; i < _particleVertex.size(); i++) {
        Vector3 prevAccel = _accelerations[i];
        Vector3 prevVel = _velocities[i];
        //float prevTemp = _temperatures[i];
        //float newTemp = _temperatures[i] + _dTemp[i];

        //_accelerations[i] = _f[i] / _densities[i];
        //Vector3 position = _particleVertex[i].position + _velocities[i] * timestep + _accelerations[i] * timestep * timestep;
        //Vector3 velocity = (position - _particleVertex[i].position) / timestep;
        //_velocities[i] = velocity;
        //_particleVertex[i].position = position;

        _temperatures[i] += _dTemp[i];// * timestep;

        _minTemp = std::min(_minTemp, _temperatures[i]);
        _maxTemp = std::max(_maxTemp, _temperatures[i]);

        _minDensity = std::min(_minDensity, _temperatures[i]);
        _maxDensity = std::max(_maxDensity, _temperatures[i]);

        //std::cout << "temp[" << i << "] "; std::cout << _temperatures[i] << std::endl;
        //_temperatures[i] += (prevTemp + newTemp) / 2.0 * timestep;

        _accelerations[i] = _f[i] / _densities[i];
        _velocities[i] += (prevAccel + _accelerations[i]) / 2.0 * timestep;
        _particleVertex[i].position += prevVel * timestep + prevAccel / 2.0 * timestep * timestep;

        //_accelerations[i] = _f[i] / _densities[i];
        //_velocities[i] += _accelerations[i] * timestep;
        //_particleVertex[i].position += _velocities[i] * timestep;

        checkCollisions(i);
    }
}

void SPH::checkCollisions(size_t i)
{
    Vector3 position = _particleVertex[i].position;
    Vector3 velocity = _velocities[i];
    Float r = _particleRadius;

    Float force = -0.1f;

    if (position[0] > 1.0f - r)
    {
        position[0] = 1.0f - r;
        velocity[0] *= force;
    }
    if (position[0] < r)
    {
        position[0] = r;
        velocity[0] *= force;
    }

    if (position[1] > 1.0f - r)
    {
        position[1] = 1.0f - r;
        velocity[1] *= force;
    }

    if (position[1] < r)
    {
        position[1] = r;
        velocity[1] *= force;
    }

    if (position[2] > 1.0f - r)
    {
        position[2] = 1.0f - r;
        velocity[2] *= force;
    }

    if (position[2] < r)
    {
        position[2] = r;
        velocity[2] *= force;
    }

    _particleVertex[i].position = position;

    /*
    if (_temperatures[i] == 25.0f) {
        _particleVertex[i].color = 0x0000ff_rgbf;
    }
    else {
        _particleVertex[i].color = 0x00ff00_rgbf;
    }
    */

    float velLength = velocity.length();
    float accLength = _accelerations[i].length();

    _minVel = std::min(velLength, _minVel);
    _maxVel = std::max(velLength, _maxVel);

    _minAcc = std::min(velLength, _minAcc);
    _maxAcc = std::max(velLength, _maxAcc);

    if (_showDensity) _particleVertex[i].color = floatToColor(_densities[i], _minDensity, _maxDensity);
    if (_showTemperature) _particleVertex[i].color = floatToColor(_temperatures[i], _minTemp, _maxTemp);
    if (_showAcceleration) _particleVertex[i].color = floatToColor(accLength, _minAcc, _maxAcc);
    if (_showVelocity) _particleVertex[i].color = floatToColor(velLength, _minVel, _maxVel);
    if (_showPosition) _particleVertex[i].color = positionToColor(_particleVertex[i].position);
    _velocities[i] = velocity;
}

void SPH::step(Float timestep)
{
    computeUncorrectedDensities();
    computeDensities();
    computeForces();
    integrate(timestep);
}

void SPH::showTemperature()
{
    _showTemperature = true;
    _showDensity = false;
    _showVelocity = false;
    _showAcceleration = false;
    _showPosition = false;
}

void SPH::showDensity()
{
    _showTemperature = false;
    _showDensity = true;
    _showVelocity = false;
    _showAcceleration = false;
    _showPosition = false;
}

void SPH::showVelocity()
{
    _showTemperature = false;
    _showDensity = false;
    _showVelocity = true;
    _showAcceleration = false;
    _showPosition = false;
}

void SPH::showAcceleration()
{
    _showTemperature = false;
    _showDensity = false;
    _showVelocity = false;
    _showAcceleration = true;
    _showPosition = false;
}

void SPH::showPosition()
{
    _showTemperature = false;
    _showDensity = false;
    _showVelocity = false;
    _showAcceleration = false;
    _showPosition = true;
}

}

