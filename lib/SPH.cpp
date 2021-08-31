#include "SPH.h"

namespace Magnum {

SPH::SPH(Float particleRadius, size_t n) :
    _accelerations(n, {0.0f, 0.0f, 0.0f}),
    _velocities(n, {0.0f, 0.0f, 0.0f}),
    _densities(n, 0.0f),
    _densityCorrections(n, 0.0f),
    _temperatures(n, 0.0f),
    _f(n, {0.0f, 0.0f, 0.0f}),
    _particleRadius(particleRadius),
    _particleRadiusSqr(particleRadius * particleRadius),
    _deltaParticleRadius(_delta * particleRadius),
    _virtualDistance(_particleRadius * 0.4),
    _virtualDistanceSqr(_virtualDistance * _virtualDistance),
    _kernels(_particleRadius),
    _maxParticles(n) {}

void SPH::addParticle(Vector3 position, Vector3 initialVelocity)
{
    if (_numParticles < _maxParticles) {
        _particleVertex.push_back({position, 0x00ff00_rgbf});
        _velocities[_numParticles] = initialVelocity;
        _temperatures[_numParticles] = 25.0f;
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

        Float v0 = aux.length() / _kernels.poly6GradW(virtualParticle).length();
        _densityCorrections[i] = _densities[i] * (1.0f + v0 * _kernels.W(virtualParticle));
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
                    + (_particleMass / _densityCorrections[i]) * ((pi + pk) / 2.0f) * _kernels.spikyGradW(virtualParticle);

                fViscosity += _particleMass * (_velocities[j] - _velocities[i]) / _densities[j] * _kernels.laplW(diff);

                dTemperature += (_particleMass / (pi * pj)) * _thermalConductivity * (_temperatures[i] - _temperatures[j])
                    * ((dot(diff, _kernels.spikyGradW(diff))) / (diff.dot() * _smallPositive));
            }
        }

        fViscosity *= _particleViscosity;

        if (fAtmosfericPressure.length() > _dampingThreshold) {
            //_temperatures[i] = dTemperature * _radiationHalfLife * -1.0f;
            dTemperature -= _temperatures[i] / _radiationHalfLife;
            damping = -_dampingCoefficient * _velocities[i];
        }

        _temperatures[i] += dTemperature;
        //std::cout << "d temp " << dTemperature << std::endl;

        //if (isnan(dTemperature)) exit(0);
        //std::cout << "temp " << _temperatures[i] << std::endl;

        //buoyancy = _buoyancyCoefficient * _temperatures[i] * _buoyancyDirection;
        buoyancy = _buoyancyCoefficient * 50.0f * _buoyancyDirection;

        //std::cout << buoyancy[0] << " " << buoyancy[1] << " " << buoyancy[2] << std::endl;

        _f[i] = (fPressure + 1.0f * fAtmosfericPressure) + fViscosity + _densities[i] * (_gravity);// +  buoyancy);

        //std::cout << _f[0] << " " << _f[1] << " " << _f[2] << std::endl;
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

        //_accelerations[i] = _f[i] / _densities[i];
        //Vector3 position = _particleVertex[i].position + _velocities[i] * timestep + _accelerations[i] * timestep * timestep;
        //Vector3 velocity = (position - _particleVertex[i].position) / timestep;
        //_velocities[i] = velocity;
        //_particleVertex[i].position = position;

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
    if (i % 2 == 0) {
        _particleVertex[i].color = 0x0000ff_rgbf;
    }
    else {
        _particleVertex[i].color = 0x00ff00_rgbf;
    }
    _velocities[i] = velocity;
}

void SPH::step(Float timestep)
{
    computeUncorrectedDensities();
    computeDensities();
    computeForces();
    integrate(timestep);
}

}

