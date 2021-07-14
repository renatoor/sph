#include "SPH.h"

namespace Magnum {

SPH::SPH(Float particleRadius, size_t n) :
    _accelerations(n, {0.0f, 0.0f, 0.0f}),
    _velocities(n, {0.0f, 0.0f, 0.0f}),
    _densities(n, 0.0f),
    _densityCorrections(n, 0.0f),
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
        _positions.push_back(position);
        _velocities[_numParticles] = initialVelocity;
        _numParticles++;
    }
}

void SPH::computeUncorrectedDensities()
{
    for (size_t i = 0; i < _positions.size(); i++) {
        _densities[i] = _kernels.W0();

        for (size_t j = 0; j < _positions.size(); j++) {
            if (i == j) continue;

            Vector3 diff = _positions[i] - _positions[j];

            Float r = diff.length();

            if (r > 0.0f && r <= _particleRadius) {
                _densities[i] += _particleMass * _kernels.W(diff);
            }
        }
    }
}

void SPH::computeDensities()
{
    for (size_t i = 0; i < _positions.size(); i++) {
        Vector3 aux (0.0f, 0.0f, 0.0f);

        for (size_t j = 0; j < _positions.size(); j++) {
            if (i == j) continue;

            Vector3 diff = _positions[i] - _positions[j];

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
    for (size_t i = 0; i < _positions.size(); i++) {
        Vector3 fAtmosfericPressure {0.0f, 0.0f, 0.0f};
        Vector3 fPressure {0.0f, 0.0f, 0.0f};
        Vector3 fViscosity {0.0f, 0.0f, 0.0f};
        //Vector3 fSurfaceTension {0.0f, 0.0f, 0.0f};

        //int count = 0;
        //Vector3 surfaceTensionGrad {0.0f, 0.0f, 0.0f};
        //Float surfaceTensionLapl = 0.0f;

        for (size_t j = 0; j < _positions.size(); j++) {
            if (i == j) continue;

            Float pi = computePressure(_densities[i]);
            Float pj = computePressure(_densities[j]);
            Float pk = computePressure(_densityCorrections[i]);
            //Float pj = computePressure(_densityCorrections[j]);

            Vector3 diff = _positions[i] - _positions[j];

            Float r = diff.dot();
            if (r > 0.0f)
                fAtmosfericPressure += _particleMass / pj * _kernels.spikyGradW(diff);

            if (r > 0.0f && r <= _particleRadiusSqr) {
                fPressure -= _particleMass * (pi + pj) / (2.0f * _densities[j]) * _kernels.spikyGradW(diff)
                    + _particleMass * (pi + pk) / (2.0f * _densityCorrections[i]) * _kernels.spikyGradW(virtualParticle);
                fViscosity += _particleMass * (_velocities[j] - _velocities[i]) / _densities[j] * _kernels.laplW(diff);
                //surfaceTensionGrad += _kernels.poly6GradW(diff) / _densities[j];
                //surfaceTensionLapl += _kernels.poly6LaplW(diff) / _densities[j];
                //count++;
            }
        }

        //if (surfaceTensionGrad.length() >= Math::sqrt(_restDensity / count)) {
        //    fSurfaceTension = -surfaceTensionGrad / surfaceTensionGrad.length() * surfaceTensionLapl * _surfaceTension;
        //}

        fViscosity *= _particleViscosity;
        //_f[i] = fPressure + fViscosity + _densities[i] * _gravity;// + fSurfaceTension;
        _f[i] = (fPressure + 100.0f * fAtmosfericPressure) / _densities[i] + fViscosity / _densities[i] * _gravity;// + fSurfaceTension;
    }
}

Float SPH::computePressure(Float density)
{
    return _gasConstant * (density - _restDensity);
}

void SPH::integrate(Float timestep)
{
    timestep = 0.01;

    for (size_t i = 0; i < _positions.size(); i++) {
        Vector3 prevAccel = _accelerations[i];
        Vector3 prevVel = _velocities[i];

        //_accelerations[i] = _f[i] / _densities[i];
        //Vector3 position = _positions[i] + _velocities[i] * timestep + _accelerations[i] * timestep * timestep;
        //Vector3 velocity = (position - _positions[i]) / timestep;
        //_velocities[i] = velocity;
        //_positions[i] = position;

        _accelerations[i] = _f[i] / _densities[i];
        _velocities[i] += (prevAccel + _accelerations[i]) / 2.0 * timestep;
        _positions[i] += prevVel * timestep + prevAccel / 2.0 * timestep * timestep;

        //_accelerations[i] = _f[i] / _densities[i];
        //_velocities[i] += _accelerations[i] * timestep;
        //_positions[i] += _velocities[i] * timestep;

        checkCollisions(i);
    }
}

void SPH::checkCollisions(size_t i)
{
    Vector3 position = _positions[i];
    Vector3 velocity = _velocities[i];
    Float r = _particleRadius;

    Float force = -0.3f;

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

    _positions[i] = position;
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

