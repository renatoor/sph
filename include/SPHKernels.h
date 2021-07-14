#pragma once

#include <Magnum/Magnum.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/Math/Vector3.h>

namespace Magnum {

class Poly6Kernel {
public:

    Poly6Kernel(Float radius) {
        _radius = radius;
        _radiusSqr = _radius * _radius;
        _k = 315.0f / (64.0f * Constants::pi() * Math::pow(_radius, 9.0f));
        _l = -945.0f / (32.0f * Constants::pi() * Math::pow(_radius, 9.0f));
        _W0 = W(0.0f);
    }

    Float W(const Float r) const {
        const Float r2 = r * r;
        return _k * Math::pow(_radiusSqr - r2, 3.0f);
    }

    Float W(const Vector3& r) const {
        const auto r2 = r.dot();
        return _k * Math::pow(_radiusSqr - r2, 3.0f);
    }

    Vector3 gradW(const Vector3& r) const {
        const auto r2 = r.dot();
        const auto hr = _radiusSqr - r2;
        const auto hr2 = hr * hr;
        return r * _l * hr2;
    }

    Float laplW(const Vector3& r) const {
        const auto r2 = r.dot();
        return _l * (_radiusSqr - r2) * (3.0 * _radiusSqr - 7.0 * r2);
    }

    Float W0() const { return _W0; }

private:

    Float _radius;
    Float _radiusSqr;
    Float _k;
    Float _l;
    Float _W0;

};

class SpikyKernel {
public:

    SpikyKernel(Float radius) {
        _radius = radius;
        _radiusSqr = _radius * _radius;
        _l = -45.0f / (Constants::pi() * Math::pow(_radius, 6.0f));
    }

    Vector3 gradW(const Vector3& r) const {
        const Float rl = r.length();
        const Float hr = _radius - rl;
        const Float hr2 = hr * hr;
        return _l * hr2 * (r / rl);
    }

private:

    Float _radius;
    Float _radiusSqr;
    Float _l;

};

class ViscosityKernel {
public:
    ViscosityKernel(Float radius) :
        _radius(radius),
        _l(45.0f / (Constants::pi() * Math::pow(_radius, 6.0f))) {}

    Float laplW(const Vector3& r) const {
        const Float rl = r.length();
        return _l * (_radius * rl);
    }

private:

    Float _radius;
    Float _l;

};

class SPHKernels {
public:

    SPHKernels(Float radius) : _poly6(radius), _spiky(radius), _viscosity(radius) {}

    Float W0() const { return _poly6.W0(); }

    Float W(const Vector3& r) const { return _poly6.W(r); }

    Vector3 poly6GradW(const Vector3& r) const { return _poly6.gradW(r); }

    Float poly6LaplW(const Vector3& r) const { return _poly6.laplW(r); }

    Vector3 spikyGradW(const Vector3& r) const { return _spiky.gradW(r); }

    Float laplW(const Vector3& r) const { return _viscosity.laplW(r); }

private:

    Poly6Kernel _poly6;
    SpikyKernel _spiky;
    ViscosityKernel _viscosity;

};

}

