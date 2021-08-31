#pragma once

#include <Magnum/Math/Vector3.h>
#include <Magnum/Math/Color.h>

namespace Magnum {

using namespace Math::Literals;

struct ParticleVertex {
    Vector3 position;
    Color3 color;
};

}

