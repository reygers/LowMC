#pragma once

#include "lowmc/FourVector.h"

namespace lowmc {

struct Photon {
    FourVector momentum;
    double weight = 1.0;
};

}  // namespace lowmc
