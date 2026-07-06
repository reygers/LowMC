#pragma once

#include "lowmc/FourVector.h"

namespace lowmc {

enum class InOutSign : int {
    Incoming = -1,
    Outgoing = 1,
};

[[nodiscard]] constexpr int signValue(InOutSign sign) {
    return static_cast<int>(sign);
}

struct ChargedParticle {
    FourVector momentum;
    double charge = 0.0;
    InOutSign inOutSign = InOutSign::Outgoing;
};

}  // namespace lowmc
