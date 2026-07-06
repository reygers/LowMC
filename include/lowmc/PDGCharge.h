#pragma once

#include <cstdlib>

namespace lowmc {

struct PDGChargeLookupResult {
    bool known = false;
    int chargeTimesThree = 0;
};

[[nodiscard]] constexpr int quarkChargeTimesThree(int quarkDigit) {
    switch (quarkDigit) {
    case 1:
    case 3:
    case 5:
        return -1;
    case 2:
    case 4:
    case 6:
        return +2;
    default:
        return 0;
    }
}

[[nodiscard]] constexpr bool isSupportedQuarkDigit(int quarkDigit) {
    return quarkDigit >= 1 && quarkDigit <= 6;
}

[[nodiscard]] inline PDGChargeLookupResult chargeTimesThreeFromPDG(int pdgId) {
    if (pdgId == 0) {
        return {};
    }

    const int sign = pdgId > 0 ? 1 : -1;
    const int absId = std::abs(pdgId);

    if (absId >= 1'000'000'000 && absId < 1'100'000'000) {
        const int protonNumber = (absId / 10'000) % 1000;
        return {true, sign * 3 * protonNumber};
    }

    switch (absId) {
    case 1:
    case 3:
    case 5:
        return {true, sign * -1};
    case 2:
    case 4:
    case 6:
        return {true, sign * +2};
    case 11:
    case 13:
    case 15:
        return {true, sign * -3};
    case 12:
    case 14:
    case 16:
    case 21:
    case 22:
    case 23:
    case 25:
        return {true, 0};
    case 24:
        return {true, sign * +3};
    default:
        break;
    }

    if (absId >= 10'000'000) {
        return {};
    }

    const int spinDigit = absId % 10;
    const int q3 = (absId / 10) % 10;
    const int q2 = (absId / 100) % 10;
    const int q1 = (absId / 1000) % 10;

    if (spinDigit == 0) {
        return {};
    }

    if (q1 == 0 && isSupportedQuarkDigit(q2) && isSupportedQuarkDigit(q3)) {
        int charge = 0;
        if (q2 % 2 == 1 && q3 % 2 == 0) {
            charge = quarkChargeTimesThree(q3) - quarkChargeTimesThree(q2);
        } else {
            charge = quarkChargeTimesThree(q2) - quarkChargeTimesThree(q3);
        }
        return {true, sign * charge};
    }

    if (isSupportedQuarkDigit(q1) &&
        isSupportedQuarkDigit(q2) &&
        isSupportedQuarkDigit(q3)) {
        return {
            true,
            sign * (
                quarkChargeTimesThree(q1) +
                quarkChargeTimesThree(q2) +
                quarkChargeTimesThree(q3)),
        };
    }

    return {};
}

}  // namespace lowmc
