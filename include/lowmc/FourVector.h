#pragma once

#include <cmath>

namespace lowmc {

struct FourVector {
    double E = 0.0;
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;

    constexpr FourVector() = default;
    constexpr FourVector(double energy, double x, double y, double z)
        : E(energy), px(x), py(y), pz(z) {}

    [[nodiscard]] constexpr double dot(const FourVector& other) const {
        return E * other.E - px * other.px - py * other.py - pz * other.pz;
    }

    [[nodiscard]] constexpr double mass2() const {
        return dot(*this);
    }

    [[nodiscard]] double pt() const {
        return std::hypot(px, py);
    }

    constexpr FourVector& operator+=(const FourVector& other) {
        E += other.E;
        px += other.px;
        py += other.py;
        pz += other.pz;
        return *this;
    }

    constexpr FourVector& operator-=(const FourVector& other) {
        E -= other.E;
        px -= other.px;
        py -= other.py;
        pz -= other.pz;
        return *this;
    }

    constexpr FourVector& operator*=(double scale) {
        E *= scale;
        px *= scale;
        py *= scale;
        pz *= scale;
        return *this;
    }
};

[[nodiscard]] constexpr FourVector operator+(FourVector lhs, const FourVector& rhs) {
    lhs += rhs;
    return lhs;
}

[[nodiscard]] constexpr FourVector operator-(FourVector lhs, const FourVector& rhs) {
    lhs -= rhs;
    return lhs;
}

[[nodiscard]] constexpr FourVector operator*(FourVector vector, double scale) {
    vector *= scale;
    return vector;
}

[[nodiscard]] constexpr FourVector operator*(double scale, FourVector vector) {
    vector *= scale;
    return vector;
}

[[nodiscard]] constexpr FourVector operator/(FourVector vector, double scale) {
    vector *= 1.0 / scale;
    return vector;
}

}  // namespace lowmc
