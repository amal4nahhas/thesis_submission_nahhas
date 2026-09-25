/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/

#include "MathUtils.h"

#include <cmath>

bool Math::isPrime(unsigned int n)
{
    if (n < 2) return false;
    for (int i = 2; i * i <= static_cast<int>(n); i++)
    {
        if (n % i == 0) return false;
    }
    return true;
}

unsigned int Math::nextPrime(unsigned int m)
{
    auto p = m;
    while (!isPrime(p)) p++;
    return p;
}

bool Math::areFloatsEqual(float a, float b, float epsilon) { return std::fabs(a - b) < epsilon; }

bool Math::areVector2fEqual(const V2f &a, const V2f &b, float epsilon)
{
    return areFloatsEqual(a.x, b.x, epsilon) && areFloatsEqual(a.y, b.y, epsilon);
}

std::uint8_t Math::lerp8(std::uint8_t c0, std::uint8_t c1, float t)
{
    float c = c0 + (c1 - c0) * t;
    return static_cast<std::uint8_t>(c + 0.5f);
}

sf::Color Math::lerpColor(sf::Color color0, sf::Color color1, float t)
{
    sf::Color result;
    result.r = lerp8(color0.r, color1.r, t);
    result.g = lerp8(color0.g, color1.g, t);
    result.b = lerp8(color0.b, color1.b, t);
    result.a = 255;
    return result;
}
