/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include "../simulation/Particle.h"

class Math
{
   public:
    static bool isPrime(unsigned int n);

    static unsigned int nextPrime(unsigned int m);

    static bool areFloatsEqual(float a, float b, float epsilon = 0.0001f);

    static bool areVector2fEqual(const V2f &a, const V2f &b, float epsilon = 0.0001f);

    static std::uint8_t lerp8(std::uint8_t a, std::uint8_t b, float t);

    static sf::Color lerpColor(sf::Color color0, sf::Color color1, float t);
};
