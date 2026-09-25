/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once
#include "Particle.h"

class KernelFunction
{
   public:
    /// @brief Compute the kernel function value between two particles.
    /// @param a
    /// @param b
    /// @param h
    /// @return The kernel function value between the two particles.
    static float Function(Particle &a, Particle &b, float h);

    /// @brief Compute the gradient of the kernel function between two particles.
    /// @param a
    /// @param b
    /// @param h
    /// @return
    static V2f Gradient(const Particle &a, const Particle &b, float h);
};