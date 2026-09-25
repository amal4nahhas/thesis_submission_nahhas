/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#include "KernelFunction.h"

#include <vector>

#include "../utilities/MathUtils.h"

// Define the Pi value because it is not available in c++17.
constexpr float PI = 3.14159265358979323846f;

// ____________________________________________________________________________________________________________________________________
float KernelFunction::Function(Particle &a, Particle &b, float h)
{
    float kernelValue = 0.0;
    float alpha = 5.0 / (14.0 * PI * h * h);       // 2D normalization factor
    float dist = Particle::computeDistance(b, a);  // reversed distance for the convolution.
    float q = dist / h;
    float term1 = std::pow(2 - q, 3);
    float term2 = -4 * std::pow(1 - q, 3);
    if (q < 2.f)
    {
        kernelValue += term1;
        if (q < 1.f) kernelValue += term2;
        kernelValue *= alpha;
    }
    return kernelValue;
}

// ____________________________________________________________________________________________________________________________________
V2f KernelFunction::Gradient(const Particle &a, const Particle &b, float h)
{
    V2f kernelGradient(0.0f, 0.0f);
    float alpha = 5.f / (14.f * PI * h * h);       // 2D normalization factor.
    float dist = Particle::computeDistance(b, a);  // reversed distance for the convolution.
    V2f displacement(a.position - b.position);
    float q = dist / h;
    V2f nablaQ = !Math::areFloatsEqual(dist, 0.f)
                     ? (displacement / (dist * h))
                     : V2f(0.0f, 0.0f);  // gradient of the normalized distance q
    float term1 = -3 * std::pow(2 - q, 2);
    float term2 = 12 * std::pow(1 - q, 2);
    if (Math::areFloatsEqual(q, 0.0f))
        return kernelGradient;  // Return zero vector if particles are at the same position.
    else if (q >= 2.f)
        return kernelGradient;  // adjacent particle is out of local support.
    else if (q >= 1.f)
        kernelGradient = nablaQ * alpha * term1;
    else
        kernelGradient = nablaQ * alpha * (term1 + term2);
    // Correct negative zeros caused by float computations.
    if (kernelGradient.x == -0.0f) kernelGradient.x = 0.0f;
    if (kernelGradient.y == -0.0f) kernelGradient.y = 0.0f;
    return kernelGradient;
}
