/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/

#include "ParticleSimulation.h"

#include <algorithm>
#include <cmath>

#include "KernelFunction.h"
#include "NeighbourSearch.h"
#include "../scene/ParticleSettings.h"
#include "../utilities/SimulationStats.h"

// ____________________________________________________________________________________________________________________________________
float CURRENT_GRAVITY_SWING_ANGLE_DEG = 0.f;

// ____________________________________________________________________________________________________________________________________
float ParticleSimulation::computeDensityFactor(Particle &particle, bool isStatic)
{
    float factor;
    if (!isStatic)
    {
        float kernelSumFluid = 0;
        float kernelSumBoundary = 0;
        for (int neighbourIndex : particle.neighbours)
        {
            Particle &neighbour = vectorParticles[neighbourIndex];
            if (neighbour.type == ParticleType::Fluid)
                kernelSumFluid += KernelFunction::Function(particle, neighbour, particle.spacing);
            else if (neighbour.type == ParticleType::Boundary)
                kernelSumBoundary +=
                    KernelFunction::Function(particle, neighbour, particle.spacing);
        }
        float invertedVolume = 1.0f / PARTICLE_VOLUME;
        factor = (invertedVolume - kernelSumFluid) / kernelSumBoundary;
    }
    else if (isStatic)
    {
        factor = DENSITY_BOUNDARY_FACTOR;
    }
    return factor;
}

// ____________________________________________________________________________________________________________________________________
float ParticleSimulation::computePressureFactor(Particle &particle, bool isStatic)
{
    float factor;
    if (!isStatic)
    {
        V2f gradientSumFluid = V2f(0.f, 0.f);
        V2f gradientSumBoundary = V2f(0.f, 0.f);
        for (int neighbourIndex : particle.neighbours)
        {
            Particle &neighbour = vectorParticles[neighbourIndex];
            if (neighbour.type == ParticleType::Fluid)
                gradientSumFluid += KernelFunction::Gradient(particle, neighbour, particle.spacing);
            else if (neighbour.type == ParticleType::Boundary)
                gradientSumBoundary +=
                    KernelFunction::Gradient(particle, neighbour, particle.spacing);
        }
        float factor = (gradientSumFluid.x * gradientSumBoundary.x +
                        gradientSumFluid.y * gradientSumBoundary.y) /
                       gradientSumBoundary.lengthSquared();
    }
    else if (isStatic)
    {
        factor = PRESSURE_BOUNDARY_FACTOR;
    }
    factor = PRESSURE_BOUNDARY_FACTOR;
    return factor;
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::resetForces(Particle &particle) { particle.force = V2f(0.0f, 0.0f); }

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::resetAccelerations(Particle &particle)
{
    particle.accelerationExternal = V2f(0.0f, 0.0f);
    particle.accelerationPressure = V2f(0.0f, 0.0f);
    particle.accelerationViscosity = V2f(0.0f, 0.0f);
}


// ____________________________________________________________________________________________________________________________________
V2f ParticleSimulation::ComputeSwingingGravity(float maxAngleDeg, float periodSeconds)
{
    // Single-threaded static clock - see the header doc: this must be called exactly once per
    // step, never from inside a #pragma omp parallel for over particles.
    static float elapsedTime = 0.f;
    elapsedTime += TIMESTEP;

    constexpr float kPi = 3.14159265358979323846f;
    float maxAngleRad = maxAngleDeg * kPi / 180.f;
    // A sine wave naturally does exactly the back-and-forth the tsunami scene wants (0 -> +max
    // -> 0 -> -max -> 0 -> ...), and unlike a linear ramp that instantly reverses direction at
    // +/-max, its angular velocity eases to zero at the peaks instead of snapping - no sudden
    // jolt in gravity direction (and therefore in the fluid's forcing) at the turnaround.
    float angle = maxAngleRad * std::sin(2.f * kPi * elapsedTime / periodSeconds);
    CURRENT_GRAVITY_SWING_ANGLE_DEG = angle * 180.f / kPi;

    const float g = 9.81f;
    // angle=0 -> (0, g) straight down (this simulation's y-down convention); tilting by
    // `angle` rotates that vector toward +/-x.
    return V2f(g * std::sin(angle), g * std::cos(angle));
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::applyExternalForces(Particle &particle, V2f gravity)
{
    // F_g = m * g;
    resetForces(particle);
    V2f externalForce = gravity * particle.mass;
    particle.force += externalForce;
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::computeExternalAcceleration(Particle &particle)
{
    particle.accelerationExternal += particle.force / particle.mass;  // a = F / m.
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::computeSurfaceTension(Particle &particle, int numLayers)
{
    V2f surfaceTension = V2f(0.f, 0.f);
    // Sum over weighted mass nof neigbhours particles multiplied by a factor - k_sf / m_i
    for (int neighbourIndex : particle.neighbours)
    {
        Particle &neighbour = vectorParticles[neighbourIndex];
        if (neighbour.type == ParticleType::Fluid)
        {
            float kernelFunction = KernelFunction::Function(particle, neighbour, particle.spacing);
            V2f distance = particle.position - neighbour.position;
            surfaceTension += distance * kernelFunction * neighbour.mass;
        }
    }
    particle.accelerationExternal += surfaceTension * -SURFACE_TENSION_COEFFICIENT / particle.mass;
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::computeDensity(Particle &particle, int numLayers)
{
    float densityAcc = 0.f;
    float boundaryFactor;
    for (int neighbourIndex : particle.neighbours)
    {
        Particle &neighbour = vectorParticles[neighbourIndex];
        if (neighbour.type == ParticleType::Fluid)
            densityAcc +=
                KernelFunction::Function(particle, neighbour, particle.spacing) * neighbour.mass;
        else if (neighbour.type == ParticleType::Boundary)
            densityAcc += KernelFunction::Function(particle, neighbour, particle.spacing) *
                          neighbour.mass * DENSITY_BOUNDARY_FACTOR;
    }
    particle.density = densityAcc;
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::computePressure(Particle &particle)
{
    // If the density deviation is negative, don't generate pressure.
    float densityDev = particle.density - REST_DENSITY;
    particle.pressure = densityDev < 0 ? 0 : GAS_STIFFNESS * densityDev;
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::computePressureAcceleration(Particle &particle, int numLayers)
{
    V2f accelerationAcc = V2f(0, 0);
    V2f kernelGradient;
    for (int neighbourIndex : particle.neighbours)
    {
        Particle &neighbour = vectorParticles[neighbourIndex];
        kernelGradient = KernelFunction::Gradient(particle, neighbour, particle.spacing);
        float pressureTermParticle =  // p_i / rho_i^2
            (particle.pressure / (particle.density * particle.density));
        float pressureTermNeighbour =
            (neighbour.pressure / (neighbour.density * neighbour.density));  // p_j / rho_j^2
        float pressureTermMirrored =
            (particle.pressure / (neighbour.density * neighbour.density));  // p_i / rho_j^2.
        if (neighbour.type == ParticleType::Fluid)
        {
            // Addition to not let them go through the boundary particles.
            accelerationAcc +=
                kernelGradient * (neighbour.mass * (pressureTermParticle + pressureTermNeighbour));
        }
        else if (neighbour.type == ParticleType::Boundary)
        {
            // Single boundary layer (no scaling for now)
            accelerationAcc += kernelGradient * (PRESSURE_BOUNDARY_FACTOR * neighbour.mass *
                                                 (pressureTermParticle + pressureTermMirrored));
        }
    }
    particle.accelerationPressure = -accelerationAcc;
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::computeViscosityAcceleration(Particle &particle)
{
    V2f accelerationAcc = V2f(0, 0);
    V2f kernelGradient;
    for (int neighbourIndex : particle.neighbours)
    {
        Particle &neighbour = vectorParticles[neighbourIndex];
        kernelGradient = KernelFunction::Gradient(particle, neighbour, particle.spacing);
        float sphFactor = neighbour.mass / neighbour.density;           // m / rho
        V2f velocityRelative = particle.velocity - neighbour.velocity;  // v_ij
        V2f displacement = particle.position - neighbour.position;      // x_ij
        float num = displacement.x * velocityRelative.x +
                    displacement.y * velocityRelative.y;  // x_ij . v_ij
        float denom = displacement.x * displacement.x + displacement.y * displacement.y +
                      (0.01 * TIMESTEP * TIMESTEP);  // |x_ij|^2 + 0.01h^2
        accelerationAcc += kernelGradient * (sphFactor * (num / denom));
    }
    accelerationAcc = accelerationAcc * 2.f * VISCOSITY;
    particle.accelerationViscosity = accelerationAcc;
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::computeTotalAcceleration(Particle &particle)
{
    particle.accelerationTotal = particle.accelerationExternal + particle.accelerationPressure +
                                 particle.accelerationViscosity;
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::adaptTimeStep(float cflNumber, VectorParticles &vectorParticles,
                                       float dtMax, float dtMin)
{
    // Adaptive timestep control based on CFL condition (lambda = 0.8)
    float maxSpeed = SimulationStats::findMaxSpeed(vectorParticles);
    float spacing = vectorParticles[0].spacing;
    if (maxSpeed > 1e-6f)  // Avoid division by zero
    {
        float dtSafe = cflNumber * spacing / maxSpeed;
        dtSafe = std::clamp(dtSafe, dtMin, dtMax);
        SetTimeStep(dtSafe);
    }
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::integrate(Particle &particle)
{
    particle.velocity += particle.accelerationTotal * TIMESTEP;
    particle.position += particle.velocity * TIMESTEP;
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::isParticleInContainer(Particle &particle, std::vector<float> minsAndMaxs)
{
    float minX = minsAndMaxs[0];
    float maxX = minsAndMaxs[1];
    float minY = minsAndMaxs[2];
    float maxY = minsAndMaxs[3];
    float x = particle.position.x;
    float y = particle.position.y;
    if (x < minX || x >= maxX || y < minY || y >= maxY)
    {
        particle.type = ParticleType::Boundary;
    }
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::omitLostParticles(VectorParticles &vectorParticles,
                                           std::vector<float> minsAndMaxs, int numLayers)
{
    float spacing = vectorParticles[0].spacing;
    float offset = numLayers * spacing * 0.f;
    float minX = minsAndMaxs[0] + offset;
    float maxX = minsAndMaxs[1] - offset;
    float minY = minsAndMaxs[2] + offset;
    float maxY = minsAndMaxs[3] - offset;

#pragma omp parallel for
    for (size_t i = 0; i < vectorParticles.size(); ++i)
    {
        Particle &particle = vectorParticles[i];
        if (particle.type == ParticleType::Fluid)
        {
            float x = particle.position.x;
            float y = particle.position.y;

            if (x <= minX || x >= maxX || y <= minY || y >= maxY)
            {
                particle.type = ParticleType::Inactive;
                particle.updateColor(sf::Color(255, 255, 255, 0));
                NUM_PARTICLES_FLUID_CURRENT--;
            }
        }
    }
}

void ParticleSimulation::spatialHash(SpatialHashing &spatialHashing,
                                     VectorParticles &vectorParticles)
{
    spatialHashing.construct(vectorParticles);
    spatialHashing.query(vectorParticles);
}

void ParticleSimulation::compactHash(CompactHashing &compactHashing,
                                     VectorParticles &vectorParticles)
{
    compactHashing.construct(vectorParticles);
    compactHashing.query(vectorParticles);
}

void ParticleSimulation::indexSorting(IndexSort &indexSort, VectorParticles &vectorParticles,
                                      SortingMethod sortingMethod)
{
    indexSort.construct(vectorParticles, IndexingScheme::LinearCurve, sortingMethod);
    indexSort.query(vectorParticles, IndexingScheme::LinearCurve);
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::zIndexSorting(IndexSort &zIndexSort, VectorParticles &vectorParticles,
                                       SortingMethod sortingMethod)
{
    zIndexSort.construct(vectorParticles, IndexingScheme::ZCurve, sortingMethod);
    zIndexSort.query(vectorParticles, IndexingScheme::ZCurve);
}

void ParticleSimulation::quadraticSearch(QuadraticNeighbourSearch &quadraticSearch,
                                         VectorParticles &vectorParticles)
{
    quadraticSearch.query(vectorParticles);
}

// ____________________________________________________________________________________________________________________________________
void ParticleSimulation::simulationStep(VectorParticles &vectorParticles,
                                        CompactHashing &compactHashing,
                                        SpatialHashing &spatialHashing, IndexSort &indexSort,
                                        IndexSort &zIndexSort, std::vector<float> minsAndMaxs,
                                        int numLayers, int neighbourSearch, bool swingGravity)
{
    // Computed once, single-threaded, before the per-particle parallel loop below - see
    // ComputeSwingingGravity's doc comment for why it must never be called per-particle.
    V2f gravity = swingGravity ? ComputeSwingingGravity() : V2f(0.f, 9.81f);

// Clear all forces and acceleration.
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(vectorParticles.size()); ++i)
    {
        Particle &particle = vectorParticles[i];
        if (particle.type == ParticleType::Fluid)
        {
            applyExternalForces(particle, gravity);
            resetAccelerations(particle);
        }
    }

// Compute density and pressure for fluid particles.
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(vectorParticles.size()); ++i)
    {
        Particle &particle = vectorParticles[i];
        if (particle.type == ParticleType::Fluid)
        {
            computeDensity(particle, numLayers);
            computePressure(particle);
        }
    }

// Compute accelerations.
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(vectorParticles.size()); ++i)
    {
        Particle &particle = vectorParticles[i];
        if (particle.type == ParticleType::Fluid)
        {
            computeExternalAcceleration(particle);
            computePressureAcceleration(particle, numLayers);
            computeViscosityAcceleration(particle);
        }
    }

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(vectorParticles.size()); ++i)
    {
        Particle &particle = vectorParticles[i];
        if (particle.type == ParticleType::Fluid)
        {
            computeTotalAcceleration(particle);
        }
    }

// Perform an integration step.
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(vectorParticles.size()); ++i)
    {
        Particle &particle = vectorParticles[i];
        if (particle.type == ParticleType::Fluid) integrate(particle);
    }

    // Omit the lost particles with a specifc world size.
    omitLostParticles(vectorParticles, minsAndMaxs, numLayers);
}