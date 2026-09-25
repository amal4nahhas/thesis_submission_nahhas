/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include "NeighbourSearch.h"
#include "Particle.h"

extern const unsigned int WINDOW_WIDTH;
extern const unsigned int WINDOW_HEIGHT;

/// Current gravity swing angle in degrees (0 = straight down), updated once per call to
/// ComputeSwingingGravity; read this instead of calling that function again, since calling it
/// advances the swing's internal clock (see its own doc comment).
extern float CURRENT_GRAVITY_SWING_ANGLE_DEG;

class ParticleSimulation
{
   public:
    /// @brief Compute the correcting coeffecient for the density computation.
    /// @param particle
    /// @return
    float computeDensityFactor(Particle &particle, bool isStatic = true);

    /// @brief Compute the correcting coeffecient for the pressure acceleration computation.
    /// @param particle
    /// @return
    float computePressureFactor(Particle &particle, bool isStatic = true);

    /// @brief Reset the forces to zero in all directions.
    /// @param particle
    void resetForces(Particle &particle);

    /// @brief Reset the acceleration in all directioons.
    /// @param particle
    void resetAccelerations(Particle &particle);

    /// @brief Add external or gravitational forces (one earth) to each particle.
    /// @param particle
    /// @param gravity current gravity vector (see ComputeSwingingGravity for the tsunami
    /// scene's oscillating tilt; every other scene just passes straight-down gravity).
    void applyExternalForces(Particle &particle, V2f gravity);

    /// @brief Compute external acceleration.
    /// @param particle
    void computeExternalAcceleration(Particle &particle);

    void computeSurfaceTension(Particle &particle, int numLayers);

    /// @brief Compute the density for a given particle.
    /// @param particle
    void computeDensity(Particle &particle, int numLayers);

    /// @brief Compute pressure for a given particles.
    /// @param particle
    void computePressure(Particle &particle);

    /// @brief Compute the pressure acceleration.
    /// @param particle
    void computePressureAcceleration(Particle &particle, int numLayers);

    /// @brief Compute viscosity acceleration.
    /// @param particle
    void computeViscosityAcceleration(Particle &particle);

    /// @brief Compute the total acceleration by summing all acceleration contributions.
    /// @param particle
    void computeTotalAcceleration(Particle &particle);

    /// @brief Adapt the time step based on the CFL condition, the maximum velocity and the maximum
    /// acceleration
    /// @param cflNumber
    /// @param vectorParticles
    /// @param dtMax
    /// @param dtMin
    void adaptTimeStep(float cflNumber, VectorParticles &vectorParticles, float dtMax, float dtMin);

    /// @brief Perform an integration step using the explicit Euler method.
    /// @param particle
    void integrate(Particle &particle);

    /// @brief Check whether a particle is within the container defined by the boundary particles.
    /// @param particle
    /// @param minsAndMaxs
    void isParticleInContainer(Particle &particle, std::vector<float> minsAndMaxs);

    /// @brief Omit particles that pass the boundary particles and float freely, assuming that the
    /// container is not rotated.
    /// @param vectorParticles
    /// @param minsAndMaxs
    void omitLostParticles(VectorParticles &vectorParticles, std::vector<float> minsAndMaxs,
                           int numLayers);

    void spatialHash(SpatialHashing &spatialHashing, VectorParticles &vectorParticles);

    void compactHash(CompactHashing &compactHashing, VectorParticles &vectorParticles);

    void indexSorting(IndexSort &indexSort, VectorParticles &vectorParticles,
                      SortingMethod sortingMethod = SortingMethod::InsertionSort);

    void zIndexSorting(IndexSort &zIndexSort, VectorParticles &vectorParticles,
                       SortingMethod sortingMethod = SortingMethod::RadixSort);

    void quadraticSearch(QuadraticNeighbourSearch &quadraticSearch,
                         VectorParticles &vectorParticles);

    /// @brief Perform a simulation step.
    /// @param vectorParticles
    /// @param swingGravity when true, tilts gravity back and forth like a pendulum instead of
    /// straight down - see ComputeSwingingGravity. Used to keep otherwise-static-looking
    /// scenes (e.g. the tsunami scene) visually active.
    void simulationStep(VectorParticles &vectorParticles, CompactHashing &compactHashing,
                        SpatialHashing &spatialHashing, IndexSort &indexSort, IndexSort &zIndexSort,
                        std::vector<float> minsAndMaxs, int numLayers, int neighbourSearch = 1,
                        bool swingGravity = false);

    /// @brief Gravity vector that oscillates like a pendulum between +/-maxAngleDeg from
    /// straight down, completing one full swing (0 -> +max -> 0 -> -max -> 0) every
    /// periodSeconds. Call exactly once per simulation step (not per particle): it advances
    /// its own internal clock by TIMESTEP each call, so calling it per-particle would
    /// fast-forward the swing and race across threads.
    /// @param maxAngleDeg swing amplitude in degrees from straight down.
    /// @param periodSeconds time for one full swing cycle - SMALLER is FASTER.
    static V2f ComputeSwingingGravity(float maxAngleDeg = 30.f, float periodSeconds = 6.f);
};