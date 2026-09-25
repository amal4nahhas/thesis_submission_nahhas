/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include <string>

#include "../simulation/Particle.h"

// Variables related to container.
// These numbers state how many particles are needed to fully fill the container with
extern unsigned int CONTAINER_SIZE_X;
extern unsigned int CONTAINER_SIZE_Y;

// Various variables to describe the number of particles in the simulation.
// Variables related to the fluid particles.
extern unsigned int NUM_FLUID_PER_ROW;
extern unsigned int NUM_FLUID_PER_COLUMN;
extern unsigned int NUM_PARTICLES_FLUID;
extern unsigned int NUM_PARTICLES_FLUID_CURRENT;

// Variables related to the boundary particles.
extern unsigned int NUM_BOUNDARY_PER_ROW;
extern unsigned int NUM_BOUNDARY_PER_COLUMN;
extern unsigned int NUM_PARTICLES_BOUNDARY;

// Time step, gazz stiffnes and viscosity paramteres.
extern float TIMESTEP;
extern float GAS_STIFFNESS;
extern float VISCOSITY;

// Settings of the fluid simulation.
extern float REST_DENSITY;
extern float SPACING;
extern float PARTICLE_VOLUME;
extern float FLUID_VOLUME;
extern float PARTICLE_MASS;
extern float DENSITY_BOUNDARY_FACTOR;
extern float PRESSURE_BOUNDARY_FACTOR;
extern float SURFACE_TENSION_COEFFICIENT;

// An array to store all particles.
extern VectorParticles vectorParticles;

class Scene
{
   public:
    unsigned int numParticlesInitF;
    unsigned int numParticlesCurrentF;  // In case particles get omitted.
    unsigned int numParticlesB;
};

class RectalScene : Scene
{
   public:
    // Fluid Block
    unsigned int numParticlesPerRowF;
    unsigned int numParticlesPerColumnF;
    // Boundary Frame
    unsigned int numParticlesPerRowB;
    unsigned int numParticlesPerColumnB;

    // Simulation parameters
    float timestep;
    float gasStiffness;
    float viscosity;
};

// ************************************************************************************************************************************
// Setting and computing global variables:
// ************************************************************************************************************************************

/// @brief Set the container size in number of particles per axis.
/// @param containerSizeX
/// @param containerSizeY
void SetContainerSize(unsigned int containerSizeX, unsigned int containerSizeY);

/// @brief Set the number of fluid per row and column.
/// @param numFluidPerRow
/// @param numFluidPerColumn
void SetNumFluidParticles(unsigned int numFluidPerRow, unsigned int numFluidPerColumn);

/// @brief Compute the number of boundary particles per row and col using the container size.
void ComputeNumBoundaryParticles(int numLayers = 1);

/// @brief Compute the fluid volume from the number of fluid particles and the predefined particle
/// volume. This is for the sake of integrity.
void ComputeFluidVolume();

/// @brief Compute the particles mass from the volume and the rest density.
void ComputeParticleMass();

/// @brief Compute the spacing using the particle's volume.
void ComputeSpacing();

/// @brief Set the time step.
/// @param timeStep
void SetTimeStep(float timeStep);

/// @brief Set the gazz stiffnes.
/// @param gatStiffnes
void SetGasStiffnes(float gatStiffnes);

/// @brief Set the viscoisty coeffecient.
/// @param viscosity
void SetViscosity(float viscosity);

/// @brief Set the rest density of the fluid.
/// @param restDensity
void SetRestDensity(float restDensity = 1111.f);

// ************************************************************************************************************************************
// Instantiation:
// ************************************************************************************************************************************

/// @brief Instantiate fluid particles.
/// @param vectorParticles
void InstantiateFluidParticles(VectorParticles &vectorParticles);

/// @brief Instantiate the boundary particles.
/// @param vectorParticles
void InstantiateBoundaryParticles(VectorParticles &vectorParticles);

/// @brief Instantaite one fluid particle.
/// @param vectorParticles
Particle &InstantiateOneParticle(VectorParticles &vectorParticles);

/// @brief Instantiate one boundary line with length in number of particles.
/// @param vectorParticles
/// @param length
void InstantiateBoundaryLine(VectorParticles &vectorParticles, int length, int numLayers = 1);

/// @brief Instantiate fluid and boundary particles.
/// @param vectorParticles
void InstantiateParticles(VectorParticles &vectorParticles);

// ************************************************************************************************************************************
// Positioning:
// ************************************************************************************************************************************

/// @brief Position particle with the given position.
/// @param particle
/// @param position
void PositionParticle(Particle &particle, float x, float y);

/// @brief Position fluid particles locally with a fix spacing.
/// @param vectorParticles
void PositionFluidParticlesLocal(VectorParticles &vectorParticles);

/// @brief Position the boundary particles locally.
/// @param vectorParticles
void PositionBoundaryParticlesLocal(VectorParticles &vectorParticles, int numLayers = 1);

/// @brief Position fluid and boundary particles.
/// @param vectorParticles
void PositionParticles(VectorParticles &vectorParticles, int numLayers = 1);

/// @brief Translate the fluid particles by the given values.
/// @param vectorParticles
/// @param dx
/// @param dy
void TranslateParticlesFluid(VectorParticles &vectorParticles, float dx, float dy);

/// @brief Rotate the fluid particles grid by the given angle in degrees.
/// @param vectorParticles
/// @param angleDegrees
void RotateParticles(VectorParticles &vectorParticles, float angleDegrees);

void PositionLineLocal(VectorParticles &vectorParticles, int numLayers = 1);

// ************************************************************************************************************************************
// Initialization of different scenarios:
// ************************************************************************************************************************************

/// @brief Set up the scene with container size, number of fluid particles and rest density.
/// density.
/// @param vectorParticles
/// @param containerSizeX
/// @param containerSizeY
/// @param numFluidParticlesRow
/// @param numFluidParticlesColumn
/// @param restDensity
void SetupScene(VectorParticles &vectorParticles, unsigned int containerSizeX,
                unsigned int containerSizeY, int numFluidParticlesRow, int numFluidParticlesColumn,
                int numLayers = 1, float restDensity = 1111.f);

/// @brief Initialize a water column.
/// @param vectorParticles
/// @param numLayers
void InitializeWaterColumn(VectorParticles &vectorParticles, int numLayers = 1);

/// @brief Initialize a water column for the sake of testing.
/// @param vectorParticles
/// @param numLayers
/// @param containerSizeX
/// @param containerSizeY
void InitializeStaticWaterColumn(VectorParticles &vectorParticles, int numLayers,
                                 int containerSizeX, int containerSizeY);

/// @brief Initialize the breaking dam.
/// @param vectorParticles
/// @param numLayers
void InitializeBreakingDamSimple(VectorParticles &vectorParticles, int numLayers = 1, int run = 1);

/// @brief Initialize a wide, low breaking-dam scene (box width:height = 10:3) mostly filled
/// with fluid (occupying 8:2 of the box on the same unit scale, left-anchored, resting on the
/// floor, not filled to the top) - exactly 1M fluid particles. Used for the thesis's large-N
/// dynamic-scene benchmark.
/// @param vectorParticles
/// @param numLayers
void InitializeBreakingDamWideLow(VectorParticles &vectorParticles, int numLayers = 1);

/// @brief Same wide-low breaking dam scene as above, but with the box/fluid unit size u given
/// explicitly instead of hardcoded to 250 (which lands on exactly 1M fluid particles) - lets a
/// dynamic-scaling benchmark reuse the same scene shape at other particle counts. Fluid particle
/// count is 16*u*u, so u = round(sqrt(targetFluidCount / 16.0)) to hit a given target.
/// @param vectorParticles
/// @param numLayers
/// @param u
void InitializeBreakingDamWideLow(VectorParticles &vectorParticles, int numLayers, int u);

/// @brief Initialize the fluid simulation with the given scene number.
/// @param vectorParticles
/// @param numLayers
/// @param scene
void InitializeFluidSimulation(VectorParticles &vectorParticles, int numLayers, int scene);

/// @brief Initialize a scene from a position file produced by tools/scene_import (one
/// "x y type" line per particle, type in {fluid, boundary}). Bypasses the rectangular
/// row/col scene setup used elsewhere, since a hand-drawn scene isn't a rectangle.
/// @param vectorParticles
/// @param path path to the position file, resolved relative to the working directory.
/// @param numLayers only used to size the neighbour-search grid's cell-reference table;
///        does not thicken the boundary itself — the file is trusted as-is.
/// @param restDensity rest density for the loaded fluid.
/// @param timeStep simulation time step.
/// @param gasStiffness gas stiffness constant (k).
/// @param viscosity viscosity coefficient.
void InitializeFromSceneFile(VectorParticles &vectorParticles, const std::string &path,
                             int numLayers = 1, float restDensity = 1111.f,
                             float timeStep = 0.001f, float gasStiffness = 7500.f,
                             float viscosity = 0.1f);
