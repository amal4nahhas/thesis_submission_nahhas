/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#include "ParticleSettings.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "../simulation/Particle.h"

// ____________________________________________________________________________________________________________________________________
unsigned int NUM_FLUID_PER_ROW;
unsigned int NUM_FLUID_PER_COLUMN;
unsigned int NUM_PARTICLES_FLUID;
unsigned int NUM_PARTICLES_FLUID_CURRENT;

// ____________________________________________________________________________________________________________________________________
unsigned int CONTAINER_SIZE_X;
unsigned int CONTAINER_SIZE_Y;

// ____________________________________________________________________________________________________________________________________
unsigned int NUM_BOUNDARY_PER_ROW;
unsigned int NUM_BOUNDARY_PER_COLUMN;
unsigned int NUM_PARTICLES_BOUNDARY;

// ____________________________________________________________________________________________________________________________________
float TIMESTEP;
float GAS_STIFFNESS;
float VISCOSITY;

// ____________________________________________________________________________________________________________________________________
float REST_DENSITY;
float PARTICLE_VOLUME = 0.01f;  // [mm^2]
float FLUID_VOLUME;
float SPACING;
float PARTICLE_MASS;

float DENSITY_BOUNDARY_FACTOR = 1.2;
float PRESSURE_BOUNDARY_FACTOR = 1.1f;

float SURFACE_TENSION_COEFFICIENT = 0.8f;

// ____________________________________________________________________________________________________________________________________
VectorParticles vectorParticles;

// ************************************************************************************************************************************
// Set global variables:
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
void SetNumFluidParticles(unsigned int numFluidPerRow, unsigned int numFluidPerColumn)
{
    NUM_FLUID_PER_ROW = numFluidPerRow;
    NUM_FLUID_PER_COLUMN = numFluidPerColumn;
    NUM_PARTICLES_FLUID = NUM_FLUID_PER_ROW * NUM_FLUID_PER_COLUMN;
    NUM_PARTICLES_FLUID_CURRENT = NUM_PARTICLES_FLUID;
}

// ____________________________________________________________________________________________________________________________________
void SetContainerSize(unsigned int containerSizeX, unsigned int containerSizeY)
{
    CONTAINER_SIZE_X = containerSizeX;
    CONTAINER_SIZE_Y = containerSizeY;
}

// ____________________________________________________________________________________________________________________________________
void ComputeNumBoundaryParticles(int numLayers)
{
    NUM_BOUNDARY_PER_ROW = CONTAINER_SIZE_X + 2 * numLayers;
    NUM_BOUNDARY_PER_COLUMN = CONTAINER_SIZE_Y + 2 * numLayers;
    NUM_PARTICLES_BOUNDARY =
        NUM_BOUNDARY_PER_ROW * NUM_BOUNDARY_PER_COLUMN - CONTAINER_SIZE_Y * CONTAINER_SIZE_X;
}

// ____________________________________________________________________________________________________________________________________
void ComputeFluidVolume()
{
    FLUID_VOLUME = NUM_BOUNDARY_PER_COLUMN * NUM_BOUNDARY_PER_ROW * PARTICLE_VOLUME;
}

// ____________________________________________________________________________________________________________________________________
void ComputeSpacing() { SPACING = sqrt(PARTICLE_VOLUME); }

// ____________________________________________________________________________________________________________________________________
void ComputeParticleMass() { PARTICLE_MASS = REST_DENSITY * PARTICLE_VOLUME; }

// ____________________________________________________________________________________________________________________________________
void SetTimeStep(float timeStep) { TIMESTEP = timeStep; }

// ____________________________________________________________________________________________________________________________________
void SetGasStiffnes(float gatStiffnes) { GAS_STIFFNESS = gatStiffnes; }

// ____________________________________________________________________________________________________________________________________
void SetViscosity(float viscosity) { VISCOSITY = viscosity; }

// ____________________________________________________________________________________________________________________________________
void SetRestDensity(float restDensity) { REST_DENSITY = restDensity; }

// ************************************************************************************************************************************
// Instantiation:
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
void InstantiateFluidParticles(VectorParticles &vectorParticles)
{
    for (int i = 0; i < NUM_PARTICLES_FLUID; i++)
    {
        Particle particle(PARTICLE_MASS, REST_DENSITY, SPACING, ParticleType::Fluid, i);
        vectorParticles.push_back(particle);
    }
}

// ____________________________________________________________________________________________________________________________________
void InstantiateBoundaryParticles(VectorParticles &vectorParticles)
{
    for (int i = 0; i < NUM_PARTICLES_BOUNDARY; i++)
    {
        unsigned id = NUM_PARTICLES_FLUID + i;
        Particle particle(PARTICLE_MASS, REST_DENSITY, SPACING, ParticleType::Boundary, id);
        vectorParticles.push_back(particle);
    }
}

// ____________________________________________________________________________________________________________________________________
Particle &InstantiateOneParticle(VectorParticles &vectorParticles)
{
    Particle particle(PARTICLE_MASS, REST_DENSITY, SPACING, ParticleType::Fluid, 1);
    vectorParticles.push_back(particle);
    return vectorParticles[0];
}

// ____________________________________________________________________________________________________________________________________
void InstantiateBoundaryLine(VectorParticles &vectorParticles, int length, int numLayers)
{
    for (int i = 0; i < length * numLayers; i++)
    {
        Particle particle(PARTICLE_MASS, REST_DENSITY, SPACING, ParticleType::Boundary, vectorParticles.size() + i);
        vectorParticles.push_back(particle);
    }
}

// ____________________________________________________________________________________________________________________________________
void InstantiateParticles(VectorParticles &vectorParticles)
{
    InstantiateFluidParticles(vectorParticles);
    InstantiateBoundaryParticles(vectorParticles);
}

// ************************************************************************************************************************************
// Positioning:
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
void PositionParticle(Particle &particle, float x, float y) { particle.setPosition(x, y); }

// ____________________________________________________________________________________________________________________________________
void PositionFluidParticlesLocal(VectorParticles &vectorParticles)
{
    float offset = SPACING / 2.f;
    unsigned int particleIndex = 0;
    int numRows = NUM_FLUID_PER_COLUMN;
    int numCols = NUM_FLUID_PER_ROW;
    for (int row = 0; row < numRows; row++)
    {
        for (int col = 0; col < numCols; col++)
        {
            float x, y;
            x = col * SPACING + offset;
            y = row * SPACING + offset;
            V2f position = V2f(x, y);
            Particle &particle = vectorParticles[particleIndex];
            particle.position = position;
            particleIndex++;
        }
    }
    printf("√ Positioned %u fluid particles in a grid.\n", NUM_PARTICLES_FLUID);
}

// ____________________________________________________________________________________________________________________________________
void PositionBoundaryParticlesLocal(VectorParticles &vectorParticles, int numLayers)
{
    float offset = SPACING / 2.f;
    unsigned int particleIndex = NUM_PARTICLES_FLUID;
    int numRows = NUM_BOUNDARY_PER_COLUMN;
    int numCols = NUM_BOUNDARY_PER_ROW;
    for (int row = 0; row < numRows; row++)
    {
        for (int col = 0; col < numCols; col++)
        {
            // Skip the inner fluid particles.
            if (row >= numLayers && row < numRows - numLayers && col >= numLayers &&
                col < numCols - numLayers)
            {
                continue;
            }
            float x, y;
            x = col * SPACING + offset;
            y = row * SPACING + offset;
            V2f position = V2f(x, y);
            Particle &particle = vectorParticles[particleIndex];
            particle.position = position;
            particleIndex++;
        }
    }
}

// ____________________________________________________________________________________________________________________________________
void PositionLineLocal(VectorParticles &vectorParticles, int numLayers)
{
    float offset = SPACING / 2.f;
    unsigned int particleIndex = 0;
    int length = vectorParticles.size() / numLayers;
    for (int row = 0; row < numLayers; row++)
    {
        for (int i = 0; i < length; i++)
        {
            float x, y;
            x = i * SPACING + offset;
            y = row * SPACING + offset;
            V2f position = V2f(x, y);
            Particle &particle = vectorParticles[particleIndex];
            particle.position = position;
            particleIndex++;
        }
    }
    printf("√ Positioned %d boundary particles in a line.\n", length);
}

// ____________________________________________________________________________________________________________________________________
void TranslateParticlesFluid(VectorParticles &vectorParticles, float dx, float dy)
{
    for (auto &particle : vectorParticles)
    {
        if (particle.type == ParticleType::Fluid)
        {
            particle.translate(dx, dy);
        }
    }
}

// ____________________________________________________________________________________________________________________________________
void TranslateParticles(VectorParticles &vectorParticles, float dx, float dy)
{
    for (auto &particle : vectorParticles)
    {
        particle.translate(dx, dy);
    }
}

// ____________________________________________________________________________________________________________________________________
void RotateParticles(VectorParticles &vectorParticles, float angleDegrees)
{
    for (auto &particle : vectorParticles)
    {
        particle.rotate(angleDegrees);
    }
}

// ____________________________________________________________________________________________________________________________________
void PositionParticles(VectorParticles &vectorParticles, int numLayers)
{
    PositionFluidParticlesLocal(vectorParticles);
    PositionBoundaryParticlesLocal(vectorParticles, numLayers);
    TranslateParticlesFluid(vectorParticles, SPACING * numLayers, SPACING * numLayers);
}

// ************************************************************************************************************************************
// Setup and Initialization.
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
void SetupScene(VectorParticles &vectorParticles, unsigned int containerSizeX,
                unsigned int containerSizeY, int numFluidParticlesRow, int numFluidParticlesColumn,
                int numLayers, float restDensity)
{
    SetContainerSize(containerSizeX, containerSizeY);
    SetNumFluidParticles(numFluidParticlesRow, numFluidParticlesColumn);
    ComputeNumBoundaryParticles(numLayers);

    // Set rest density.
    SetRestDensity(restDensity);
    ComputeFluidVolume();  // For the sake of integrity.

    // Compute a particle's volume, the spacing and the particle's mass.
    ComputeSpacing();
    ComputeParticleMass();
}

void InitializeScene(VectorParticles &vectorParticles, int numLayers, int containerSizeX,
                     int containerSizeY, int fluidX, int fluidY, bool empty = false)
{
    SetupScene(vectorParticles, containerSizeX, containerSizeY, fluidX, fluidY, numLayers);

    if (!empty)
    {
        InstantiateParticles(vectorParticles);
        PositionParticles(vectorParticles, numLayers);
        TranslateParticlesFluid(vectorParticles, 0, 0);
    }
    else
    {
        NUM_PARTICLES_FLUID_CURRENT = 0;
        InstantiateBoundaryParticles(vectorParticles);
        PositionBoundaryParticlesLocal(vectorParticles, numLayers);
    }

    // Set the desired k, dt, and visc here:
    float k = 250;
    float dt = 0.005f;
    float visc = 0.1f;

    SetGasStiffnes(k);
    SetTimeStep(dt);
    SetViscosity(visc);
}

// ____________________________________________________________________________________________________________________________________
void InitializeWaterColumn(VectorParticles &vectorParticles, int numLayers)
{
    int containerSizeX = 1000;
    int containerSizeY = 1100;

    int fluidX = containerSizeX;
    int fluidY = containerSizeX;

    SetupScene(vectorParticles, containerSizeX, containerSizeY, fluidX, fluidY, numLayers);

    InstantiateParticles(vectorParticles);
    PositionParticles(vectorParticles, numLayers);

    TranslateParticlesFluid(vectorParticles, 0, 100 * SPACING);

    float k0 = 500.f;
    float k1 = 100.f;
    float k2 = 1000.f;
    float k3 = 10000.f;

    // Set the desired k, dt, and visc here:
    float k = 5000.f;
    float dt = 0.0002f;
    float visc = 0.1f;

    SetGasStiffnes(k);
    SetTimeStep(dt);
    SetViscosity(visc);
}

// ____________________________________________________________________________________________________________________________________
void InitializeStaticWaterColumn(VectorParticles &vectorParticles, int numLayers,
                                 int containerSizeX, int containerSizeY)
{
    int fluidX = containerSizeX;
    int fluidY = containerSizeY;

    SetupScene(vectorParticles, containerSizeX, containerSizeY, fluidX, fluidY, numLayers);

    InstantiateParticles(vectorParticles);
    PositionParticles(vectorParticles, numLayers);

    TranslateParticlesFluid(vectorParticles, 0, 0);

    float k0 = 500.f;
    float k1 = 100.f;
    float k2 = 1000.f;
    float k3 = 10000.f;

    // Set the desired k, dt, and visc here:
    float k = k0;
    float dt = 0.005f;
    float visc = 0.1f;

    SetGasStiffnes(k);
    SetTimeStep(dt);
    SetViscosity(visc);
}

// ____________________________________________________________________________________________________________________________________
void InitializeBreakingDamSimple(VectorParticles &vectorParticles, int numLayers, int run)
{
    int factor = 1;
    if (run == 2)
        factor = 2;
    else if (run == 3)
        factor = 4;

    int containerSizeX = 120 * factor;
    int containerSizeY = 70 * factor;

    int fluidX = containerSizeX * 5 / 6;
    int fluidY = containerSizeY * factor * 5 / 6;

    SetupScene(vectorParticles, containerSizeX, containerSizeY, fluidX, fluidY, numLayers);

    InstantiateParticles(vectorParticles);
    PositionParticles(vectorParticles, numLayers);

    TranslateParticlesFluid(vectorParticles, 0, SPACING * (containerSizeY - fluidY));

    float k1 = 100.f;
    float dt11 = 0.01f;
    float dt12 = 0.02f;

    float k2 = 1000.f;
    float dt21 = 0.002f;
    float dt22 = 0.007f;

    float k3 = 250.f;
    float dt31 = 0.005f;
    float dt32 = 0.0019f;

    // Set the desired k and dt here:
    float k = 5000.f * 1.0f;
    float dt = 0.001f;
    float visc = 0.1f;

    SetGasStiffnes(k);
    SetTimeStep(dt);
    SetViscosity(visc);
}

// ____________________________________________________________________________________________________________________________________
void InitializeBreakingDamWideLow(VectorParticles &vectorParticles, int numLayers)
{
    // u=250 lands on exactly 1,000,000 fluid particles (fluidX*fluidY = 2000*500) - see the
    // parametrized overload below for the general case.
    InitializeBreakingDamWideLow(vectorParticles, numLayers, 250);
}

// ____________________________________________________________________________________________________________________________________
void InitializeBreakingDamWideLow(VectorParticles &vectorParticles, int numLayers, int u)
{
    // Box: width:height = 10:3 (wide, low). Fluid occupies 8:2 of the box on the SAME unit
    // scale as the box itself (not an independent shape ratio) - a lot of water, resting on
    // the floor, left-anchored, not filled all the way to the top (leaves 1/3 of the box
    // height as headspace and 1/5 of its width empty on the right for the dam-break motion).
    int containerSizeX = 10 * u;  // 2500
    int containerSizeY = 3 * u;   // 750

    int fluidX = 8 * u;  // 2000
    int fluidY = 2 * u;  // 500

    SetupScene(vectorParticles, containerSizeX, containerSizeY, fluidX, fluidY, numLayers);

    InstantiateParticles(vectorParticles);
    PositionParticles(vectorParticles, numLayers);

    // Fluid is placed at the top-left by PositionParticles (row 0 = top); shift it down so it
    // rests on the floor instead.
    TranslateParticlesFluid(vectorParticles, 0, SPACING * (containerSizeY - fluidY));

    float k = 7500.f;
    float dt = 0.001f;
    float visc = 0.1f;

    SetGasStiffnes(k);
    SetTimeStep(dt);
    SetViscosity(visc);
}

// ____________________________________________________________________________________________________________________________________
static void AddBoundaryParticle(VectorParticles &vectorParticles, float x, float y)
{
    Particle particle(PARTICLE_MASS, REST_DENSITY, SPACING, ParticleType::Boundary, vectorParticles.size() + 1);
    particle.setPosition(x, y);
    vectorParticles.push_back(particle);
    NUM_PARTICLES_BOUNDARY++;
}

// ____________________________________________________________________________________________________________________________________
static void PositionFluidParticlesTwoBlocks(VectorParticles &vectorParticles, int numLayers,
                                            int containerSizeX, int containerSizeY,
                                            float sideClearance)
{
    const int totalCols = static_cast<int>(NUM_FLUID_PER_ROW);
    const int totalRows = static_cast<int>(NUM_FLUID_PER_COLUMN);
    const int leftCols = totalCols / 2;
    const int rightCols = totalCols - leftCols;

    const float inset = (numLayers + 2.0f) * SPACING;
    const float blockY0 = (numLayers + containerSizeY - totalRows) * SPACING;

    const float leftX0 = inset + sideClearance;
    const float rightBlockWidth = rightCols * SPACING;
    const float rightX0 =
        (numLayers + containerSizeX) * SPACING - inset - rightBlockWidth - sideClearance;

    unsigned int particleIndex = 0;

    for (int row = 0; row < totalRows; ++row)
    {
        for (int col = 0; col < leftCols; ++col)
        {
            Particle &particle = vectorParticles[particleIndex++];
            particle.setPosition(leftX0 + (col + 0.5f) * SPACING, blockY0 + row * SPACING);
        }
    }

    for (int row = 0; row < totalRows; ++row)
    {
        for (int col = 0; col < rightCols; ++col)
        {
            Particle &particle = vectorParticles[particleIndex++];
            particle.setPosition(rightX0 + (col + 0.5f) * SPACING, blockY0 + row * SPACING);
        }
    }
}

// ____________________________________________________________________________________________________________________________________
static void PositionFluidBlock(VectorParticles &vectorParticles, int cols, int rows, float startX,
                               float startY)
{
    unsigned int particleIndex = 0;
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            Particle &particle = vectorParticles[particleIndex++];
            particle.setPosition(startX + (col + 0.5f) * SPACING, startY + (row + 0.5f) * SPACING);
        }
    }
}

// ____________________________________________________________________________________________________________________________________
static void AddBoundaryLineHorizontal(VectorParticles &vectorParticles, float xStart, float y,
                                      int count)
{
    for (int i = 0; i < count; ++i)
    {
        AddBoundaryParticle(vectorParticles, xStart + i * SPACING, y);
    }
}

// ____________________________________________________________________________________________________________________________________
static void AddBoundaryLineVertical(VectorParticles &vectorParticles, float x, float yStart,
                                    int count)
{
    for (int i = 0; i < count; ++i)
    {
        AddBoundaryParticle(vectorParticles, x, yStart + i * SPACING);
    }
}

// ___________________________________________________________________________________________________________________________________
void InitializeFluidSimulation(VectorParticles &vectorParticles, int numLayers, int scene)
{
    if (scene == 1)
    {
        InitializeWaterColumn(vectorParticles, numLayers);
    }
    else if (scene == 2)
    {
        InitializeBreakingDamSimple(vectorParticles, numLayers, 1);
    }
    else if (scene == 5)
    {
        // Sketch-imported scene: tools/scene_import/example_scenes/tsunami_realistic_wider.txt
        InitializeFromSceneFile(vectorParticles,
                               "tools/scene_import/example_scenes/tsunami_realistic_wider.txt", numLayers
                            ,1111.f, 0.001f, 1000.f);
    }
    else if (scene == 6)
    {
        // Sketch-imported scene: tools/scene_import/example_scenes/breaking_dam_wide.txt
        InitializeFromSceneFile(vectorParticles,
                               "tools/scene_import/example_scenes/breaking_dam_wide.txt",
                               numLayers);
    }
    else if (scene == 7)
    {
        InitializeFromSceneFile(vectorParticles,
                               "tools/scene_import/example_scenes/downwards.txt", numLayers,
                               1111.f, 0.001f, 10000.f);
    }
    else if (scene == 8)
    {
        InitializeFromSceneFile(vectorParticles,
                               "tools/scene_import/example_scenes/underground.txt", numLayers,
                               1111.f, 0.001f, 10000.f);
    }
    else if (scene == 9)
    {
        InitializeFromSceneFile(vectorParticles,
                               "tools/scene_import/example_scenes/tsunami.txt", numLayers,
                               1111.f, 0.001f, 7500.f);
    }
    else if (scene == 10)
    {
        InitializeBreakingDamWideLow(vectorParticles, numLayers);
    }
}

// ____________________________________________________________________________________________________________________________________
void InitializeFromSceneFile(VectorParticles &vectorParticles, const std::string &path,
                             int numLayers, float restDensity, float timeStep,
                             float gasStiffness, float viscosity)
{
    std::ifstream file(path);
    if (!file.is_open())
    {
        throw std::runtime_error("InitializeFromSceneFile: failed to open '" + path + "'");
    }

    struct RawParticle
    {
        float x, y;
        bool isFluid;
    };

    std::vector<RawParticle> raw;
    float minX = std::numeric_limits<float>::max();
    float minY = std::numeric_limits<float>::max();
    float maxX = std::numeric_limits<float>::lowest();
    float maxY = std::numeric_limits<float>::lowest();

    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty()) continue;
        std::istringstream iss(line);
        float x, y;
        std::string type;
        if (!(iss >> x >> y >> type)) continue;

        bool isFluid = (type == "fluid");
        raw.push_back({x, y, isFluid});

        minX = std::min(minX, x);
        maxX = std::max(maxX, x);
        minY = std::min(minY, y);
        maxY = std::max(maxY, y);
    }

    if (raw.empty())
    {
        throw std::runtime_error("InitializeFromSceneFile: no particles read from '" + path + "'");
    }

    unsigned int numFluid = 0;
    for (const RawParticle &p : raw)
        if (p.isFluid) numFluid++;
    unsigned int numBoundary = static_cast<unsigned int>(raw.size()) - numFluid;

    SetRestDensity(restDensity);
    SetTimeStep(timeStep);
    SetGasStiffnes(gasStiffness);
    SetViscosity(viscosity);
    ComputeSpacing();       // SPACING = sqrt(PARTICLE_VOLUME); every other scene gets this via
                            // SetupScene() - skip it and SPACING stays 0, corrupting every
                            // cell-index computation (division by zero) down the line.
    ComputeParticleMass();  // depends on REST_DENSITY, set just above.

    // This scene isn't a rectangle, so SetupScene's row/col setup doesn't apply. But
    // IndexSort/SpatialHashing still need CONTAINER_SIZE_X/Y and
    // NUM_BOUNDARY_PER_ROW/COLUMN sized sensibly to size their cell-reference table, so
    // derive them from the scene's actual bounding box, in particle-spacing units.
    unsigned int containerSizeX =
        std::max(1u, static_cast<unsigned int>(std::ceil((maxX - minX) / SPACING)));
    unsigned int containerSizeY =
        std::max(1u, static_cast<unsigned int>(std::ceil((maxY - minY) / SPACING)));
    SetContainerSize(containerSizeX, containerSizeY);
    ComputeNumBoundaryParticles(numLayers);

    // Override with the real counts from the file - ComputeNumBoundaryParticles only sized
    // the lookup table above; it doesn't know this scene's actual (non-rectangular) count.
    NUM_PARTICLES_FLUID = numFluid;
    NUM_PARTICLES_FLUID_CURRENT = numFluid;
    NUM_PARTICLES_BOUNDARY = numBoundary;

    // Elsewhere in the codebase (SceneGeometry::computeMinsAndMaxs, in particular) assumes
    // vectorParticles is laid out fluid-first, boundary-second - indices [0, numFluid) are
    // fluid, [numFluid, end) are boundary - regardless of what order the file lists them in.
    // Violating that silently breaks anything that slices by index instead of by type.
    vectorParticles.clear();
    vectorParticles.reserve(raw.size());
    for (const RawParticle &p : raw)
    {
        if (!p.isFluid) continue;
        Particle particle(PARTICLE_MASS, REST_DENSITY, SPACING, ParticleType::Fluid,
                          static_cast<unsigned int>(vectorParticles.size()));
        particle.setPosition(p.x, p.y);
        vectorParticles.push_back(particle);
    }
    for (const RawParticle &p : raw)
    {
        if (p.isFluid) continue;
        Particle particle(PARTICLE_MASS, REST_DENSITY, SPACING, ParticleType::Boundary,
                          static_cast<unsigned int>(vectorParticles.size()));
        particle.setPosition(p.x, p.y);
        vectorParticles.push_back(particle);
    }
}
