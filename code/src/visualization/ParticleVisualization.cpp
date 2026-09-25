/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#include "ParticleVisualization.h"

#include <SFML/Graphics.hpp>
#include <SFML/Graphics/PrimitiveType.hpp>
#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>  // for srand(), rand()
#include <ctime>    // for time()
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

#include "BoundaryFill.h"
#include "../utilities/MathUtils.h"
#include "../simulation/NeighbourSearch.h"
#include "../scene/ParticleSettings.h"
#include "../scene/SceneGeometry.h"
#include "../utilities/SimulationStats.h"

// The vertex array for all particles, initialized once and reused for all frames.
sf::VertexArray particlesVA;

bool particleVAInitialized = false;

std::optional<sf::Shader> particleShader;

// ____________________________________________________________________________________________________________________________________
void UpdateParticleColorLerp(Particle &particle)
{
    float v = particle.computeVelocityMagnitude();
    const float v0 = 2.f;   // below -> blue
    const float v1 = 15.f;  // above -> white

    sf::Color color;
    sf::Color color0 = blue;
    sf::Color color1 = white;
    if (v <= v0)
    {
        color = color0;
    }
    else if (v >= v1)
    {
        color = color1;
    }
    else
    {
        float t = (v - v0) / (v1 - v0);
        color = Math::lerpColor(color0, color1, t);
    }
    particle.updateColor(color);
}

// ___________________________________________________________________________________________________________________________________
void UpdateColorLerp(VectorParticles &vectorParticles)
{
    // Each particle's update is independent, so this parallelizes safely.
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(vectorParticles.size()); ++i)
    {
        Particle &particle = vectorParticles[i];
        if (particle.type == ParticleType::Fluid) UpdateParticleColorLerp(particle);
    }
}

// ____________________________________________________________________________________________________________________________________
// Thesis-figure helper: color every particle by its cell index, normalized against the range
// of cell indices actually occupied by particles (not the grid's full theoretical range, which
// would crush the whole gradient into a sliver since the scene never fills the whole grid).
// Uses its own grid resolution (VIZ_CELL_SIZE_FACTOR below) for a smoother gradient than the
// real neighbor-search cell size would give, and computes cell indices directly via
// UniformGrid::MapPosToCellIndices rather than IndexSort::construct, since IndexSort writes
// into fields the real simulation's own neighbor search depends on.
void UpdateColorByCellIndex(VectorParticles &vectorParticles, IndexingScheme scheme,
                            sf::Color colorLow, sf::Color colorHigh)
{
    if (vectorParticles.empty()) return;

    // Fraction of the real neighbor-search cell size to use for this coloring's own grid -
    // smaller means more, smaller color blocks. 0.5 doubles the grid resolution in each axis
    // (4x as many cells) compared to the real simulation grid.
    constexpr float VIZ_CELL_SIZE_FACTOR = 0.5f;
    float cellSize = VIZ_CELL_SIZE_FACTOR * 2.f * vectorParticles[0].spacing;
    int vizRows = static_cast<int>(std::round(NUM_BOUNDARY_PER_ROW / VIZ_CELL_SIZE_FACTOR));

    UniformGrid grid;
    std::vector<int> cellIndexPerParticle(vectorParticles.size());
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(vectorParticles.size()); ++i)
    {
        CellIndices cellIndices = grid.MapPosToCellIndices(vectorParticles[i].position, cellSize);
        // NOT Indexing::computeCellIndex/computeLinearIndex for LinearCurve: that hardcodes
        // K = NUM_BOUNDARY_PER_ROW, the real (coarse) simulation grid's width, not this finer
        // viz grid's, so fine-grid (k,l) pairs would collide. Z-curve's Morton code has no such
        // dependency, so it's unaffected and can still go through computeCellIndex.
        cellIndexPerParticle[i] =
            scheme == IndexingScheme::ZCurve
                ? static_cast<int>(Indexing::computeZIndex(cellIndices))
                : cellIndices.k + cellIndices.l * vizRows;
    }

    int rangeMin = std::numeric_limits<int>::max();
    int rangeMax = std::numeric_limits<int>::min();
    for (int cellIndex : cellIndexPerParticle)
    {
        rangeMin = std::min(rangeMin, cellIndex);
        rangeMax = std::max(rangeMax, cellIndex);
    }
    if (rangeMax <= rangeMin) return;

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(vectorParticles.size()); ++i)
    {
        float t = static_cast<float>(cellIndexPerParticle[i] - rangeMin) /
                 static_cast<float>(rangeMax - rangeMin);
        vectorParticles[i].updateColor(Math::lerpColor(colorLow, colorHigh, t));
    }
}

// Recomputed every frame so a particle's color always tracks whatever cell it's currently in -
// as it crosses a cell boundary it releases the old cell's color and picks up the new one.
void UpdateColorByLinearCellIndex(VectorParticles &vectorParticles)
{
    UpdateColorByCellIndex(vectorParticles, IndexingScheme::LinearCurve, sf::Color(130, 30, 190),
                          sf::Color(255, 220, 0));
}

void UpdateColorByZCurveCellIndex(VectorParticles &vectorParticles)
{
    UpdateColorByCellIndex(vectorParticles, IndexingScheme::ZCurve, sf::Color(130, 30, 190),
                          sf::Color(255, 220, 0));
}

// ____________________________________________________________________________________________________________________________________
sf::Shader CreateParticleShader()
{
    // gl_PointSize is set here (in pixels); the caller must enable GL_PROGRAM_POINT_SIZE on
    // the active context once, or this has no effect (point size silently stays at its
    // fixed default). Uses the legacy gl_Vertex/gl_ModelViewProjectionMatrix/gl_Color
    // built-ins, same compatibility profile SFML's own fixed-function drawing relies on -
    // no manual attribute/uniform-matrix wiring needed.
    // #version must be the shader's literal first line (no leading whitespace) per the GLSL
    // spec - this GL implementation doesn't expose gl_PointCoord without it.
    static const char *vertexSrc =
        "#version 120\n"
        R"glsl(
        uniform float pointSize;

        void main()
        {
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
            gl_FrontColor = gl_Color;
            gl_PointSize = pointSize;
        }
    )glsl";

    // Circular falloff computed analytically per-pixel from gl_PointCoord instead of
    // sampling a pre-baked texture - no texture minification, so no moire from that source
    // when many particles are drawn small (e.g. 100k particles zoomed out).
    static const char *fragmentSrc =
        "#version 120\n"
        R"glsl(
        uniform float softness;  // fraction of the radius the edge fades over, 0-1

        void main()
        {
            vec2 coord = gl_PointCoord * 2.0 - vec2(1.0);
            float dist = length(coord);
            float alpha = 1.0 - smoothstep(1.0 - softness, 1.0, dist);
            if (alpha <= 0.0) discard;
            gl_FragColor = vec4(gl_Color.rgb, gl_Color.a * alpha);
        }
    )glsl";

    sf::Shader shader;
    if (!shader.loadFromMemory(vertexSrc, fragmentSrc))
    {
        std::cerr << "Failed to compile particle point-sprite shader.\n";
    }
    return shader;
}

// ____________________________________________________________________________________________________________________________________
void DrawParticlesPoints(sf::RenderTarget &target, VectorParticles &vectorParticles,
                         sf::VertexArray &particlesVA, std::size_t numParticles,
                         std::size_t startIndex)
{
    if (particlesVA.getPrimitiveType() != sf::PrimitiveType::Points)
        particlesVA.setPrimitiveType(sf::PrimitiveType::Points);
    if (particlesVA.getVertexCount() != numParticles) particlesVA.resize(numParticles);

    // One vertex per particle instead of six - less CPU work building the array, less data
    // uploaded to the GPU per frame. Each write lands in a disjoint slot, safe to parallelize.
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(numParticles); ++i)
    {
        const Particle &particle = vectorParticles[startIndex + i];
        particlesVA[i].position = particle.position;
        particlesVA[i].color = particle.color;
    }

    // All particles share one radius (SPACING / 2, see Particle.h), so point size in pixels
    // is uniform across the draw call - compute it once rather than per vertex.
    sf::View view = target.getView();
    float worldToPixel = view.getSize().x > 0.f
                             ? static_cast<float>(target.getSize().x) / view.getSize().x
                             : 1.f;
    float pointSizePixels =
        numParticles > 0 ? 2.f * vectorParticles[startIndex].radius * worldToPixel : 0.f;

    particleShader->setUniform("pointSize", pointSizePixels);
    particleShader->setUniform("softness", 0.35f);

    sf::RenderStates states;
    states.shader = &particleShader.value();
    target.draw(particlesVA, states);
}

// ____________________________________________________________________________________________________________________________________
void DrawParticles(sf::RenderTarget &target, VectorParticles &vectorParticles,
                   sf::VertexArray &particlesVA)
{
    const std::size_t numParticles = vectorParticles.size();
    UpdateColorLerp(vectorParticles);
    DrawParticlesPoints(target, vectorParticles, particlesVA, numParticles, 0);
}

// ____________________________________________________________________________________________________________________________________
void DrawParameters(sf::RenderTarget &target, const sf::Font &font, sf::View simView,
                    int windowWidth, int windowHeight, bool paused, bool recording, bool logging,
                    float densityError, float cflNumber, int frameCounter)
{
    const char *pausedChar = paused ? "[PAUSED]" : "";
    const char *recordingChar = recording ? "[RECORDING]" : "";
    const char *loggingChar = logging ? "[LOGGING]" : "";
    std::ostringstream oss;

    oss << "n = " << NUM_PARTICLES_FLUID_CURRENT << "\n"
        << "dt = " << TIMESTEP << "\n"
        << "k = " << GAS_STIFFNESS << "\n"
        << "visc = " << VISCOSITY << "\n"
        << "h = " << SPACING << "\n"
        << "error = " << std::fixed << std::setprecision(3) << densityError << "\n"
        << "CFL = " << std::fixed << std::setprecision(3) << cflNumber << "\n"
        << "frame = " << frameCounter << "\n"
        << "gamma1 = " << DENSITY_BOUNDARY_FACTOR << "\n"
        << "gamma2 = " << PRESSURE_BOUNDARY_FACTOR << "\n"
        << pausedChar << "\n"
        << recordingChar << "\n"
        << loggingChar << "\n";

    sf::Text text(font, oss.str(), 10);
    text.setFillColor(sf::Color::White);
    text.setPosition(V2f(10, 5));
    target.setView(target.getDefaultView());  // draw in screen coordinates
    target.draw(text);
    target.setView(simView);  // restore simulation view for next frame's drawing
}

// ____________________________________________________________________________________________________________________________________
void DrawScene(sf::RenderTarget &target, VectorParticles &vectorParticles,
               sf::VertexArray &particlesVA, const sf::Font &font, sf::View simView,
               int windowWidth, int windowHeight, bool paused, bool recording, bool logging,
               int frameCounter, bool showParameters)
{
    BoundaryFill::Draw(target);  // drawn first so particles render on top of it
    DrawParticles(target, vectorParticles, particlesVA);
    if (showParameters)
    {
        float densityError = SimulationStats::computeDensityError(
            SimulationStats::computeAverageDensity(vectorParticles));
        float cflNumber =
            SimulationStats::computeCflNumber(SimulationStats::findMaxSpeed(vectorParticles));
        DrawParameters(target, font, simView, windowWidth, windowHeight, paused, recording, logging,
                       densityError, cflNumber, frameCounter);
    }
}

// ____________________________________________________________________________________________________________________________________
void InitializeVisualization(int numParticlesFluid, int numParticlesBoundary, sf::Font &font,
                            const std::string &fontPath)
{
    particleShader = CreateParticleShader();
    particlesVA.setPrimitiveType(sf::PrimitiveType::Points);
    // pre-allocate max or a sensible upper bound - one vertex per particle now, not six
    particlesVA.resize(numParticlesFluid + numParticlesBoundary);

    const std::array<std::string, 4> fontPaths = {
        fontPath,
        "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/System/Library/Fonts/Times.ttc",
        "fonts/arial.ttf",
    };

    bool loaded = false;
    for (const std::string &candidatePath : fontPaths)
    {
        if (font.openFromFile(candidatePath))
        {
            loaded = true;
            break;
        }
    }

    if (!loaded)
    {
        std::cerr << "Failed to load visualization font.\n";
    }
}

// ____________________________________________________________________________________________________________________________________
void ShutdownParticleVisualization()
{
    // Resetting this (rather than letting the C++ runtime destroy it as an ordinary static
    // after main() returns) releases its GL resources while the window's context is still
    // current. See the declaration above / ParticleVisualization.h for why this matters.
    particleShader.reset();
}

// ____________________________________________________________________________________________________________________________________
void SetSimulationView(sf::RenderWindow &window, VectorParticles &vectorParticles)
{
    V2f worldSize =
        SceneGeometry::computeWorldSize(vectorParticles, NUM_PARTICLES_FLUID_CURRENT, true);
    float worldSizeMax = worldSize.x > worldSize.y ? worldSize.x : worldSize.y;
    sf::View simulationView(sf::FloatRect(V2f(0, 0), V2f(worldSizeMax, worldSizeMax)));
    V2f center =
        SceneGeometry::computeCenterPosition(vectorParticles, NUM_PARTICLES_FLUID_CURRENT, true);
    simulationView.setCenter(center);
    simulationView.zoom(1.2f);  // v > 0 zoomout, v < 0 zoom in
    window.setView(simulationView);
}

// ____________________________________________________________________________________________________________________________________
sf::View ComputeFittedSimulationView(VectorParticles &vectorParticles, float zoomMargin)
{
    V2f worldSize =
        SceneGeometry::computeWorldSize(vectorParticles, NUM_PARTICLES_FLUID_CURRENT, true);
    V2f center =
        SceneGeometry::computeCenterPosition(vectorParticles, NUM_PARTICLES_FLUID_CURRENT, true);
    sf::View fittedView(sf::FloatRect(V2f(0, 0), worldSize));
    fittedView.setCenter(center);
    fittedView.zoom(zoomMargin);
    return fittedView;
}

// ____________________________________________________________________________________________________________________________________
sf::Vector2u ComputeFittedRenderSize(VectorParticles &vectorParticles, unsigned int maxDimension)
{
    V2f worldSize =
        SceneGeometry::computeWorldSize(vectorParticles, NUM_PARTICLES_FLUID_CURRENT, true);
    if (worldSize.x <= 0.f || worldSize.y <= 0.f) return {maxDimension, maxDimension};

    float aspect = worldSize.x / worldSize.y;  // width / height
    unsigned int width, height;
    if (aspect >= 1.f)
    {
        width = maxDimension;
        height = static_cast<unsigned int>(std::round(maxDimension / aspect));
    }
    else
    {
        height = maxDimension;
        width = static_cast<unsigned int>(std::round(maxDimension * aspect));
    }
    return {std::max(width, 1u), std::max(height, 1u)};
}
