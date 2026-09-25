/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/

#define GL_SILENCE_DEPRECATION  // legacy OpenGL (used below via glEnable) is deprecated on
                                // macOS in favor of Metal, but still functions
#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <SFML/Window.hpp>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>

#include "../visualization/BoundaryFill.h"
#include "../utilities/InputManager.h"
#include "../utilities/MathUtils.h"
#include "Particle.h"
#include "../scene/ParticleSettings.h"
#include "ParticleSimulation.h"
#include "../visualization/ParticleVisualization.h"
#include "../scene/SceneGeometry.h"
#include "../utilities/Timer.h"
#include "../utilities/FrameCapture.h"

constexpr unsigned int WINDOW_WIDTH = 900;
constexpr unsigned int WINDOW_HEIGHT = 900;

int NUM_BOUNDARY_LAYERS = 1;
static bool SHOW_PARAMETERS = 1;
static bool PAUSED = 0;
// Frame capture and GIF export
static bool RECORDING_ENABLED = 0;
static bool LOGGING_ENABLED = 0;
static float LOGGING_TIME = 120.f * 4000.f;
int FRAME_COUNTER = 0;
static int MAX_NUM_FRAMES = 1000;
static int CURRENT_SCENE = 5;
static int NEIGHBOUR_SEARCH = 2;
static bool TESTING = 0;
static bool framesDirCleared = false;
int capturedFrames = 1;

int main()
{
    // Initialize the simulation scene; press 0-9 to switch scenes live.
    InitializeFluidSimulation(vectorParticles, NUM_BOUNDARY_LAYERS, CURRENT_SCENE);
    BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);

    // Create the main window first - font/particleTexture below hold GL resources, and local
    // variables are destroyed in reverse declaration order, so declaring them after window
    // ensures they're released while the window's GL context still exists.
    sf::ContextSettings contextSettings;
    contextSettings.antiAliasingLevel = 4;
    sf::RenderWindow window(sf::VideoMode({WINDOW_WIDTH, WINDOW_HEIGHT}), "SPH Simulation",
                            sf::Style::Resize, sf::State::Windowed, contextSettings);

    // Needed for the particle point-sprite shader's gl_PointSize to take effect - state is
    // per-context, so this must run after the window (and its GL context) exists.
    // 0x8642 = GL_PROGRAM_POINT_SIZE / GL_VERTEX_PROGRAM_POINT_SIZE - macOS's legacy
    // <OpenGL/gl.h> doesn't define either name, so the raw value is used directly.
    glEnable(0x8642);

    sf::Font font;
    InitializeVisualization(NUM_PARTICLES_FLUID, NUM_PARTICLES_BOUNDARY, font,
                            "/System/Library/Fonts/Supplemental/Times New Roman.ttf");

    SetSimulationView(window, vectorParticles);

    auto minsAndMAx =
        SceneGeometry::computeMinsAndMaxs(vectorParticles, NUM_PARTICLES_FLUID_CURRENT);

    ParticleSimulation particleSimulation;
    // Instantiate the neighbourhood search algorithms. Only indexSort is actually queried per
    // frame below - spatialHashing/compactHashing are still constructed since simulationStep's
    // signature takes them by reference, but sized trivially, not at full particle count.
    IndexSort indexSort(IndexingScheme::LinearCurve, SortingMethod::InsertionSort);
    IndexSort zIndexSort(IndexingScheme::ZCurve, SortingMethod::RadixSort);
    SpatialHashing spatialHashing(1);
    CompactHashing compactHashing;

    // Start the main loop
    while (window.isOpen())
    {
        ProcessUserInput(window, RECORDING_ENABLED, SHOW_PARAMETERS, PAUSED, vectorParticles,
                         particlesVA, NUM_BOUNDARY_LAYERS, CURRENT_SCENE, indexSort, zIndexSort);

        window.clear(transparent);

        if (!PAUSED && !TESTING)
        {
            indexSort.construct(vectorParticles, IndexingScheme::ZCurve, SortingMethod::RadixSort);
            indexSort.query(vectorParticles, IndexingScheme::ZCurve);
            particleSimulation.simulationStep(
                vectorParticles, compactHashing, spatialHashing, indexSort, zIndexSort,
                SceneGeometry::computeMinsAndMaxs(vectorParticles, NUM_PARTICLES_FLUID_CURRENT),
                NUM_BOUNDARY_LAYERS, NEIGHBOUR_SEARCH, /*swingGravity=*/CURRENT_SCENE == 5);

            // Rotate the camera to match the current gravity tilt (see ComputeSwingingGravity) -
            // the boundary/fluid geometry itself never moves, but rotating the view around it
            // makes the tank look like it's physically tilting while gravity stays visually
            // "down". Only touches rotation - center/size/zoom stay whatever SetSimulationView
            // last set.
            if (CURRENT_SCENE == 5)
            {
                sf::View simView = window.getView();
                simView.setRotation(sf::degrees(CURRENT_GRAVITY_SWING_ANGLE_DEG));
                window.setView(simView);
            }
        }

        DrawScene(window, vectorParticles, particlesVA, font, window.getView(), WINDOW_WIDTH,
                  WINDOW_HEIGHT, PAUSED, RECORDING_ENABLED, LOGGING_ENABLED, FRAME_COUNTER,
                  SHOW_PARAMETERS);

        window.display();

        if (!PAUSED && RECORDING_ENABLED && FRAME_COUNTER % 10 == 0)
        {
            FrameCapture::record(window, "recordings/frames/tsunami", FRAME_COUNTER);
        }

        FRAME_COUNTER++;
    }

    // Release the particle shaders and render textures explicitly while
    // `window` (and its GL context) is still alive.
    ShutdownParticleVisualization();
}