/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#include "InputManager.h"

#include <optional>

#include "../visualization/BoundaryFill.h"
#include "FrameCapture.h"
#include "../scene/ParticleSettings.h"
#include "../visualization/ParticleVisualization.h"

void ProcessUserInput(sf::RenderWindow &window, bool &record, bool &showParameters, bool &paused,
                      VectorParticles &vectorParticles, sf::VertexArray &particlesVA,
                      int &numLayers, int &currentScene, IndexSort &indexSort,
                      IndexSort &zIndexSort)
{
    while (const std::optional event = window.pollEvent())
    {
        if (event->is<sf::Event::Closed>())
        {
            window.close();
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Q))
        {
            window.clear();
            vectorParticles.clear();
            particlesVA.clear();
            window.close();
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::P))
        {
            showParameters = !showParameters;
            sf::sleep(sf::milliseconds(150));
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Space))
        {
            paused = !paused;
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num1))
        {
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            currentScene = 1;
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num2))
        {
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            currentScene = 2;
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num3))
        {
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            currentScene = 3;
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num4))
        {
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            currentScene = 4;
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num5))
        {
            // Sketch-imported scene: tools/scene_import/example_scenes/tsunami_realistic_wider.txt
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            currentScene = 5;
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num6))
        {
            // Sketch-imported scene: tools/scene_import/example_scenes/breaking_dam_wide.txt
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            currentScene = 6;
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num7))
        {
            // Sketch-imported scene: tools/scene_import/example_scenes/downwards.txt
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            currentScene = 7;
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num8))
        {
            // Sketch-imported scene: tools/scene_import/example_scenes/underground.txt
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            currentScene = 8;
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num9))
        {
            // Sketch-imported scene: tools/scene_import/example_scenes/tsunami.txt (~100k fluid)
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            currentScene = 9;
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num0))
        {
            // Procedural wide-low breaking dam (InitializeBreakingDamWideLow), ~1M fluid
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            currentScene = 10;
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::R))
        {
            record = !record;
            sf::sleep(sf::milliseconds(150));
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::S))
        {
            FrameCapture::screenShot(window, "screenshots", "s");
        }
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::L))
        {
            numLayers = (numLayers % 3) + 1;  // Cycle between 1, 2, 3 layers.
            sf::sleep(sf::milliseconds(150));
            paused = true;
            vectorParticles.clear();
            indexSort.reset();
            zIndexSort.reset();
            InitializeFluidSimulation(vectorParticles, numLayers, currentScene);
            BoundaryFill::Rebuild(vectorParticles, NUM_PARTICLES_FLUID, SPACING);
            SetSimulationView(window, vectorParticles);
        }
    }
}