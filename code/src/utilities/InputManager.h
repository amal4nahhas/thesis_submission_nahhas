/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFML/Window.hpp>

#include "../simulation/NeighbourSearch.h"
#include "../simulation/Particle.h"

/// @brief Process user input events.
/// @param window
/// @param record
/// @param showParameters
/// @param paused
/// @param indexSort reset on every scene switch, since it otherwise keeps stale handles
///        referencing the previous scene's (now-replaced) vectorParticles.
/// @param zIndexSort same as indexSort.
void ProcessUserInput(sf::RenderWindow &window, bool &record, bool &showParameters, bool &paused,
                      VectorParticles &vectorParticles, sf::VertexArray &particlesVA,
                      int &numLayers, int &currentScene, IndexSort &indexSort,
                      IndexSort &zIndexSort);