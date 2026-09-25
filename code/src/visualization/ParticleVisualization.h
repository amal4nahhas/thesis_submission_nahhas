/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <optional>

#include "../simulation/NeighbourSearch.h"
#include "../simulation/Particle.h"

extern bool particleVAInitialized;
extern sf::VertexArray particlesVA;

// Point-sprite shader used by DrawParticlesPoints: renders each particle as one GL point
// instead of a 6-vertex textured quad, with the circular falloff computed analytically in
// the fragment shader (no texture sampling, so no minification/moire from that source).
// std::optional (not a plain sf::Shader) so ShutdownParticleVisualization() can release its
// GL resources explicitly before the window's GL context is destroyed - see that function.
extern std::optional<sf::Shader> particleShader;

/// @brief Build the point-sprite shader used by DrawParticlesPoints.
sf::Shader CreateParticleShader();

/// @brief Draw particles as GL points with `particleShader`, one vertex per particle
/// instead of six. Requires GL_PROGRAM_POINT_SIZE enabled on the active context (done once
/// in SPH.cpp after window creation) so the vertex shader's gl_PointSize takes effect.
/// Takes a RenderTarget rather than a RenderWindow specifically so the same drawing code
/// also works against an off-screen RenderTexture (see FrameCapture::screenShotFitted).
/// @param target
/// @param vectorParticles
/// @param particlesVA
/// @param numParticles
/// @param startIndex draws vectorParticles[startIndex, startIndex+numParticles)
void DrawParticlesPoints(sf::RenderTarget &target, VectorParticles &vectorParticles,
                         sf::VertexArray &particlesVA, std::size_t numParticles,
                         std::size_t startIndex = 0);

/// @brief Update the color of a particle based on its speed using linear interpolation.
/// @param particle
void UpdateParticleColorLerp(Particle &particle);

/// @brief Thesis-figure helper: color every particle by its cell index under linear indexing,
/// purple (lowest) to yellow (highest), normalized against the uniform grid's own theoretical
/// index range. Swap into DrawParticles in place of UpdateColorLerp to use - see the doc
/// comment above its definition for why this isn't a runtime-toggleable mode.
/// @param vectorParticles
void UpdateColorByLinearCellIndex(VectorParticles &vectorParticles);

/// @brief Same as UpdateColorByLinearCellIndex, but under Z-curve (Morton) indexing instead.
/// @param vectorParticles
void UpdateColorByZCurveCellIndex(VectorParticles &vectorParticles);

/// @brief Draw the particles at their given locations, one GL point per particle via
/// `particleShader`.
/// @param target
/// @param vectorParticles
/// @param particlesVA
void DrawParticles(sf::RenderTarget &target, VectorParticles &vectorParticles,
                   sf::VertexArray &particlesVA);

/// @brief Release the GL resources held by particleShader. It has translation-unit-static
/// storage duration, so left to the C++ runtime it'd be destroyed AFTER main() returns - i.e.
/// after the window's GL context is already gone, which crashes on exit. Call this from
/// main() while the window is still alive, after the main loop ends.
void ShutdownParticleVisualization();

/// @brief Draw a text box with simulation parameters in the top-left corner.
/// @param target
/// @param font
/// @param simView
/// @param windowWidth
/// @param windowHeight
void DrawParameters(sf::RenderTarget &target, const sf::Font &font, sf::View simView,
                    int windowWidth, int windowHeight, bool paused, bool recording, bool logging,
                    float densityError, float cflNumber, int frameCounter);

/// @brief Draw the entire scene: particles and parameters.
/// @param target
/// @param vectorParticles
/// @param particlesVA
/// @param font
/// @param simView
/// @param windowWidth
/// @param windowHeight
/// @param paused
/// @param recording
/// @param logging
/// @param showParameters
void DrawScene(sf::RenderTarget &target, VectorParticles &vectorParticles,
               sf::VertexArray &particlesVA, const sf::Font &font, sf::View simView,
               int windowWidth, int windowHeight, bool paused, bool recording, bool logging,
               int frameCounter, bool showParameters = true);

/// @brief Initialize the visualization components such as the particle shader and vertex array.
/// @param numParticlesFluid
/// @param numParticlesBoundary
/// @param font
/// @param fontPath
void InitializeVisualization(int numParticlesFluid, int numParticlesBoundary, sf::Font &font,
                            const std::string &fontPath);

/// @brief Set the simulation view based on particle positions.
/// @param window
/// @param vectorParticles
void SetSimulationView(sf::RenderWindow &window, VectorParticles &vectorParticles);

/// @brief Compute a simulation view sized to the particle bounding box's own aspect ratio,
/// instead of the fixed square view SetSimulationView uses. Meant for off-screen capture
/// (see FrameCapture::screenShotFitted), where the render target's own aspect ratio is chosen
/// to match this view, so the captured frame isn't padded with empty space or stretched.
/// @param vectorParticles
/// @param zoomMargin fraction of extra space left around the bounding box (matches
/// SetSimulationView's own 1.2x margin by default)
sf::View ComputeFittedSimulationView(VectorParticles &vectorParticles, float zoomMargin = 1.2f);

/// @brief Compute an off-screen render size whose aspect ratio matches the particle bounding
/// box, so a ComputeFittedSimulationView capture isn't stretched. The longer side is set to
/// maxDimension pixels.
/// @param vectorParticles
/// @param maxDimension
sf::Vector2u ComputeFittedRenderSize(VectorParticles &vectorParticles,
                                     unsigned int maxDimension = 1600);
