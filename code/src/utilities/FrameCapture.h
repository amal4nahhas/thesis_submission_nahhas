/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include <SFML/Graphics.hpp>

#include <filesystem>

#include "../simulation/Particle.h"

class FrameCapture
{
   public:
    static void screenShot(sf::RenderWindow &window, std::string dirName, std::string fileName);

    /// @brief Render the scene into an off-screen texture sized and viewed to exactly match
    /// the particle bounding box's own aspect ratio (see ComputeFittedSimulationView /
    /// ComputeFittedRenderSize in ParticleVisualization.h), and save it as a PNG. Independent
    /// of any on-screen window's size, so it never disturbs what the user is looking at (and
    /// unlike screenShot(), doesn't depend on the compositor having processed a frame first).
    /// @param vectorParticles
    /// @param particlesVA reused across calls, same as the live rendering path
    /// @param font
    /// @param dirName
    /// @param fileName written as \<dirName\>/\<fileName\>.png
    /// @param maxDimension the bounding box's longer side, in pixels
    /// @param backgroundColor cleared onto the off-screen texture before drawing - pass
    /// sf::Color::Transparent for a PNG with a real alpha channel instead of a solid fill,
    /// since this render target (unlike an on-screen window) supports it natively.
    static void screenShotFitted(VectorParticles &vectorParticles, sf::VertexArray &particlesVA,
                                 const sf::Font &font, const std::string &dirName,
                                 const std::string &fileName, unsigned int maxDimension = 1600,
                                 sf::Color backgroundColor = sf::Color(235, 235, 235));

    static void record(sf::RenderWindow &window, std::string dirName, int &frameCounter);

    static void recordSimulation(sf::RenderWindow &window, int &frameCounter);

    static void exportGif(const std::string &frameDir, const std::string &outputPath, int delay);

    static void clearFrameDirectory(const std::string &dirName);

   private:
    static constexpr int SCREENSHOT_UPSCALE = 2;

    /// @brief Apply a high-quality Lanczos upscale when ImageMagick is available.
    static bool TryEnhanceScreenshot(const std::string &path);

    static int ReadCounterFromFile(const std::filesystem::path &counterPath);
    static void WriteCounterToFile(const std::filesystem::path &counterPath, int counter);
};
