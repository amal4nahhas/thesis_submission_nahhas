/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/

#include "FrameCapture.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "../visualization/ParticleVisualization.h"

// ___________________________________________________________________________________________________________________________________
bool FrameCapture::TryEnhanceScreenshot(const std::string &path)
{
    const std::string checkMagick = "command -v magick >/dev/null 2>&1";
    if (std::system(checkMagick.c_str()) != 0) return false;

    std::ostringstream cmd;
    cmd << "magick \"" << path << "\" "
        << "-filter Lanczos -resize " << SCREENSHOT_UPSCALE * 100
        << "% -unsharp 0x0.75+0.75+0.008 \"" << path << "\"";

    return std::system(cmd.str().c_str()) == 0;
}

// ___________________________________________________________________________________________________________________________________
int FrameCapture::ReadCounterFromFile(const std::filesystem::path &counterPath)
{
    std::ifstream in(counterPath);
    if (!in.is_open()) return 1;

    int counter = 1;
    in >> counter;
    if (!in.good() && !in.eof()) return 1;

    return std::max(1, counter);
}

// ___________________________________________________________________________________________________________________________________
void FrameCapture::WriteCounterToFile(const std::filesystem::path &counterPath, int counter)
{
    std::ofstream out(counterPath, std::ios::trunc);
    if (!out.is_open()) return;

    out << counter;
}

// ___________________________________________________________________________________________________________________________________
void FrameCapture::screenShot(sf::RenderWindow &window, std::string dirName, std::string fileName)
{
    std::filesystem::create_directories(dirName);

    const std::filesystem::path counterPath = std::filesystem::path(dirName) / "counter.txt";
    const int counter = ReadCounterFromFile(counterPath);

    sf::Texture tex(window.getSize());
    tex.update(window);

    sf::Image screenshot = tex.copyToImage();

    std::ostringstream name;
    name << fileName << "_" << std::setfill('0') << std::setw(6) << counter << ".png";
    std::string path = (std::filesystem::path(dirName) / name.str()).string();

    if (!screenshot.saveToFile(path)) return;

    WriteCounterToFile(counterPath, counter + 1);

    (void)TryEnhanceScreenshot(path);
}

// ___________________________________________________________________________________________________________________________________
void FrameCapture::screenShotFitted(VectorParticles &vectorParticles, sf::VertexArray &particlesVA,
                                    const sf::Font &font, const std::string &dirName,
                                    const std::string &fileName, unsigned int maxDimension,
                                    sf::Color backgroundColor)
{
    std::filesystem::create_directories(dirName);

    sf::Vector2u size = ComputeFittedRenderSize(vectorParticles, maxDimension);
    sf::View fittedView = ComputeFittedSimulationView(vectorParticles);

    sf::RenderTexture renderTexture(size);
    renderTexture.setView(fittedView);
    renderTexture.clear(backgroundColor);
    DrawScene(renderTexture, vectorParticles, particlesVA, font, fittedView,
             static_cast<int>(size.x), static_cast<int>(size.y), /*paused=*/true,
             /*recording=*/false, /*logging=*/false, /*frameCounter=*/0, /*showParameters=*/false);
    renderTexture.display();

    sf::Image screenshot = renderTexture.getTexture().copyToImage();
    std::string path = (std::filesystem::path(dirName) / (fileName + ".png")).string();
    if (!screenshot.saveToFile(path)) return;

    (void)TryEnhanceScreenshot(path);
}

// ___________________________________________________________________________________________________________________________________
void FrameCapture::record(sf::RenderWindow &window, std::string dirName, int &frameCounter)
{
    std::filesystem::create_directories(dirName);
    sf::Texture tex(sf::Vector2u(window.getSize().x, window.getSize().y));
    tex.update(window);
    sf::Image screenshot = tex.copyToImage();

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(6) << frameCounter << ".png";
    std::string path = dirName + "/" + oss.str();

    (void)screenshot.saveToFile(path);
    ++frameCounter;
}

// ___________________________________________________________________________________________________________________________________
void FrameCapture::recordSimulation(sf::RenderWindow &window, int &frameCounter)
{
    std::string dirName = "recordings/scene/frames";
    std::filesystem::create_directories(dirName);
    sf::Texture tex(sf::Vector2u(window.getSize().x, window.getSize().y));
    tex.update(window);
    sf::Image screenshot = tex.copyToImage();

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(6) << frameCounter << ".png";
    std::string path = dirName + "/" + oss.str();

    (void)screenshot.saveToFile(path);
}

// ___________________________________________________________________________________________________________________________________
void FrameCapture::exportGif(const std::string &frameDir, const std::string &outputPath, int delay)
{
    std::string cmd = "magick -delay " + std::to_string(delay) + " -loop 0 -dispose Background " +
                      frameDir + "/*.png " + outputPath;
    std::cout << "Exporting GIF: " << cmd << "\n";
    int result = system(cmd.c_str());
    if (result == 0)
    {
        std::cout << "GIF exported successfully to " << outputPath << "\n";
        std::string notifCmd =
            "osascript -e 'display notification \"GIF export completed\" with title \"SPH "
            "Simulation\"'";
        system(notifCmd.c_str());
        std::cout << "\a";
    }
    else
    {
        std::cerr << "Failed to export GIF. Make sure ImageMagick is installed.\n";
    }
}

// ___________________________________________________________________________________________________________________________________
void FrameCapture::clearFrameDirectory(const std::string &dirName)
{
    for (const auto &entry : std::filesystem::directory_iterator(dirName))
    {
        if (entry.path().extension() == ".png")
        {
            std::filesystem::remove(entry.path());
            std::cout << "Deleted: " << entry.path().filename() << "\n";
        }
    }
}
