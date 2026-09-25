/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include <filesystem>
#include <fstream>
#include <string>

#include "../simulation/Particle.h"
#include "../scene/ParticleSettings.h"
#include "../visualization/ParticleVisualization.h"

class SimulationStats
{
   public:
    static bool logStats(clock_t begin, clock_t &lastRecord, float loggingTime,
                         VectorParticles &vectorParticles, std::ofstream &Time,
                         std::ofstream &Error, std::ofstream &CFL, std::ofstream &Param,
                         bool loggingEnabled);

    static void logTimeForNeighbourSearch(clock_t begin, int currentSimStep, int maxNumSteps,
                                          std::ofstream &time);

    static void updateColorTest(VectorParticles &vectorParticles);

    static float computeAverageVelocity(VectorParticles &vectorParticles);

    static float computeBoundaryNeighbourRatio(const Particle &particle,
                                               const VectorParticles &vectorParticles);

    static float findMaxSpeed(VectorParticles &vectorParticles);

    static float computeAverageDensity(VectorParticles &vectorParticles);

    static float computeDensityError(float averageDensity);

    static float computeCflNumber(float maxSpeed);

    static float computeCflThreshold(float cflNumber, float maxSpeed);

    static bool isEquilibrium(VectorParticles &vectorParticles, int frameCounter,
                              float threshold = 0.2f);

    /// @brief Static-scene neighbour-search benchmark: for a sweep of fluid particle counts,
    /// times a single construct()+query() for every neighbour-search algorithm (all IndexSort
    /// indexing/sorting combinations, plus SpatialHashing and CompactHashing), saves a
    /// screenshot of each scene, and writes results to performance/static/4thTry/.
    /// `window`/`particlesVA`/`font` must already be set up by the caller as for normal
    /// interactive drawing.
    static void staticSceneAnalysis(VectorParticles &vectorParticles, sf::RenderWindow &window,
                                    sf::VertexArray &particlesVA, sf::Font &font);

    /// @brief Dynamic-scene counterpart to staticSceneAnalysis: for each particle count and
    /// the two radix-sorted index sort combos, builds a wide-low breaking dam scene at that
    /// size and runs a real simulation loop, timing neighbour-search construct/query
    /// separately from the rest of the physics update on every step. Headless - only the
    /// initial scene screenshot uses the window. Writes performance/dynamic/scaling/results.csv.
    static void dynamicSceneAnalysis(VectorParticles &vectorParticles, sf::RenderWindow &window,
                                     sf::VertexArray &particlesVA, sf::Font &font);

    /// @brief For a sweep of fluid particle counts, builds a static water column, runs a
    /// single construct()+query() with index sort (Z-curve indexing, radix sort) to populate
    /// every fluid particle's neighbours, then compresses those lists with
    /// NeighbourCompression and reports the achieved compression ratio against the raw
    /// 4-bytes-per-neighbour representation. Headless. Writes
    /// performance/static/compression/results.csv.
    static void neighbourCompressionAnalysis(VectorParticles &vectorParticles);

    /// @brief Memory-only comparison of SpatialHashing vs CompactHashing at the same particle
    /// counts, using each structure's own analytic estimatedMemoryBytes() (see
    /// NeighbourSearch.h) rather than OS-reported RSS, which becomes unreliable once particle
    /// counts exceed physical RAM. Headless. Writes
    /// performance/static/4thTry/memoryComparison/results.csv.
    static void memoryComparisonAnalysis(VectorParticles &vectorParticles);

   private:
    struct StaticAlgResult
    {
        std::string name;
        std::string indexingScheme;
        std::string sortingMethod;
        double constructSec;
        double querySec;
        bool skipped = false;  // no measurement taken - write '-' instead of a time
    };

    static std::filesystem::path findProjectRoot();
    static std::ofstream openOutputFile(const std::filesystem::path &relativePath);
    static void writeStaticPerAlgorithmFile(const std::string &runDir, const StaticAlgResult &r,
                                            long particleCount);
    static void appendStaticCsvRow(std::ofstream &csv, long particleCount, unsigned int nFluid,
                                   unsigned int nBoundary, const StaticAlgResult &r,
                                   double maxRssMb);
};
