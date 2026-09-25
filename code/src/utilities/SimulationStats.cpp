/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/

#include "SimulationStats.h"

#include <sys/resource.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>

#include "../utilities/MathUtils.h"
#include "../simulation/NeighbourCompression.h"
#include "../simulation/NeighbourSearch.h"
#include "../scene/ParticleSettings.h"
#include "../simulation/ParticleSimulation.h"
#include "../scene/SceneGeometry.h"
#include "../utilities/Timer.h"

extern int NUM_BOUNDARY_LAYERS;

std::filesystem::path SimulationStats::findProjectRoot()
{
    std::error_code errorCode;
    auto currentPath = std::filesystem::current_path(errorCode);
    if (errorCode)
    {
        throw std::runtime_error("Failed to read current working directory.");
    }

    for (auto path = currentPath; !path.empty(); path = path.parent_path())
    {
        if (std::filesystem::exists(path / "performance"))
        {
            return path;
        }

        if (path == path.root_path())
        {
            break;
        }
    }

    throw std::runtime_error("Failed to locate project root containing the performance directory.");
}

std::ofstream SimulationStats::openOutputFile(const std::filesystem::path &relativePath)
{
    const auto outputPath = findProjectRoot() / relativePath;

    std::error_code errorCode;
    std::filesystem::create_directories(outputPath.parent_path(), errorCode);
    if (errorCode)
    {
        throw std::runtime_error("Failed to create output directory: " +
                                 outputPath.parent_path().string());
    }

    std::ofstream output(outputPath);
    if (!output.is_open())
    {
        throw std::runtime_error("Failed to open output file: " + outputPath.string());
    }

    return output;
}

void SimulationStats::writeStaticPerAlgorithmFile(const std::string &runDir,
                                                  const StaticAlgResult &r, long particleCount)
{
    if (r.skipped) return;  // nothing measured, nothing to log for this algorithm/size
    std::filesystem::path dir = std::filesystem::path(runDir) / r.name;
    std::filesystem::create_directories(dir);
    std::ofstream out(dir / (std::to_string(particleCount) + ".txt"));
    double total = r.constructSec + r.querySec;
    out << r.constructSec << ' ' << r.querySec << ' ' << total << "\n";
}

void SimulationStats::appendStaticCsvRow(std::ofstream &csv, long particleCount,
                                         unsigned int nFluid, unsigned int nBoundary,
                                         const StaticAlgResult &r, double maxRssMb)
{
    csv << particleCount << ',' << nFluid << ',' << nBoundary << ',' << r.name << ','
        << r.indexingScheme << ',' << r.sortingMethod << ',';
    if (r.skipped)
    {
        csv << "-,-,-,";
    }
    else
    {
        double total = r.constructSec + r.querySec;
        csv << r.constructSec << ',' << r.querySec << ',' << total << ',';
    }
    csv << maxRssMb << "\n";
}

bool SimulationStats::logStats(clock_t begin, clock_t &lastRecord, float loggingTime,
                               VectorParticles &vectorParticles, std::ofstream &Time,
                               std::ofstream &Error, std::ofstream &CFL, std::ofstream &Param,
                               bool loggingEnabled)
{
    if (!loggingEnabled) return true;

    clock_t current = clock() - begin;

    if (loggingTime >= 0.f && current / CLOCKS_PER_SEC > loggingTime)
    {
        return false;
    }

    if (current - lastRecord >= CLOCKS_PER_SEC / 200)
    {
        lastRecord = current;
        float densityError = computeDensityError(computeAverageDensity(vectorParticles));
        float cflNumber = computeCflNumber(findMaxSpeed(vectorParticles));

        Time << (float)current / CLOCKS_PER_SEC << "\n";
        Error << densityError << "\n";
        CFL << cflNumber << "\n";

        Param << NUM_PARTICLES_FLUID << "\n"
              << TIMESTEP << "\n"
              << GAS_STIFFNESS << "\n"
              << VISCOSITY << "\n"
              << SPACING << "\n"
              << densityError << "\n";
    }
    return true;
}

void SimulationStats::logTimeForNeighbourSearch(clock_t begin, int currentSimStep, int maxNumSteps,
                                                std::ofstream &time)
{
    (void)currentSimStep;
    (void)maxNumSteps;
    clock_t current = clock() - begin;
    time << (float)current / CLOCKS_PER_SEC << "\n";
    printf("Logged the time");
}

void SimulationStats::updateColorTest(VectorParticles &vectorParticles)
{
    Particle &p1 = vectorParticles[0];
    p1.updateColorParticleNeighbours(vectorParticles, focusedParticleColor, focusedNeighboursColor);

    Particle &p2 = vectorParticles[44];
    p2.updateColorParticleNeighbours(vectorParticles, focusedParticleColor, focusedNeighboursColor);

    Particle &p3 = vectorParticles[99];
    p3.updateColorParticleNeighbours(vectorParticles, focusedParticleColor, focusedNeighboursColor);
}

float SimulationStats::computeAverageVelocity(VectorParticles &vectorParticles)
{
    float speedSum = 0.f;
    for (Particle &particle : vectorParticles)
    {
        if (particle.type == ParticleType::Fluid)
        {
            speedSum += particle.computeVelocityMagnitude();
        }
    }
    float averageSpeed = speedSum / NUM_PARTICLES_FLUID_CURRENT;
    return averageSpeed;
}

float SimulationStats::computeBoundaryNeighbourRatio(const Particle &particle,
                                                     const VectorParticles &vectorParticles)
{
    if (particle.neighbours.empty()) return 0.f;

    int boundaryNeighbours = 0;
    for (int neighbourIndex : particle.neighbours)
    {
        if (vectorParticles[neighbourIndex].type == ParticleType::Boundary) boundaryNeighbours++;
    }

    return static_cast<float>(boundaryNeighbours) / static_cast<float>(particle.neighbours.size());
}

float SimulationStats::findMaxSpeed(VectorParticles &vectorParticles)
{
    float maxSpeed = 0.f;
    for (Particle &particle : vectorParticles)
    {
        if (particle.type == ParticleType::Fluid)
        {
            float speed = particle.computeVelocityMagnitude();
            if (speed > maxSpeed) maxSpeed = speed;
        }
    }
    return maxSpeed;
}

float SimulationStats::computeAverageDensity(VectorParticles &vectorParticles)
{
    float densitySum = 0.f;
    for (const Particle &particle : vectorParticles)
    {
        if (particle.type == ParticleType::Fluid)
        {
            densitySum += particle.density;
        }
    }
    float averageDensity = densitySum / (float)NUM_PARTICLES_FLUID_CURRENT;
    return averageDensity;
}

float SimulationStats::computeDensityError(float averageDensity)
{
    float densityErrorPercentage = 0.f;
    float densityError = (averageDensity - REST_DENSITY) / REST_DENSITY;
    densityError = densityError < 0 ? -densityError : densityError;
    densityErrorPercentage = densityError * 100.f;
    return densityErrorPercentage;
}

float SimulationStats::computeCflNumber(float maxSpeed)
{
    if (maxSpeed == 0.f) return 0.f;
    float cflNumber = (maxSpeed * TIMESTEP) / SPACING;
    return cflNumber;
}

float SimulationStats::computeCflThreshold(float cflNumber, float maxSpeed)
{
    if (maxSpeed == 0.f) return 0.f;
    float cflBoundary = cflNumber * (SPACING / maxSpeed);
    return cflBoundary;
}

bool SimulationStats::isEquilibrium(VectorParticles &vectorParticles, int frameCounter,
                                    float threshold)
{
    float speedMax = findMaxSpeed(vectorParticles);
    printf("Max speed: %.6f\n", speedMax);
    bool reached;
    if (frameCounter < 1000)
    {
        reached = false;
    }
    else
    {
        reached = speedMax < threshold;
    }
    return reached;
}

void SimulationStats::staticSceneAnalysis(VectorParticles &vectorParticles,
                                          sf::RenderWindow &window, sf::VertexArray &particlesVA,
                                          sf::Font &font)
{
    std::vector<long> particleCounts = {
        1000000, 2000000, 3000000, 4000000, 5000000,
        6000000, 7000000, 8000000, 9000000, 10000000,
    };

    const std::string runDir = "performance/static/4thTry";
    const std::string screenshotDir = runDir + "/screenshots";
    std::filesystem::create_directories(runDir);
    std::filesystem::create_directories(screenshotDir);

    std::ofstream csv(std::filesystem::path(runDir) / "results.csv");
    csv << "particleCount,nFluid,nBoundary,algorithm,indexingScheme,sortingMethod,constructSec,"
           "querySec,totalSec,maxRSS_MB\n";

    for (long targetN : particleCounts)
    {
        // One try/catch per particle count, so a failure at one size doesn't lose the rest of
        // the sweep. This can't catch a hard crash (SIGSEGV/SIGABRT) - only actual C++
        // exceptions like std::bad_alloc.
        try
        {
            int side = static_cast<int>(std::round(std::sqrt(static_cast<double>(targetN))));

            vectorParticles.clear();
            InitializeStaticWaterColumn(vectorParticles, NUM_BOUNDARY_LAYERS, side, side);

            SetSimulationView(window, vectorParticles);

            // Draw + display twice and pump events once before capturing: the first display()
            // issues the GL draw calls, but on macOS the compositor needs the window to have
            // processed at least one event loop pass before copyToImage() reliably reads back
            // a non-blank framebuffer.
            for (int i = 0; i < 2; ++i)
            {
                window.clear(lightGray);
                DrawScene(window, vectorParticles, particlesVA, font, window.getView(),
                         static_cast<int>(window.getSize().x),
                         static_cast<int>(window.getSize().y), /*paused=*/true,
                         /*recording=*/false, /*logging=*/false, /*frameCounter=*/0,
                         /*showParameters=*/false);
                window.display();
                while (const std::optional event = window.pollEvent()) {}
            }

            sf::Texture tex(window.getSize());
            tex.update(window);
            sf::Image screenshot = tex.copyToImage();
            (void)screenshot.saveToFile(screenshotDir + "/" + std::to_string(targetN) + ".png");

            std::vector<StaticAlgResult> results;
            Timer timer;

            // IndexSortZCurveInsertion is deliberately excluded here - far slower than every
            // other combo at these particle counts and not worth the wait. A '-' placeholder is
            // logged for it below instead of a real measurement.
            struct IndexCombo
            {
                std::string label;
                IndexingScheme scheme;
                std::string schemeLabel;
                SortingMethod method;
                std::string methodLabel;
            };
            std::vector<IndexCombo> combos = {
                {"IndexSortLinearInsertion", IndexingScheme::LinearCurve, "LinearCurve",
                 SortingMethod::InsertionSort, "InsertionSort"},
                {"IndexSortLinearRadix", IndexingScheme::LinearCurve, "LinearCurve",
                 SortingMethod::RadixSort, "RadixSort"},
                {"IndexSortZCurveRadix", IndexingScheme::ZCurve, "ZCurve",
                 SortingMethod::RadixSort, "RadixSort"},
            };

            for (const auto &combo : combos)
            {
                IndexSort indexSort(combo.scheme, combo.method);

                timer.reset();
                indexSort.construct(vectorParticles, combo.scheme, combo.method);
                double constructSec = timer.elapsedSeconds();

                timer.reset();
                indexSort.query(vectorParticles, combo.scheme);
                double querySec = timer.elapsedSeconds();

                results.push_back(
                    {combo.label, combo.schemeLabel, combo.methodLabel, constructSec, querySec});
                std::cout << "Timed " << combo.label << " for " << targetN << " particles\n";
            }

            results.push_back(
                {"IndexSortZCurveInsertion", "ZCurve", "InsertionSort", 0.0, 0.0, /*skipped=*/true});
            std::cout << "Skipped IndexSortZCurveInsertion for " << targetN << " particles\n";

            // Spatial hashing: table size m = 2 * total particle count (fluid + boundary).
            unsigned int totalParticles = NUM_PARTICLES_FLUID + NUM_PARTICLES_BOUNDARY;
            SpatialHashing spatialHashing(Math::nextPrime(totalParticles * 2));
            timer.reset();
            spatialHashing.construct(vectorParticles);
            double shConstruct = timer.elapsedSeconds();
            timer.reset();
            spatialHashing.query(vectorParticles);
            double shQuery = timer.elapsedSeconds();
            results.push_back({"SpatialHashing", "-", "-", shConstruct, shQuery});
            std::cout << "Timed SpatialHashing for " << targetN << " particles\n";

            // Compact hashing: unordered_map-backed, no fixed table size to reason about.
            CompactHashing compactHashing;
            timer.reset();
            compactHashing.construct(vectorParticles);
            double chConstruct = timer.elapsedSeconds();
            timer.reset();
            compactHashing.query(vectorParticles);
            double chQuery = timer.elapsedSeconds();
            results.push_back({"CompactHashing", "-", "-", chConstruct, chQuery});
            std::cout << "Timed CompactHashing for " << targetN << " particles\n";

            struct rusage usage;
            getrusage(RUSAGE_SELF, &usage);
            // macOS reports ru_maxrss in bytes (Linux reports KB) - this is built/run on macOS.
            double maxRssMb = usage.ru_maxrss / (1024.0 * 1024.0);

            for (const auto &r : results)
            {
                writeStaticPerAlgorithmFile(runDir, r, targetN);
                appendStaticCsvRow(csv, targetN, NUM_PARTICLES_FLUID, NUM_PARTICLES_BOUNDARY, r,
                                   maxRssMb);
            }
            csv.flush();

            std::cout << "Done: " << targetN << " particles (" << NUM_PARTICLES_FLUID
                      << " fluid, " << NUM_PARTICLES_BOUNDARY << " boundary), maxRSS "
                      << maxRssMb << " MB\n";
        }
        catch (const std::bad_alloc &e)
        {
            csv << targetN << ",,,FAILED,bad_alloc,,,,,\n";
            csv.flush();
            std::cout << "FAILED,bad_alloc," << targetN << "\n";
        }
        catch (const std::exception &e)
        {
            csv << targetN << ",,,FAILED,exception,,,,," << e.what() << "\n";
            csv.flush();
            std::cout << "FAILED,exception," << targetN << "," << e.what() << "\n";
        }
    }

    std::cout << "All particle counts processed. Exiting static scene analysis.\n";
}

void SimulationStats::dynamicSceneAnalysis(VectorParticles &vectorParticles,
                                           sf::RenderWindow &window, sf::VertexArray &particlesVA,
                                           sf::Font &font)
{
    (void)particlesVA;
    (void)font;

    // Only the two radix-sorted index sort combos are measured - insertion sort and both
    // hashing methods were already shown to be uncompetitive at these scales in the static
    // sweep.
    std::vector<long> particleCounts = {1000000, 2000000, 3000000, 5000000};
    const int numSteps = 1000;

    struct DynCombo
    {
        std::string label;
        IndexingScheme scheme;
        SortingMethod method;
    };
    std::vector<DynCombo> combos = {
        {"IndexSortLinearRadix", IndexingScheme::LinearCurve, SortingMethod::RadixSort},
        {"IndexSortZCurveRadix", IndexingScheme::ZCurve, SortingMethod::RadixSort},
    };

    const std::string runDir = "performance/dynamic/scaling";
    std::filesystem::create_directories(runDir);

    std::ofstream csv(std::filesystem::path(runDir) / "results.csv");
    csv << "particleCount,algorithm,step,constructSec,querySec,physicsSec,fullStepSec\n";

    ParticleSimulation particleSimulation;

    for (long targetN : particleCounts)
    {
        // u=250 gives exactly 1,000,000 fluid particles (16*u*u) - see
        // InitializeBreakingDamWideLow's parametrized overload.
        int u = static_cast<int>(std::round(std::sqrt(targetN / 16.0)));

        for (const auto &combo : combos)
        {
            try
            {
                vectorParticles.clear();
                InitializeBreakingDamWideLow(vectorParticles, NUM_BOUNDARY_LAYERS, u);

                IndexSort indexSort(combo.scheme, combo.method);
                IndexSort zIndexSortUnused(combo.scheme, combo.method);
                SpatialHashing spatialHashingUnused(1);
                CompactHashing compactHashingUnused;

                Timer timer;
                std::cout << "Starting " << combo.label << " for " << targetN
                          << " fluid particles (u=" << u << ", NUM_PARTICLES_FLUID=" << NUM_PARTICLES_FLUID
                          << ")\n";

                for (int step = 0; step < numSteps; ++step)
                {
                    timer.reset();
                    indexSort.construct(vectorParticles, combo.scheme, combo.method);
                    double constructSec = timer.elapsedSeconds();

                    timer.reset();
                    indexSort.query(vectorParticles, combo.scheme);
                    double querySec = timer.elapsedSeconds();

                    timer.reset();
                    particleSimulation.simulationStep(
                        vectorParticles, compactHashingUnused, spatialHashingUnused, indexSort,
                        zIndexSortUnused,
                        SceneGeometry::computeMinsAndMaxs(vectorParticles, NUM_PARTICLES_FLUID_CURRENT),
                        NUM_BOUNDARY_LAYERS, /*neighbourSearch=*/1, /*swingGravity=*/false);
                    double physicsSec = timer.elapsedSeconds();

                    double fullStepSec = constructSec + querySec + physicsSec;

                    csv << targetN << ',' << combo.label << ',' << step << ',' << constructSec
                        << ',' << querySec << ',' << physicsSec << ',' << fullStepSec << "\n";

                    if (step % 100 == 0)
                    {
                        csv.flush();
                        std::cout << "  step " << step << "/" << numSteps << " construct="
                                  << constructSec << " query=" << querySec << " physics="
                                  << physicsSec << "\n";
                    }
                }
                csv.flush();
                std::cout << "Done: " << combo.label << " for " << targetN << " particles\n";
            }
            catch (const std::exception &e)
            {
                csv << targetN << ',' << combo.label << ",FAILED,,,,\n";
                csv.flush();
                std::cout << "FAILED," << combo.label << "," << targetN << "," << e.what()
                          << "\n";
            }
        }
    }

    std::cout << "All particle counts processed. Exiting dynamic scene analysis.\n";
}

void SimulationStats::neighbourCompressionAnalysis(VectorParticles &vectorParticles)
{
    std::vector<long> particleCounts = {
        100000, 500000, 1000000, 3000000, 5000000, 10000000,
    };

    const std::string runDir = "performance/static/compression";
    std::filesystem::create_directories(runDir);
    std::ofstream csv(std::filesystem::path(runDir) / "results.csv");
    csv << "particleCount,nFluid,nBoundary,rawBytes,compressedBytes,compressionRatio\n";

    for (long targetN : particleCounts)
    {
        try
        {
            int side = static_cast<int>(std::round(std::sqrt(static_cast<double>(targetN))));

            vectorParticles.clear();
            InitializeStaticWaterColumn(vectorParticles, NUM_BOUNDARY_LAYERS, side, side);

            IndexSort indexSort(IndexingScheme::ZCurve, SortingMethod::RadixSort);
            indexSort.construct(vectorParticles, IndexingScheme::ZCurve, SortingMethod::RadixSort);
            indexSort.query(vectorParticles, IndexingScheme::ZCurve);

            std::vector<NeighbourCompression::CompressedList> compressed =
                NeighbourCompression::CompressAll(vectorParticles, indexSort);
            NeighbourCompression::CompressionStats stats =
                NeighbourCompression::ComputeStats(compressed);

            double ratio = stats.rawBytes > 0
                              ? static_cast<double>(stats.compressedBytes) / stats.rawBytes
                              : 0.0;

            csv << targetN << ',' << NUM_PARTICLES_FLUID << ',' << NUM_PARTICLES_BOUNDARY << ','
                << stats.rawBytes << ',' << stats.compressedBytes << ',' << ratio << '\n';
            csv.flush();

            std::cout << "N=" << targetN << " (fluid=" << NUM_PARTICLES_FLUID
                      << "): raw=" << stats.rawBytes << "B compressed=" << stats.compressedBytes
                      << "B ratio=" << ratio << " (" << (100.0 * (1.0 - ratio))
                      << "% smaller)\n";
        }
        catch (const std::exception &e)
        {
            csv << targetN << ",,,,,\n";
            csv.flush();
            std::cout << "FAILED at " << targetN << " particles: " << e.what() << "\n";
        }
    }

    std::cout << "All particle counts processed. Exiting neighbour compression analysis.\n";
}

// Memory-only comparison of SpatialHashing vs CompactHashing at the same particle counts, using
// each structure's own analytic estimatedMemoryBytes() (SpatialHashing/CompactHashing in
// NeighbourSearch.h) rather than OS-reported RSS, which becomes unreliable once the working set
// exceeds physical RAM (paged-out memory undercounts the real peak). Writes
// performance/static/4thTry/memoryComparison/results.csv.
void SimulationStats::memoryComparisonAnalysis(VectorParticles &vectorParticles)
{
    std::vector<long> particleCounts = {
        1000000, 2000000, 3000000, 4000000, 5000000,
        6000000, 7000000, 8000000, 9000000, 10000000, 20000000, 30000000,
    };

    const std::string runDir = "performance/static/4thTry/memoryComparison";
    std::filesystem::create_directories(runDir);
    std::ofstream csv(std::filesystem::path(runDir) / "results.csv");
    csv << "particleCount,algorithm,estimatedMemoryMB\n";

    for (long targetN : particleCounts)
    {
        try
        {
            int side = static_cast<int>(std::round(std::sqrt(static_cast<double>(targetN))));

            vectorParticles.clear();
            InitializeStaticWaterColumn(vectorParticles, NUM_BOUNDARY_LAYERS, side, side);

            unsigned int totalParticles = NUM_PARTICLES_FLUID + NUM_PARTICLES_BOUNDARY;

            SpatialHashing spatialHashing(Math::nextPrime(totalParticles * 2));
            spatialHashing.construct(vectorParticles);
            double spatialMb = spatialHashing.estimatedMemoryBytes() / (1024.0 * 1024.0);
            csv << targetN << ",SpatialHashing," << spatialMb << "\n";
            csv.flush();
            std::cout << "SpatialHashing " << targetN << " particles: " << spatialMb << " MB\n";

            CompactHashing compactHashing;
            compactHashing.construct(vectorParticles);
            double compactMb = compactHashing.estimatedMemoryBytes() / (1024.0 * 1024.0);
            csv << targetN << ",CompactHashing," << compactMb << "\n";
            csv.flush();
            std::cout << "CompactHashing " << targetN << " particles: " << compactMb << " MB\n";
        }
        catch (const std::exception &e)
        {
            std::cout << "FAILED," << targetN << "," << e.what() << "\n";
        }
    }

    std::cout << "Memory comparison sweep complete.\n";
}
