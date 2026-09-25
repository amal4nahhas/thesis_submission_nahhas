/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once
#include <unordered_map>

#include "Particle.h"
#include "../scene/ParticleSettings.h"

extern const float PI;

class NeighbourSearch
{
   public:
    /// @brief Test whether two particles are neighbours according to the distance and the local
    /// radius support.
    /// @param particle1
    /// @param particle2
    /// @return
    static bool AreParticlesNeighbours(const Particle &particle1, const Particle &particle2);
};

class QuadraticNeighbourSearch : public NeighbourSearch
{
   public:
    /// @brief Find the neighbours of a particle by testing all other particles in the simulation.
    /// @param particle
    /// @param vectorParticles
    void findNeighbours(Particle &particle, VectorParticles &vectorParticles);

    void query(VectorParticles &vectorParticles);
};

class UniformGrid : public NeighbourSearch
{
   public:
    /// @brief Compute cell size, equal to or greater than the local radius support of the kernel
    /// function.
    /// @param spacing
    /// @return
    float ComputeCellSize(float spacing);

    /// @brief Map float coordinates to cell indices assuming a grid of cells with a given cell
    /// size.
    /// @param position
    /// @param cellSize
    /// @return
    CellIndices MapPosToCellIndices(const V2f &position, float cellSize = SPACING * 2.f);
};

struct Cell
{
    int cellIndex;
    CellIndices cellIndices;
    VectorParticles particlesInCell;
};

enum IndexingScheme
{
    LinearCurve,
    ZCurve
};

enum SortingMethod
{
    InsertionSort,
    RadixSort,
    stdSort
};

class Indexing
{
   public:
    static unsigned int computeCellIndex(CellIndices cellIndices, IndexingScheme indexingScheme);
    static unsigned int computeLinearIndex(CellIndices cellIndices);
    static unsigned int computeZIndex(CellIndices cellIndices);
    static unsigned int mortonCode(CellIndices cellIndices);
    static unsigned int spreadBits(unsigned int x);
};

struct Handle
{
   public:
    unsigned int particleIndex;
    int cellIndex;

    Handle(unsigned int particleIndex, int cellIndex) : particleIndex(particleIndex), cellIndex(cellIndex) {}
};

using Handles = std::vector<Handle>;

class Sorting
    {
    public:
        /// @brief Sort the handles by their cell indices using insertion sort (in-place).
        /// @param handles
        static const void InsertionSort(Handles &handles);
        /// @brief Radix sort (LSD) for 32-bit unsigned cell indices. Stable, two-pass (16-bit).
        /// @param handles
        static void RadixSort(Handles &handles);
    };

class IndexSort : public UniformGrid
{
   public:
    std::vector<int> cellReferences;
    Handles handles;

    IndexingScheme indexingScheme;
    SortingMethod sortingMethod;

    // Number of handles whose cell index changed during the most recent UpdateHandles() call
    // (0 after a full rebuild via AddParticlesToHandles, which doesn't track this).
    int lastChangedCellCount = 0;

    IndexSort(IndexingScheme indexingScheme, SortingMethod sortingMethod)
        : indexingScheme(indexingScheme), sortingMethod(sortingMethod) {};
    
    void AddParticlesToHandles(VectorParticles& vectorParticles);

    void UpdateHandles(VectorParticles& vectorParticles);

    void SortHandles();

    void FindCellReferences();

    void testCandidateNeighbours(Particle &particle, VectorParticles &vectorParticles);

    void findNeighbours(VectorParticles &vectorParticles);

    void resetCellReferences();

    void construct(VectorParticles &vectorParticles,
                   IndexingScheme indexingScheme,
                   SortingMethod sortingMethod);

    void construct(VectorParticles &vectorParticles);

    void query(VectorParticles &vectorParticles, IndexingScheme indexingScheme);

    void query(VectorParticles &vectorParticles);

    void reset()
    {
        handles.clear();
        resetCellReferences();
    };
};

enum HashingMethod
{
    Spatial,
    Compact
};

class Hashing
{
   public:
    /// @brief Compute the hash key for a given cell w.r.t. the size of the hash table.
    /// @param cell
    /// @return
    static unsigned int computeHashKey(CellIndices cellIndices, HashingMethod hashingMethod,
                                       unsigned int m = 0);
};

class SpatialHashing : public UniformGrid, Hashing
{
   public:
    // Size of the hash table.
    std::vector<std::vector<int>> hashTable;

    SpatialHashing(int m) : hashTable(m) {}

    void clearCells()
    {
        for (auto &cell : hashTable) cell.clear();
    }

    /// @brief Hash the particles to the hash grid.
    /// @param particle
    /// @param particleInd
    void hashParticle(Particle &particle, int particleInd);

    /// @brief Test whether the particles laying in the surrounding cells are neighbours, and if so,
    /// add them to the neighbours vector of the particle.
    /// @param particle
    /// @param vectorParticles
    void testCandidateNeighbours(Particle &particle, const VectorParticles &vectorParticles);

    /// @brief Update the neighbours for a given particle.
    /// @param particle
    /// @param vectorParticles
    void updateNeighbours(Particle &particle, const VectorParticles &vectorParticles);

    void construct(VectorParticles &vectorParticles);
    void query(VectorParticles &vectorParticles);

    /// @brief Analytic estimate of this structure's own heap footprint in bytes: the outer
    /// vector's control blocks (allocated for every cell whether or not it holds particles)
    /// plus each cell's own heap-allocated capacity. Deterministic, unlike getrusage()-based
    /// measurement (see SimulationStats::memoryComparisonAnalysis).
    size_t estimatedMemoryBytes() const
    {
        size_t bytes = hashTable.capacity() * sizeof(std::vector<int>);
        for (const auto &cell : hashTable)
        {
            bytes += cell.capacity() * sizeof(int);
        }
        return bytes;
    }
};

class CompactHashing : public UniformGrid, Hashing
{
   public:
    std::unordered_map<int, std::vector<int>> compactHashTable;

    /// @brief Hash the particles to the hash grid.
    /// @param particle
    /// @param particleInd
    void hashParticle(Particle &particle, int particleInd);

    /// @brief Test whether the particles laying in the surrounding cells are neighbours, and if so,
    /// add them to the neighbours vector of the particle.
    /// @param particle
    /// @param vectorParticles
    void testCandidateNeighbours(Particle &particle, const VectorParticles &vectorParticles);

    /// @brief Update the neighbours for a given particle.
    /// @param particle
    /// @param vectorParticles
    void updateNeighbours(Particle &particle, const VectorParticles &vectorParticles);

    void construct(VectorParticles &vectorParticles);
    void query(VectorParticles &vectorParticles);

    /// @brief Analytic estimate of this structure's own heap footprint in bytes: the bucket
    /// array plus one approximated node per occupied cell plus each occupied cell's own
    /// heap-allocated capacity. An approximation of the real libc++/libstdc++ node layout, but
    /// deterministic, unlike getrusage()-based measurement (see
    /// SimulationStats::memoryComparisonAnalysis).
    size_t estimatedMemoryBytes() const
    {
        struct ApproxNode
        {
            size_t hash;
            void *next;
            int key;
            std::vector<int> value;
        };
        size_t bytes = compactHashTable.bucket_count() * sizeof(void *);
        bytes += compactHashTable.size() * sizeof(ApproxNode);
        for (const auto &entry : compactHashTable)
        {
            bytes += entry.second.capacity() * sizeof(int);
        }
        return bytes;
    }
};

