/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#include "NeighbourSearch.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>
#include <iostream>
#include <random>

#include "../scene/ParticleSettings.h"

// ____________________________________________________________________________________________________________________________________
bool NeighbourSearch::AreParticlesNeighbours(const Particle &particle1, const Particle &particle2)
{
    float distance = Particle::computeDistance(particle1, particle2);
    bool areNeighbours = distance < 2.f * SPACING ? true : false;
    return areNeighbours;
}

// ____________________________________________________________________________________________________________________________________
void QuadraticNeighbourSearch::findNeighbours(Particle &particle, VectorParticles &vectorParticles)
{
    particle.neighbours.clear();
    for (int i = 0; i < vectorParticles.size(); i++)
    {
        Particle &candidate = vectorParticles[i];
        if (AreParticlesNeighbours(particle, candidate))
        {
            particle.neighbours.push_back(i);
        }
    }
}

// ____________________________________________________________________________________________________________________________________
void QuadraticNeighbourSearch::query(VectorParticles &vectorParticles)
{
    for (Particle &particle : vectorParticles)
    {
        if (particle.type == ParticleType::Fluid) findNeighbours(particle, vectorParticles);
    }
}

// ************************************************************************************************************************************
// Uniform Grid
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
float UniformGrid::ComputeCellSize(float spacing) { return SPACING * 2; }

// ____________________________________________________________________________________________________________________________________
CellIndices UniformGrid::MapPosToCellIndices(const V2f &position, float cellSize)
{
    // Pairwise floor function.
    int cellX = static_cast<int>(std::floor(position.x / cellSize));
    int cellY = static_cast<int>(std::floor(position.y / cellSize));
    return CellIndices{cellX, cellY};
}

// ************************************************************************************************************************************
// Indexing
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
unsigned int Indexing::computeLinearIndex(CellIndices cellIndices)
{
    auto K = NUM_BOUNDARY_PER_ROW;
    return cellIndices.k + cellIndices.l * K;
}

// ____________________________________________________________________________________________________________________________________
unsigned int Indexing::computeZIndex(CellIndices cellIndices) { return mortonCode(cellIndices); }

// ____________________________________________________________________________________________________________________________________
unsigned int Indexing::computeCellIndex(CellIndices cellIndices, IndexingScheme indexingScheme)
{
    if (indexingScheme == IndexingScheme::ZCurve)
        return computeZIndex(cellIndices);
    else
        return computeLinearIndex(cellIndices);

}

// ____________________________________________________________________________________________________________________________________
unsigned int Indexing::spreadBits(unsigned int x)
{
    x &= 0x0000ffff;  // keep only lower 16 bits
    x = (x | (x << 8)) & 0x00ff00ff;
    x = (x | (x << 4)) & 0x0f0f0f0f;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;
    return x;
}

// ____________________________________________________________________________________________________________________________________
unsigned int Indexing::mortonCode(CellIndices cellIndices)
{
    return spreadBits(cellIndices.k) | (spreadBits(cellIndices.l) << 1);
}

// ************************************************************************************************************************************
// Sorting
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
const void Sorting::InsertionSort(Handles &handles)
{
    for (size_t i = 1; i < handles.size(); i++)
    {
        Handle key = handles[i];
        int j = i - 1;
        while (j >= 0 && handles[j].cellIndex > key.cellIndex)
        {
            handles[j + 1] = handles[j];
            j--;
        }
        handles[j + 1] = key;
    }
}

// ____________________________________________________________________________________________________________________________________
void Sorting::RadixSort(Handles &handles)
{
    if (handles.empty()) return;

    const size_t n = handles.size();
    Handles aux;
    aux.resize(n, Handle(0, 0));

    // Two-pass LSD radix sort using 16-bit buckets (65536 buckets).
    const unsigned int MASK = 0xFFFFu;
    const size_t BUCKETS = 1u << 16;
    std::vector<int> count(BUCKETS);

    // Pass 1: sort by low 16 bits into aux
    std::fill(count.begin(), count.end(), 0);
    for (const Handle &h : handles) ++count[static_cast<unsigned int>(h.cellIndex) & MASK];
    for (size_t i = 1; i < BUCKETS; ++i) count[i] += count[i - 1];
    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        unsigned int key = static_cast<unsigned int>(handles[i].cellIndex) & MASK;
        aux[--count[key]] = handles[i];
    }

    // Pass 2: sort by high 16 bits back into handles
    std::fill(count.begin(), count.end(), 0);
    for (const Handle &h : aux) ++count[(static_cast<unsigned int>(h.cellIndex) >> 16) & MASK];
    for (size_t i = 1; i < BUCKETS; ++i) count[i] += count[i - 1];
    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        unsigned int key = (static_cast<unsigned int>(aux[i].cellIndex) >> 16) & MASK;
        handles[--count[key]] = aux[i];
    }
}

// *********************************************j ***************************************************************************************
// Index Sort
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
void IndexSort::AddParticlesToHandles(VectorParticles& vectorParticles)
{
    handles.clear();
    handles.resize(vectorParticles.size(), Handle(0, 0));
#pragma omp parallel for
    for (int i = 0; i < vectorParticles.size(); i++)
    {
        Particle &particle = vectorParticles[i];
        particle.cellIndices = UniformGrid::MapPosToCellIndices(particle.position, 2 * particle.spacing);
        particle.cellIndex = Indexing::computeCellIndex(particle.cellIndices, indexingScheme);
        Handle handle{static_cast<unsigned int>(i), particle.cellIndex};
        handles[i] = handle;
    }
}

void IndexSort::UpdateHandles(VectorParticles& vectorParticles)
{
    // NOTE: a scene switch (InputManager) must call IndexSort::reset() itself, right after
    // clearing vectorParticles.
    if (handles.size() == 0) AddParticlesToHandles(vectorParticles);
    else
    {
        int changedIndices = 0;
        std::vector<int> toErase;
    // Each thread gets a private partial count via reduction; toErase is mutated from
    // multiple threads, so its push_back is serialized below.
    #pragma omp parallel for reduction(+:changedIndices)
        for (int handleIndex = 0; handleIndex < handles.size(); handleIndex++)
        {
            Handle& handle = handles[handleIndex];
            Particle& particle = vectorParticles[handle.particleIndex];
            if (particle.type == ParticleType::Inactive)
            {
                #pragma omp critical(indexSortErase)
                {
                    toErase.push_back(handleIndex);
                }
                continue;
            }
            particle.cellIndices = UniformGrid::MapPosToCellIndices(particle.position, 2 * particle.spacing);
            particle.cellIndex = Indexing::computeCellIndex(particle.cellIndices, indexingScheme);
            bool cellIndexChanged = handle.cellIndex != particle.cellIndex;
            if (cellIndexChanged)
            {
                handle.cellIndex = particle.cellIndex;
                changedIndices++;
            }
        }
        lastChangedCellCount = changedIndices;
        // toErase isn't in ascending order (filled from parallel threads), so it must be sorted
        // before erasing back-to-front, otherwise erase() can receive an index past the
        // already-shrunk vector's size.
        std::sort(toErase.begin(), toErase.end());
        for (int i = toErase.size() - 1; i>=0; i--)
        {
            handles.erase(handles.begin() + toErase[i]);
        }
    }
}

// ____________________________________________________________________________________________________________________________________
void IndexSort::SortHandles()
{
    if (sortingMethod == SortingMethod::InsertionSort)
        Sorting::InsertionSort(handles);
    else if (sortingMethod == SortingMethod::stdSort)
        std::sort(handles.begin(), handles.end(),
                  [](const Handle &a, const Handle &b) { return a.cellIndex < b.cellIndex; });
    else if (sortingMethod == SortingMethod::RadixSort)
        Sorting::RadixSort(handles);
    else
        std::sort(handles.begin(), handles.end(),
                  [](const Handle &a, const Handle &b) { return a.cellIndex < b.cellIndex; });
}

// ____________________________________________________________________________________________________________________________________
void IndexSort::resetCellReferences()
{
    if (indexingScheme == IndexingScheme::LinearCurve)
        cellReferences.assign(NUM_BOUNDARY_PER_ROW * NUM_BOUNDARY_PER_COLUMN, -1);
    else if (indexingScheme == IndexingScheme::ZCurve)
        cellReferences.assign(Indexing::computeZIndex({static_cast<int>(NUM_BOUNDARY_PER_ROW),
                                                       static_cast<int>(NUM_BOUNDARY_PER_COLUMN)}),
                              -1);
}

// ____________________________________________________________________________________________________________________________________
void IndexSort::FindCellReferences()
{
    resetCellReferences();
#pragma omp parallel for
    for (int k = 0; k < handles.size(); k++)
    {
        int currentCellIndex = handles[k].cellIndex;
        bool isFirstInCell = (k == 0 || currentCellIndex != handles[k-1].cellIndex);
        if (isFirstInCell)
            cellReferences[currentCellIndex] = k;
    }
}

// ____________________________________________________________________________________________________________________________________
void IndexSort::testCandidateNeighbours(Particle &particle, VectorParticles &vectorParticles)
{
    particle.neighbours.clear();
    auto cellIndices = particle.cellIndices;
    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            auto cellIndicesNeighbour = CellIndices{cellIndices.k + dx, cellIndices.l + dy};
            unsigned int cellIndexNeighbour =
                Indexing::computeCellIndex(cellIndicesNeighbour, indexingScheme);
            if (cellIndexNeighbour == std::numeric_limits<unsigned int>::max() ||
                cellIndexNeighbour >= cellReferences.size())
                continue;
            int ref = cellReferences.at(cellIndexNeighbour);
            if (ref == -1) continue;  // Cell is empty
            for (int i = ref; i < handles.size(); i++)
            {
                Handle &handle = handles[i];
                int candidateIndex = handle.particleIndex;
                if (handle.cellIndex != cellIndexNeighbour) break;
                const Particle &candidate = vectorParticles[candidateIndex];
                bool isNeighbour = AreParticlesNeighbours(particle, candidate);
                if (isNeighbour) particle.neighbours.push_back(candidateIndex);
            }
        }
    }
}

// ____________________________________________________________________________________________________________________________________
void IndexSort::findNeighbours(VectorParticles &vectorParticles)
{
#pragma omp parallel for
    for (Particle &particle : vectorParticles)
    {
        if (particle.type == ParticleType::Fluid)
            testCandidateNeighbours(particle, vectorParticles);
    }
}

// ____________________________________________________________________________________________________________________________________
void IndexSort::construct(VectorParticles &vectorParticles, IndexingScheme indexingScheme,
                          SortingMethod sortingMethod)
{
    if (this->indexingScheme != indexingScheme) this->indexingScheme = indexingScheme;
    if (this->sortingMethod != sortingMethod) this->sortingMethod = sortingMethod;
    UpdateHandles(vectorParticles); // n
    SortHandles(); // n^2
    FindCellReferences(); // k
}

void IndexSort::query(VectorParticles &vectorParticles, IndexingScheme indexingScheme)
{
    if (this->indexingScheme != indexingScheme) this->indexingScheme = indexingScheme;
    findNeighbours(vectorParticles);
}

// ____________________________________________________________________________________________________________________________________
void IndexSort::construct(VectorParticles &vectorParticles)
{
    UpdateHandles(vectorParticles);
    SortHandles();
    FindCellReferences();
}

void IndexSort::query(VectorParticles &vectorParticles)
{
    findNeighbours(vectorParticles);
}
// ************************************************************************************************************************************
// Hashing
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
unsigned int Hashing::computeHashKey(CellIndices cellIndices, HashingMethod hashingMethod,
                                     unsigned int m)
{
    const unsigned int prime1 = 73856093;
    const unsigned int prime2 = 19349663;
    unsigned int key;
    if (hashingMethod == HashingMethod::Spatial)
    {
        // m is the fixed hash table size, only meaningful for Spatial (Compact has no fixed
        // table and never reaches this branch), so the mod-by-zero guard belongs here rather
        // than before the branch.
        if (m == 0) return 0;
        key = ((cellIndices.k * prime1) ^ (cellIndices.l * prime2)) % m;
    }
    else
    {
        key = (cellIndices.k * prime1) ^ (cellIndices.l * prime2);
    }
    return key;
}

// ************************************************************************************************************************************
// Spatial Hashing
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
void SpatialHashing::hashParticle(Particle &particle, int particleInd)
{
    float cellSize = ComputeCellSize(particle.spacing);
    particle.cellIndices = MapPosToCellIndices(particle.position, cellSize);
    int m = hashTable.size();
    unsigned int hashKey = computeHashKey(particle.cellIndices, HashingMethod::Spatial, m);
    // Two particles from different OpenMP threads (construct() runs this under a parallel
    // for) can collide into the same bucket; push_back on that shared std::vector<int> from
    // both threads at once is a data race that can corrupt the vector, so the actual
    // insertion is serialized here. Indexing hashTable[hashKey] itself stays lock-free since
    // the outer vector is fixed-size and never reallocates.
    #pragma omp critical(spatialHashInsert)
    {
        hashTable[hashKey].push_back(particleInd);
    }
}

// ____________________________________________________________________________________________________________________________________
void SpatialHashing::testCandidateNeighbours(Particle &particle,
                                             const VectorParticles &vectorParticles)
{
    CellIndices cellIndices = particle.cellIndices;
    auto m = hashTable.size();
    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            auto cellNeighbourIndices = CellIndices{cellIndices.k + dx, cellIndices.l + dy};
            unsigned int cellNeighbourKey =
                computeHashKey(cellNeighbourIndices, HashingMethod::Spatial, m);
            for (int candidateIndex : hashTable[cellNeighbourKey])
            {
                const Particle &candidate = vectorParticles[candidateIndex];
                bool isNeighbour = AreParticlesNeighbours(particle, candidate);
                if (isNeighbour) particle.neighbours.push_back(candidateIndex);
            }
        }
    }
}

// ____________________________________________________________________________________________________________________________________
void SpatialHashing::updateNeighbours(Particle &particle, const VectorParticles &vectorParticles)
{
    particle.neighbours.clear();
    testCandidateNeighbours(particle, vectorParticles);
}

void SpatialHashing::construct(VectorParticles &vectorParticles)
{
    // Clear all buckets but keep the hash table size the same.
    clearCells();
    #pragma omp parallel for
    for (int i = 0; i < vectorParticles.size(); i++)
    {
        if (vectorParticles[i].type != ParticleType::Inactive)
        {
            hashParticle(vectorParticles[i], i);
        }
    }
}

void SpatialHashing::query(VectorParticles &vectorParticles)
{
    #pragma omp parallel for
    for (int i = 0; i < vectorParticles.size(); i++)
    {
        Particle &particle = vectorParticles[i];
        if (particle.type == ParticleType::Fluid) updateNeighbours(particle, vectorParticles);
    }
}

// ************************************************************************************************************************************
// Compact Hashing
// ************************************************************************************************************************************

// ____________________________________________________________________________________________________________________________________
void CompactHashing::hashParticle(Particle &particle, int particleInd)
{
    float cellSize = ComputeCellSize(particle.spacing);
    particle.cellIndices = UniformGrid::MapPosToCellIndices(particle.position, cellSize);
    unsigned int hashKey = computeHashKey(particle.cellIndices, HashingMethod::Compact);
    // Unlike SpatialHashing's fixed-size table, compactHashTable is an unordered_map: inserting
    // a new key from one OpenMP thread while another thread reads/writes the map (construct()
    // runs this under a parallel for) can trigger a concurrent rehash, corrupting the map and
    // aborting the process. The whole map access is serialized here, not just the push_back.
    #pragma omp critical(compactHashInsert)
    {
        compactHashTable[hashKey].push_back(particleInd);
    }
}

// ____________________________________________________________________________________________________________________________________
void CompactHashing::testCandidateNeighbours(Particle &particle,
                                             const VectorParticles &vectorParticles)
{
    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            CellIndices cellIndices = particle.cellIndices;
            auto cellNeighbourIndices = CellIndices{cellIndices.k + dx, cellIndices.l + dy};
            unsigned int cellNeighbourKey =
                computeHashKey(cellNeighbourIndices, HashingMethod::Compact);
            // query() runs this under a parallel for with no writes to compactHashTable
            // happening concurrently, so a plain find() is race-free here - unlike
            // operator[], it never inserts a new entry for an empty/unseen cell, which
            // would otherwise mutate the map from multiple threads at once during a read.
            auto cellIt = compactHashTable.find(cellNeighbourKey);
            if (cellIt == compactHashTable.end()) continue;
            for (int candidateIndex : cellIt->second)
            {
                const Particle &candidate = vectorParticles[candidateIndex];
                bool isNeighbour = AreParticlesNeighbours(particle, candidate);
                if (isNeighbour) particle.neighbours.push_back(candidateIndex);
            }
        }
    }
}

// ____________________________________________________________________________________________________________________________________
void CompactHashing::updateNeighbours(Particle &particle, const VectorParticles &vectorParticles)
{
    particle.neighbours.clear();
    testCandidateNeighbours(particle, vectorParticles);
}

void CompactHashing::construct(VectorParticles &vectorParticles)
{
    compactHashTable.clear();
    #pragma omp parallel for
    for (int i = 0; i < vectorParticles.size(); i++)
    {
        if (vectorParticles[i].type != ParticleType::Inactive)
        {
            hashParticle(vectorParticles[i], i);
        }
    }
}

void CompactHashing::query(VectorParticles &vectorParticles)
{
    #pragma omp parallel for
    for (int i = 0; i < vectorParticles.size(); i++)
    {
        Particle &particle = vectorParticles[i];
        if (particle.type == ParticleType::Fluid) updateNeighbours(particle, vectorParticles);
    }
}

// ************************************************************************************************************************************
// Bounding Volme Hierarchies
// ************************************************************************************************************************************

class BinaryRadixTree
{
};