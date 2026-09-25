/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#include "NeighbourCompression.h"

#include <algorithm>

namespace NeighbourCompression
{
std::vector<int> BuildRankOfParticle(const IndexSort &indexSort, size_t numParticles)
{
    std::vector<int> rankOfParticle(numParticles, -1);
    for (int rank = 0; rank < static_cast<int>(indexSort.handles.size()); ++rank)
    {
        rankOfParticle[indexSort.handles[rank].particleIndex] = rank;
    }
    return rankOfParticle;
}

CompressedList Compress(const std::vector<int> &neighbours, const std::vector<int> &rankOfParticle)
{
    std::vector<int> ranks;
    ranks.reserve(neighbours.size());
    for (int particleIndex : neighbours) ranks.push_back(rankOfParticle[particleIndex]);
    std::sort(ranks.begin(), ranks.end());

    CompressedList compressed;
    int previousRank = -1;
    for (int rank : ranks)
    {
        int delta = rank - previousRank - 1;
        if (delta == 0)
        {
            compressed.push_back(0);
        }
        else if (delta == 1)
        {
            compressed.push_back(1);
        }
        else if (delta < 256)
        {
            compressed.push_back(2);
            compressed.push_back(static_cast<uint8_t>(delta));
        }
        else
        {
            compressed.push_back(3);
            compressed.push_back(static_cast<uint8_t>(delta & 0xFF));
            compressed.push_back(static_cast<uint8_t>((delta >> 8) & 0xFF));
            compressed.push_back(static_cast<uint8_t>((delta >> 16) & 0xFF));
            compressed.push_back(static_cast<uint8_t>((delta >> 24) & 0xFF));
        }
        previousRank = rank;
    }
    return compressed;
}

std::vector<int> Decompress(const CompressedList &compressed)
{
    std::vector<int> ranks;
    int previousRank = -1;
    size_t i = 0;
    while (i < compressed.size())
    {
        uint8_t code = compressed[i++];
        int delta;
        if (code == 0)
        {
            delta = 0;
        }
        else if (code == 1)
        {
            delta = 1;
        }
        else if (code == 2)
        {
            delta = compressed[i];
            i += 1;
        }
        else  // code == 3
        {
            delta = static_cast<int>(compressed[i]) | (static_cast<int>(compressed[i + 1]) << 8) |
                    (static_cast<int>(compressed[i + 2]) << 16) |
                    (static_cast<int>(compressed[i + 3]) << 24);
            i += 4;
        }
        int rank = previousRank + delta + 1;
        ranks.push_back(rank);
        previousRank = rank;
    }
    return ranks;
}

std::vector<CompressedList> CompressAll(const VectorParticles &vectorParticles, const IndexSort &indexSort)
{
    std::vector<int> rankOfParticle = BuildRankOfParticle(indexSort, vectorParticles.size());

    std::vector<CompressedList> result;
    result.reserve(vectorParticles.size());
    for (const Particle &particle : vectorParticles)
    {
        result.push_back(Compress(particle.neighbours, rankOfParticle));
    }
    return result;
}

CompressionStats ComputeStats(const std::vector<CompressedList> &compressedLists)
{
    CompressionStats stats;
    for (const CompressedList &compressed : compressedLists)
    {
        stats.compressedBytes += compressed.size();
        stats.rawBytes += Decompress(compressed).size() * sizeof(int);
    }
    return stats;
}
}  // namespace NeighbourCompression
