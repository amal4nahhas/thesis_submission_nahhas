/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include <cstdint>
#include <vector>

#include "NeighbourSearch.h"
#include "Particle.h"

// Thesis measurement tool: compresses each particle's neighbour list the way
// \cite{compressedNeighbourLists} does, to report the memory it saves over the raw
// std::vector<int> neighbour lists. Built from the existing neighbour lists after a query,
// entirely separate from them - the live simulation keeps using Particle::neighbours as-is,
// since decoding a delta stream on every access would slow down the actual force computation.
// Not wired into the simulation loop; call CompressAll after a query when you want to measure it.
namespace NeighbourCompression
{
// One compressed neighbour list: a sequence of entries, each a 1-byte control code optionally
// followed by data bytes, encoding the delta-compressed ranks as described below.
//
// Simplified from the paper: there, four 2-bit control codes are packed into a single byte to
// maximize memory savings. Here, each control code gets its own full byte instead - much
// simpler to implement and verify, at the cost of some of the achievable compression.
//
// Control code -> meaning -> extra data bytes:
//   0 -> delta == 0                -> none
//   1 -> delta == 1                -> none
//   2 -> delta in [2, 256)         -> 1 byte
//   3 -> delta >= 256              -> 4 bytes (little-endian)
using CompressedList = std::vector<uint8_t>;

// Total size of a set of compressed lists vs. the raw 4-bytes-per-neighbour representation
// they were built from, for reporting the achieved compression ratio.
struct CompressionStats
{
    size_t rawBytes = 0;
    size_t compressedBytes = 0;
};

/// @brief Build the rank of every particle in the sorted handles array, i.e. the inverse of
/// indexSort.handles: rankOfParticle[particleIndex] = rank. Compression needs a neighbour's
/// RANK in the cell-sorted order (small deltas only appear in that order), not its raw storage
/// index into vectorParticles.
/// @param indexSort a constructed (sorted) IndexSort
/// @param numParticles size of vectorParticles, so untouched (non-fluid) slots resolve safely
std::vector<int> BuildRankOfParticle(const IndexSort &indexSort, size_t numParticles);

/// @brief Compress one particle's neighbour list (raw indices into vectorParticles, as stored
/// in Particle::neighbours) into the delta-encoded format above. Neighbours are first mapped
/// to their rank via rankOfParticle, then sorted ascending, then delta-encoded relative to a
/// virtual previous rank of -1 (so the first neighbour's delta equals its own rank).
/// @param neighbours raw particle indices, e.g. Particle::neighbours
/// @param rankOfParticle from BuildRankOfParticle
CompressedList Compress(const std::vector<int> &neighbours, const std::vector<int> &rankOfParticle);

/// @brief Inverse of Compress: decodes a compressed list back into the sorted ranks it was
/// built from. Returns ranks, not particle indices - resolving a rank back to a particle index
/// is indexSort.handles[rank].particleIndex. Exists to verify Compress is correct and to
/// measure decode cost, not used on the live simulation path.
std::vector<int> Decompress(const CompressedList &compressed);

/// @brief Compress every particle's current neighbour list in one call, building
/// rankOfParticle once and reusing it for all of them.
/// @param vectorParticles must already have up-to-date Particle::neighbours (i.e. query() was
/// just called)
/// @param indexSort the IndexSort used for that query, providing the sorted-rank order
std::vector<CompressedList> CompressAll(const VectorParticles &vectorParticles, const IndexSort &indexSort);

/// @brief Sum the raw (4 bytes/neighbour) and compressed sizes of a set of compressed lists,
/// for reporting the achieved compression ratio.
CompressionStats ComputeStats(const std::vector<CompressedList> &compressedLists);
}  // namespace NeighbourCompression
