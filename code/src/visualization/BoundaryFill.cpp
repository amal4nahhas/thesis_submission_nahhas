/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#include "BoundaryFill.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

#include "../scene/SceneGeometry.h"

sf::VertexArray BoundaryFill::fillVA(sf::PrimitiveType::Triangles);

void BoundaryFill::Rebuild(VectorParticles &vectorParticles, unsigned int numParticlesFluid,
                           float spacing)
{
    fillVA.clear();
    if (spacing <= 0.f || vectorParticles.empty()) return;

    // Bounding box over ALL particles (boundaryOnly=false), not just boundary ones - an
    // enclosed empty pocket can sit anywhere inside the scene, and fluid cells need to be
    // classified too (see below).
    std::vector<float> bounds = SceneGeometry::computeMinsAndMaxs(
        vectorParticles, static_cast<int>(numParticlesFluid), /*boundaryOnly=*/false);
    const float minX = bounds[0], maxX = bounds[1], minY = bounds[2], maxY = bounds[3];
    if (maxX <= minX || maxY <= minY) return;

    // Grid resolution: finer than particle spacing so that consecutive boundary particles
    // along a traced curve always splat onto edge-connected cells rather than just
    // corner-touching ones - a wall that only touches itself at cell corners doesn't reliably
    // block a 4-connected flood fill. This is a rendering-only grid with no particles in it,
    // so it can be as fine as looks good without any physics implication (unlike particle
    // spacing itself, which must stay exactly SPACING - see Appendix, Minimum Particle
    // Spacing) - finer here just means smoother-looking fill edges.
    float cellSize = spacing * 0.2f;

    // Cap the total cell count so a scene with a huge bounding box relative to spacing (e.g.
    // the 1-10M-particle rectangular static-benchmark tanks) can't blow up memory/time -
    // coarsen the grid instead of skipping the feature outright. Those scenes are plain
    // rectangles with nothing enclosed anyway, so the lost precision costs nothing there.
    constexpr long long kMaxCells = 20'000'000;
    unsigned int cellsX = 0, cellsY = 0;
    while (true)
    {
        cellsX = static_cast<unsigned int>(std::ceil((maxX - minX) / cellSize)) + 1;
        cellsY = static_cast<unsigned int>(std::ceil((maxY - minY) / cellSize)) + 1;
        if (static_cast<long long>(cellsX) * static_cast<long long>(cellsY) <= kMaxCells) break;
        cellSize *= 1.5f;
    }

    std::vector<CellState> grid(static_cast<std::size_t>(cellsX) * cellsY, CellState::Empty);
    auto cellIndex = [&](int cx, int cy) -> std::size_t
    { return static_cast<std::size_t>(cy) * cellsX + cx; };
    auto toCell = [&](float v, float minV, unsigned int count) -> int
    {
        int c = static_cast<int>((v - minV) / cellSize);
        return std::clamp(c, 0, static_cast<int>(count) - 1);
    };

    // Splat each boundary particle onto every cell within a small radius, not just its own
    // containing cell, so a regular grid of particles spaced `spacing` apart - whether a
    // traced boundary curve or a fluid fill - becomes a continuous region in the grid with no
    // gaps between neighboring particles, regardless of angle relative to the grid axes.
    const float splatRadius = std::max(cellSize * 1.5f, spacing * 0.75f);
    const int splatCells = static_cast<int>(std::ceil(splatRadius / cellSize));
    auto splat = [&](const Particle &p, CellState state, bool overwriteFluid)
    {
        int cx = toCell(p.position.x, minX, cellsX);
        int cy = toCell(p.position.y, minY, cellsY);
        for (int dy = -splatCells; dy <= splatCells; dy++)
        {
            int ny = cy + dy;
            if (ny < 0 || ny >= static_cast<int>(cellsY)) continue;
            for (int dx = -splatCells; dx <= splatCells; dx++)
            {
                int nx = cx + dx;
                if (nx < 0 || nx >= static_cast<int>(cellsX)) continue;
                float wx = minX + (nx + 0.5f) * cellSize;
                float wy = minY + (ny + 0.5f) * cellSize;
                float ddx = wx - p.position.x, ddy = wy - p.position.y;
                if (ddx * ddx + ddy * ddy > splatRadius * splatRadius) continue;
                std::size_t idx = cellIndex(nx, ny);
                if (grid[idx] == CellState::Wall) continue;  // wall always wins
                if (grid[idx] == CellState::Fluid && !overwriteFluid) continue;
                grid[idx] = state;
            }
        }
    };

    // Boundary particles first, so a wall is never overwritten by a fluid splat below.
    for (std::size_t i = numParticlesFluid; i < vectorParticles.size(); i++)
        splat(vectorParticles[i], CellState::Wall, /*overwriteFluid=*/true);

    // Fluid particles too, so they are never treated as enclosed empty space - both so that a
    // scene whose boundary runs along the whole domain edge (e.g. the rectangular
    // static-benchmark tanks) has a Fluid, not Empty, interior, and so that the grid cells
    // *between* neighboring fluid particles (cellSize is finer than their spacing) don't read
    // as empty gaps.
    for (std::size_t i = 0; i < numParticlesFluid; i++)
        splat(vectorParticles[i], CellState::Fluid, /*overwriteFluid=*/false);

    // Flood fill through anything that isn't a Wall (Empty or Fluid), starting from every
    // Fluid cell and from the grid border, marking every reached cell as visited. Any Empty
    // cell left unvisited afterwards is fully enclosed by boundary particles and unreachable
    // from both the actual fluid and the outside of the scene.
    //
    // Seeding from Fluid matters for scenes like a long winding channel that starts full of
    // fluid at one end and empty at the other (nothing has flowed there yet on frame 0): the
    // scene's bounding box hugs the channel walls tightly, so the grid border alone barely
    // touches "outside" at all, and without this, the entire not-yet-reached corridor would
    // wrongly read as an enclosed void instead of the open flow path it actually is.
    std::vector<uint8_t> visited(grid.size(), 0);
    std::vector<std::size_t> frontier;
    frontier.reserve(grid.size() / 4 + 1);
    auto tryVisit = [&](int cx, int cy)
    {
        if (cx < 0 || cx >= static_cast<int>(cellsX) || cy < 0 || cy >= static_cast<int>(cellsY))
            return;
        std::size_t idx = cellIndex(cx, cy);
        if (visited[idx] || grid[idx] == CellState::Wall) return;
        visited[idx] = 1;
        frontier.push_back(idx);
    };
    for (std::size_t idx = 0; idx < grid.size(); idx++)
        if (grid[idx] == CellState::Fluid) tryVisit(static_cast<int>(idx % cellsX),
                                                     static_cast<int>(idx / cellsX));
    for (unsigned int cx = 0; cx < cellsX; cx++)
    {
        tryVisit(static_cast<int>(cx), 0);
        tryVisit(static_cast<int>(cx), static_cast<int>(cellsY) - 1);
    }
    for (unsigned int cy = 0; cy < cellsY; cy++)
    {
        tryVisit(0, static_cast<int>(cy));
        tryVisit(static_cast<int>(cellsX) - 1, static_cast<int>(cy));
    }
    for (std::size_t head = 0; head < frontier.size(); head++)
    {
        std::size_t idx = frontier[head];
        int cx = static_cast<int>(idx % cellsX);
        int cy = static_cast<int>(idx / cellsX);
        tryVisit(cx + 1, cy);
        tryVisit(cx - 1, cy);
        tryVisit(cx, cy + 1);
        tryVisit(cx, cy - 1);
    }

    // Build one quad per non-open cell - both the enclosed Empty pockets and the Wall cells
    // themselves. Filling Wall cells too (not just the enclosed space behind them) gives the
    // boundary particles a solid, gap-free backing: individual particle dots don't tile
    // perfectly (antialiased circles a full diameter apart can leave hairline gaps between
    // them along a curve), and without this fill those gaps showed the window's clear color
    // through as a visible seam right at the wall. Fluid and open-reachable Empty cells (i.e.
    // anything `visited`) are still left alone.
    for (unsigned int cy = 0; cy < cellsY; cy++)
    {
        for (unsigned int cx = 0; cx < cellsX; cx++)
        {
            std::size_t idx = cellIndex(static_cast<int>(cx), static_cast<int>(cy));
            if (grid[idx] == CellState::Fluid || visited[idx]) continue;

            float x0 = minX + cx * cellSize, y0 = minY + cy * cellSize;
            float x1 = x0 + cellSize, y1 = y0 + cellSize;
            const sf::Vertex quad[6] = {
                {{x0, y0}, boundaryColor}, {{x1, y0}, boundaryColor}, {{x1, y1}, boundaryColor},
                {{x0, y0}, boundaryColor}, {{x1, y1}, boundaryColor}, {{x0, y1}, boundaryColor},
            };
            for (const sf::Vertex &v : quad) fillVA.append(v);
        }
    }

    if (fillVA.getVertexCount() > 0)
    {
        std::cerr << "BoundaryFill: filled " << (fillVA.getVertexCount() / 6)
                  << " wall/enclosed cell(s), " << cellSize << " world units each.\n";
    }
}

void BoundaryFill::Draw(sf::RenderTarget &target)
{
    if (fillVA.getVertexCount() > 0) target.draw(fillVA);
}
