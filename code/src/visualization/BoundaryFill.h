/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include <SFML/Graphics.hpp>

#include <cstdint>

#include "../simulation/Particle.h"

// Fills the wall itself and any area fully enclosed by boundary particles (e.g. the interior of
// a hand-drawn rock or wall outline, where the boundary tracing in tools/scene_import only
// places particles along the outline itself) with the same color used for boundary particles -
// purely a rendering effect, no particles are added there. This also gives the boundary
// particle dots a solid, gap-free backing (individual antialiased dots don't tile perfectly, so
// without this a thin seam of background color could show through between them). Boundary
// particles never move once a scene is loaded (see InitializeFluidSimulation), so this only
// needs to run once per scene switch, not per frame.
class BoundaryFill
{
   public:
    /// @brief Rasterize the current boundary particles onto a grid, flood-fill from the scene's
    /// outer edge and from every fluid-occupied cell to find which cells are reachable (open)
    /// vs. not, and build a triangle mesh covering every non-open, non-fluid cell - the wall
    /// itself plus any pocket fully enclosed by it. Call once whenever the scene changes (after
    /// InitializeFluidSimulation), not every frame.
    /// @param vectorParticles
    /// @param numParticlesFluid particles [0, numParticlesFluid) are fluid, the rest are boundary
    /// @param spacing current scene's particle spacing (SPACING)
    static void Rebuild(VectorParticles &vectorParticles, unsigned int numParticlesFluid,
                        float spacing);

    /// @brief Draw the mesh built by the last Rebuild() call. No-op if nothing is enclosed.
    /// @param target
    static void Draw(sf::RenderTarget &target);

   private:
    enum class CellState : uint8_t
    {
        Empty = 0,
        Wall = 1,
        Fluid = 2,
    };

    static sf::VertexArray fillVA;
};
