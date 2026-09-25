/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include <vector>

#include "../simulation/Particle.h"

class SceneGeometry
{
   public:
    static std::vector<float> computeMinsAndMaxs(VectorParticles &vectorParticles,
                                                 int numParticlesFluid, bool boundaryOnly = true);

    static V2f computeWorldSize(VectorParticles &vectorParticles, int numParticlesFluid,
                                bool boundaryOnly = true);

    static V2f computeCenterPosition(VectorParticles &vectorParticles, int numParticlesFluid,
                                     bool boundaryOnly = true);
};
