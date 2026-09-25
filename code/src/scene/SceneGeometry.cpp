/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/

#include "SceneGeometry.h"

#include <limits>

std::vector<float> SceneGeometry::computeMinsAndMaxs(VectorParticles &vectorParticles,
                                                     int numParticlesFluid, bool boundaryOnly)
{
    std::vector<float> minsAndMaxs(4);  // minX, maxX, minY, maxY
    minsAndMaxs[0] = std::numeric_limits<float>::max();
    minsAndMaxs[1] = std::numeric_limits<float>::lowest();
    minsAndMaxs[2] = std::numeric_limits<float>::max();
    minsAndMaxs[3] = std::numeric_limits<float>::lowest();

    int spacing = 0;
    int particleIndex = boundaryOnly ? numParticlesFluid : 0;
    for (int i = particleIndex; i < static_cast<int>(vectorParticles.size()); i++)
    {
        Particle &particle = vectorParticles[i];
        if (particle.position.x < minsAndMaxs[0]) minsAndMaxs[0] = particle.position.x;
        if (particle.position.x > minsAndMaxs[1]) minsAndMaxs[1] = particle.position.x;
        if (particle.position.y < minsAndMaxs[2]) minsAndMaxs[2] = particle.position.y;
        if (particle.position.y > minsAndMaxs[3]) minsAndMaxs[3] = particle.position.y;
    }

    minsAndMaxs[0] += spacing;
    minsAndMaxs[1] -= spacing;
    minsAndMaxs[2] += spacing;
    minsAndMaxs[3] -= spacing;
    return minsAndMaxs;
}

V2f SceneGeometry::computeWorldSize(VectorParticles &vectorParticles, int numParticlesFluid,
                                    bool boundaryOnly)
{
    std::vector<float> minsAndMaxs =
        computeMinsAndMaxs(vectorParticles, numParticlesFluid, boundaryOnly);
    float minX = minsAndMaxs[0];
    float maxX = minsAndMaxs[1];
    float minY = minsAndMaxs[2];
    float maxY = minsAndMaxs[3];
    float worldSizeX = (maxX - minX);
    float worldSizeY = (maxY - minY);
    return V2f(worldSizeX, worldSizeY);
}

V2f SceneGeometry::computeCenterPosition(VectorParticles &vectorParticles, int numParticlesFluid,
                                         bool boundaryOnly)
{
    std::vector<float> minsAndMaxs =
        computeMinsAndMaxs(vectorParticles, numParticlesFluid, boundaryOnly);
    float minX = minsAndMaxs[0];
    float maxX = minsAndMaxs[1];
    float minY = minsAndMaxs[2];
    float maxY = minsAndMaxs[3];
    float centerX = (minX + maxX) / 2.f;
    float centerY = (minY + maxY) / 2.f;
    return V2f(centerX, centerY);
}
