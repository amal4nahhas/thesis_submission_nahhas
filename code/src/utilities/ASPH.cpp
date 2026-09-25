/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 *
 * Thin driver for whichever SimulationStats sweep is currently being run - see
 * SimulationStats.cpp for what each one does (particle counts, algorithms, output layout).
 * Currently wired to neighbourCompressionAnalysis, which is headless (no window/font needed).
 * staticSceneAnalysis/dynamicSceneAnalysis draw a screenshot per particle count instead, and
 * would need the window/font/particlesVA setup they take as parameters.
 *
 * Usage: ./build/A
 ************************************************************/

#include "../scene/ParticleSettings.h"
#include "SimulationStats.h"

int NUM_BOUNDARY_LAYERS = 2;

int main()
{
    SimulationStats::neighbourCompressionAnalysis(vectorParticles);
    return 0;
}
