# Neighbor Search Acceleration Structures for SPH Fluid Simulation

Bachelor's thesis · Albert-Ludwigs-Universität Freiburg · Chair for Computer Graphics · 2026

> This repository is a snapshot of the thesis code, for reading and review only. It is not set up to be built or run.

A 2D fluid simulation (Smoothed-Particle Hydrodynamics) in C++ in which six neighbor search algorithms share one solver, so they can be compared fairly: index sort (linear or Z-curve indexing, each with insertion or radix sort), spatial hashing and compact hashing. The solver is parallelized with OpenMP.

**Main findings** (details in the [thesis](bachelorthesis_nahhas_final.pdf)):
- Index sort with radix sort was the most broadly effective variant. It clearly beat both hashing methods on static scenes (1 to 80 million particles), and index sort in general stayed cheaper than hashing in a running simulation.
- Insertion sort did not scale: 430.6 s vs. 4.0 s with radix sort for 10 million particles.
- Memory, not compute, was the limit: on the 8 GB test machine it ran out at around 36 million particles. A delta-encoding of the neighbor lists saved about 57.5 % of their memory.

**Demo videos:** [breaking dam](simulation_videos/breaking_dam_wide100k.mp4) · [winding channel](simulation_videos/downward_flow100k.mp4) · [tsunami](simulation_videos/tsunami100k.mp4) · [tsunami, swinging gravity](simulation_videos/tsunami_tilting50k.mp4) · [braided channel](simulation_videos/undergroundRed80k.mp4)

**Where to start reading:** `code/src/simulation/NeighbourSearch.cpp` (the search variants), `code/src/simulation/ParticleSimulation.cpp` (one solver step with OpenMP), and Sections 3 and 5 of the thesis.

**Tech:** C++17 · OpenMP · SFML · Python (NumPy, SciPy, OpenCV) for the sketch-to-scene tools in `code/tools/scene_import` · developed on macOS
