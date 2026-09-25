/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#pragma once

#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <vector>

extern sf::Color lightGray;
extern sf::Color gray;
extern sf::Color darkGray;
extern sf::Color black;
extern sf::Color white;
extern sf::Color blue;
extern sf::Color red;
extern sf::Color green;
extern sf::Color brown;
extern sf::Color yellow;
extern sf::Color azure;
extern sf::Color beigue;
extern sf::Color transparent;
extern sf::Color gold;
extern sf::Color fluidColor;
extern sf::Color boundaryColor;
extern sf::Color focusedParticleColor;
extern sf::Color focusedNeighboursColor;

// An enumerator of particle types.
enum ParticleType
{
    Fluid = 0,
    Boundary,
    Inactive
};

using V2f = sf::Vector2f;

struct Particle;
using VectorParticles = std::vector<Particle>;

struct CellIndices
{
   public:
    int k;
    int l;

    CellIndices(int k, int l) : k(k), l(l) {}
};

// Struct to represent a particle and its properties.
struct Particle
{
    unsigned int id; // mainly used for the neighbour search data structure.
    
    float mass;
    float density;
    float spacing;
    float radius;

    V2f position;
    V2f velocity;
    V2f accelerationPressure;
    V2f accelerationViscosity;
    V2f accelerationExternal;
    V2f accelerationTotal;
    V2f force;
    float pressure;

    std::vector<int> neighbours;
    CellIndices cellIndices;
    int cellIndex;

    sf::Color color;
    ParticleType type;

    Particle(float mass, float density, float spacing, ParticleType type, unsigned int id)
        : id(id),
          position(0.f, 0.f),
          velocity(0.f, 0.f),
          accelerationPressure(0.f, 0.f),
          accelerationViscosity(0.f, 0.f),
          accelerationExternal(0.f, 0.f),
          accelerationTotal(0.f, 0.f),
          force(0.f, 0.f),
          mass(mass),
          density(density),
          spacing(spacing),
          pressure(0.f),
          cellIndices{0, 0},
          cellIndex(0),
          type(type)
    {
        radius = spacing / 2.f;
        if (type == ParticleType::Fluid)
        {
            color = fluidColor;
        }
        else if (type == ParticleType::Boundary)
        {
            color = boundaryColor;
        }
    }

    /// @brief Translate a single particle by the given values.
    /// @param particle
    /// @param dx
    /// @param dy
    void translate(float dx, float dy);

    /// @brief Position a particle at (x, y) while accounting for particle center offset.
    /// @param x
    /// @param y
    void setPosition(float x, float y);

    /// @brief Rotate a single particle by the given angle in degrees.
    /// @param particle
    /// @param angleDegrees
    void rotate(float angleDegrees);

    /// @brief Update particle's color.
    /// @param particle
    /// @param color
    void updateColor(sf::Color color);

    /// @brief Update the neighbours' color.
    /// @param particle
    /// @param color
    void updateColorNeighbours(VectorParticles &vectorParticles, sf::Color color);

    /// @brief Update particle and neighbours' color.
    /// @param particle
    /// @param colorParticle
    /// @param colorNeighbours
    void updateColorParticleNeighbours(VectorParticles &vectorParticles, sf::Color colorParticle,
                                       sf::Color colorNeighbours);

    /// @brief Focus on a particle and its neighbours.
    /// @param vectorParticles
    /// @param particle
    void focus(VectorParticles &vectorParticles);

    /// @brief Focus on multiple particles and their neighbours.
    /// @param vectorParticles
    /// @param particlesToFocus
    void focusOnParticles(VectorParticles &vectorParticles, std::vector<int> particlesToFocus);

    /// @brief Compute the velocity magnitude of a particle.
    float computeVelocityMagnitude();

    /// @brief Compute the Euclidean distance between two particles.
    /// @param a
    /// @param b
    /// @return
    static float computeDistance(const Particle &a, const Particle &b);
};
