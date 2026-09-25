/************************************************************
 * Project: SPH Fluid Solver in 2D
 * Author: Amal Nahhas <amal4nahhas@gmail.com>
 ************************************************************/
#include "Particle.h"

// Color definitions
sf::Color lightGray(211, 211, 211);
sf::Color gray(128, 128, 128);
sf::Color darkGray(64, 64, 64);
sf::Color black(0, 0, 0);
sf::Color white(255, 255, 255);
sf::Color blue(0, 145, 250);
sf::Color red(205, 78, 68);
sf::Color green(0, 135, 0);
sf::Color brown(150, 75, 0);
sf::Color yellow(210, 190, 0);
sf::Color azure(0, 170, 170);
sf::Color beigue(245, 245, 220);
sf::Color transparent(0, 0, 0, 0);
sf::Color gold(246, 202, 96);

sf::Color fluidColor = blue;
sf::Color boundaryColor = darkGray;
sf::Color focusedParticleColor = black;
sf::Color focusedNeighboursColor = red;

// __________________________________________________________________________________________________________________________________
void Particle::translate(float dx, float dy)
{
    position.x += dx;
    position.y += dy;
}

// ____________________________________________________________________________________________________________________________________
void Particle::setPosition(float x, float y)
{
    float offset = spacing / 2.f;
    position = V2f(x + offset, y + offset);
}

// ____________________________________________________________________________________________________________________________________
void Particle::rotate(float angleDegrees)
{
    float angleRadians = angleDegrees * 3.14159265f / 180.f;
    float cosAngle = cos(angleRadians);
    float sinAngle = sin(angleRadians);
    float xNew = position.x * cosAngle - position.y * sinAngle;
    float yNew = position.x * sinAngle + position.y * cosAngle;
    position.x = xNew;
    position.y = yNew;
}

// ____________________________________________________________________________________________________________________________________
void Particle::updateColor(sf::Color color) { this->color = color; }

// ____________________________________________________________________________________________________________________________________
void Particle::updateColorNeighbours(VectorParticles &vectorParticles, sf::Color color)
{
    for (int neighbour : neighbours)
    {
        vectorParticles[neighbour].updateColor(color);
    }
}

// ____________________________________________________________________________________________________________________________________
void Particle::updateColorParticleNeighbours(VectorParticles &vectorParticles,
                                             sf::Color colorParticle, sf::Color colorNeighbours)
{
    updateColor(colorParticle);
    updateColorNeighbours(vectorParticles, colorNeighbours);
}

// ____________________________________________________________________________________________________________________________________
void Particle::focus(VectorParticles &vectorParticles)
{
    updateColorNeighbours(vectorParticles, focusedNeighboursColor);
    updateColor(focusedParticleColor);
}

// ____________________________________________________________________________________________________________________________________
void Particle::focusOnParticles(VectorParticles &vectorParticles, std::vector<int> particlesToFocus)
{
    for (int targetParticle : particlesToFocus)
    {
        Particle &particleTarget = vectorParticles[targetParticle];
        particleTarget.updateColorNeighbours(vectorParticles, focusedNeighboursColor);
        particleTarget.updateColor(focusedParticleColor);
    }
}

// ____________________________________________________________________________________________________________________________________
float Particle::computeVelocityMagnitude()
{
    float velocityMagnitude = std::hypot(velocity.x, velocity.y);
    return velocityMagnitude;
}

// ____________________________________________________________________________________________________________________________________
float Particle::computeDistance(const Particle &a, const Particle &b)
{
    float dx = b.position.x - a.position.x;
    float dy = b.position.y - a.position.y;
    return std::sqrt(dx * dx + dy * dy);
}