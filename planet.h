#pragma once
#include "raylib.h"
#include "raymath.h"
#include "Orbit.h"
#include <vector>
#include <string>
#include <cmath>

struct DVec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    DVec3 operator+(const DVec3& o) const { return { x+o.x, y+o.y, z+o.z }; }
    DVec3 operator-(const DVec3& o) const { return { x-o.x, y-o.y, z-o.z }; }

    double length() const { return std::sqrt(x*x + y*y + z*z); }

    static DVec3 zero() { return {0.0, 0.0, 0.0}; }

    Vector3 toVec3(const DVec3& origin = zero()) const {
        return { (float)(x-origin.x), (float)(y-origin.y), (float)(z-origin.z) };
    }
};

class Planet {
public:
    std::string name;
    double      mass;
    double      radius;    // meters
    double      soiRadius; // meters

    Orbit  orbit;   // ← all Keplerian elements live here now

    Color  color = WHITE;

    Planet* parent = nullptr;
    std::vector<Planet*> children;

    DVec3 position;

    static constexpr double G = 6.674e-11;

    Planet(std::string name_,
           double mass_,
           double radius_,
           double semiMajorAxis,
           double eccentricity,
           double inclination,
           double argOfPeriapsis,
           double lonAscendingNode,
           double meanAnomalyEpoch,
           Planet* parent_ = nullptr)
        : name(name_), mass(mass_), radius(radius_), parent(parent_)
    {
        orbit.sma   = semiMajorAxis;
        orbit.ecc   = eccentricity;
        orbit.inc   = inclination;
        orbit.aop   = argOfPeriapsis;
        orbit.lan   = lonAscendingNode;
        orbit.mae   = meanAnomalyEpoch;
        orbit.epoch = 0.0;
        orbit.valid = (parent_ != nullptr);

        if (parent_) {
            parent_->children.push_back(this);
            soiRadius = semiMajorAxis
                      * std::pow(mass_ / parent_->mass, 2.0/5.0);
        } else {
            soiRadius = 1e30;
        }
    }

    // Convenience accessors kept for legacy use
    double semiMajorAxis()    const { return orbit.sma; }
    double eccentricity()     const { return orbit.ecc; }
    double inclination()      const { return orbit.inc; }
    double argOfPeriapsis()   const { return orbit.aop; }
    double lonAscendingNode() const { return orbit.lan; }
    double meanAnomalyEpoch() const { return orbit.mae; }

    void update(double t) {
        if (parent)
            position = computePosition(t);
        for (auto* child : children)
            child->update(t);
    }

    Planet* getSOI(const DVec3& worldPos) {
        for (auto* child : children) {
            DVec3 d = { worldPos.x-child->position.x,
                        worldPos.y-child->position.y,
                        worldPos.z-child->position.z };
            if (d.length() < child->soiRadius)
                return child->getSOI(worldPos);
        }
        return this;
    }

    void draw(const DVec3& camPos) const {
        DrawSphereWires(position.toVec3(camPos),
                        (float)radius, 12, 12, color);
        // Draw orbit with Pe/Ap markers, using DARKGRAY for planets
        if (parent)
            orbit.draw(parent->position, camPos,
                       DARKGRAY, 512, false, 0.0);
        for (auto* child : children)
            child->draw(camPos);
    }

private:
    DVec3 computePosition(double t) const {
        double n  = std::sqrt(G * parent->mass
                              / std::pow(orbit.sma, 3.0));
        double M  = orbit.mae + n * t;
        double E  = orbit.solveKepler(M, orbit.ecc);

        double nu = 2.0 * std::atan2(
            std::sqrt(1.0+orbit.ecc) * std::sin(E/2.0),
            std::sqrt(1.0-orbit.ecc) * std::cos(E/2.0));

        double r    = orbit.sma * (1.0 - orbit.ecc*std::cos(E));
        double xOrb = r * std::cos(nu);
        double yOrb = r * std::sin(nu);

        double Rxx,Rxy,Ryx,Ryy,Rzx,Rzy;
        orbit.buildRotMatrix(Rxx,Rxy,Ryx,Ryy,Rzx,Rzy);

        return {
            parent->position.x + Rxx*xOrb + Rxy*yOrb,
            parent->position.y + Ryx*xOrb + Ryy*yOrb,
            parent->position.z + Rzx*xOrb + Rzy*yOrb
        };
    }
};