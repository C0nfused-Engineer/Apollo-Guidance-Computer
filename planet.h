#pragma once
#include "raylib.h"
#include "raymath.h"
#include <vector>
#include <string>
#include <cmath>

// Double-precision 3D vector for world-space positions in real units (meters)
struct DVec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    DVec3 operator+(const DVec3& o) const { return { x+o.x, y+o.y, z+o.z }; }
    DVec3 operator-(const DVec3& o) const { return { x-o.x, y-o.y, z-o.z }; }

    double length() const { return sqrt(x*x + y*y + z*z); }

    static DVec3 zero() { return {0.0, 0.0, 0.0}; }

    // Convert to raylib float Vector3, relative to an origin.
    // Always subtract a reference point before casting to avoid float precision loss.
    Vector3 toVec3(const DVec3& origin = zero()) const {
        return { (float)(x - origin.x), (float)(y - origin.y), (float)(z - origin.z) };
    }
};

class Planet {
public:
    std::string name;
    double      mass;
    double      radius;     // meters
    double      soiRadius;  // meters

    // Keplerian orbital elements (SI: meters, radians)
    double semiMajorAxis;
    double eccentricity;
    double inclination;
    double argOfPeriapsis;
    double lonAscendingNode;
    double meanAnomalyEpoch;

    Color color = WHITE;

    Planet* parent;
    std::vector<Planet*> children;

    // World position in meters
    DVec3 position;

    static constexpr double G = 6.674e-11;

    Planet(std::string  name,
           double       mass,
           double       radius,
           double       semiMajorAxis,
           double       eccentricity,
           double       inclination,
           double       argOfPeriapsis,
           double       lonAscendingNode,
           double       meanAnomalyEpoch,
           Planet*      parent = nullptr)
        : name(name),
          mass(mass),
          radius(radius),
          semiMajorAxis(semiMajorAxis),
          eccentricity(eccentricity),
          inclination(inclination),
          argOfPeriapsis(argOfPeriapsis),
          lonAscendingNode(lonAscendingNode),
          meanAnomalyEpoch(meanAnomalyEpoch),
          parent(parent)
    {
        if (parent) {
            parent->children.push_back(this);
            soiRadius = semiMajorAxis * pow(mass / parent->mass, 2.0 / 5.0);
        } else {
            soiRadius = 1e30;
        }
    }

    void update(double t) {
        if (parent)
            position = computePosition(t);
        for (auto* child : children)
            child->update(t);
    }

    Planet* getSOI(const DVec3& worldPos) {
        for (auto* child : children) {
            DVec3 delta = { worldPos.x - child->position.x,
                            worldPos.y - child->position.y,
                            worldPos.z - child->position.z };
            if (delta.length() < child->soiRadius)
                return child->getSOI(worldPos);
        }
        return this;
    }

    // camPos: world-space camera position in meters.
    // All geometry is shifted relative to camPos before converting to float,
    // keeping rendering near the float origin regardless of world scale.
    void draw(const DVec3& camPos) const {
        Vector3 relPos = position.toVec3(camPos);
        DrawSphereWires(relPos, (float)radius, 12, 12, color);
        drawOrbit(camPos);
        for (auto* child : children)
            child->draw(camPos);
    }

private:
    double solveKepler(double M, double e) const {
        M = fmod(M, 2.0 * M_PI);
        if (M < 0) M += 2.0 * M_PI;
        double E = (e < 0.8) ? M : M_PI;
        for (int i = 0; i < 100; i++) {
            double dE = (M - E + e * sin(E)) / (1.0 - e * cos(E));
            E += dE;
            if (fabs(dE) < 1e-10) break;
        }
        return E;
    }

    DVec3 computePosition(double t) const {
        double n  = sqrt(G * parent->mass / pow(semiMajorAxis, 3.0));
        double M  = meanAnomalyEpoch + n * t;
        double E  = solveKepler(M, eccentricity);

        double nu = 2.0 * atan2(
            sqrt(1.0 + eccentricity) * sin(E / 2.0),
            sqrt(1.0 - eccentricity) * cos(E / 2.0)
        );

        double r     = semiMajorAxis * (1.0 - eccentricity * cos(E));
        double x_orb = r * cos(nu);
        double y_orb = r * sin(nu);

        double cosO = cos(lonAscendingNode), sinO = sin(lonAscendingNode);
        double cosi = cos(inclination),      sini = sin(inclination);
        double cosw = cos(argOfPeriapsis),   sinw = sin(argOfPeriapsis);

        double x = (cosO * cosw - sinO * sinw * cosi) * x_orb
                 + (-cosO * sinw - sinO * cosw * cosi) * y_orb;
        double y = (sinO * cosw + cosO * sinw * cosi) * x_orb
                 + (-sinO * sinw + cosO * cosw * cosi) * y_orb;
        double z = (sinw * sini) * x_orb
                 + (cosw * sini) * y_orb;

        return {
            parent->position.x + x,
            parent->position.y + y,
            parent->position.z + z
        };
    }

    void drawOrbit(const DVec3& camPos) const {
        if (!parent) return;

        const int segments = 512;
        Vector3 prev = {};

        for (int i = 0; i <= segments; i++) {
            double M     = (2.0 * M_PI / segments) * i;
            double E     = solveKepler(M, eccentricity);
            double nu    = 2.0 * atan2(
                sqrt(1.0 + eccentricity) * sin(E / 2.0),
                sqrt(1.0 - eccentricity) * cos(E / 2.0)
            );
            double r     = semiMajorAxis * (1.0 - eccentricity * cos(E));
            double x_orb = r * cos(nu);
            double y_orb = r * sin(nu);

            double cosO = cos(lonAscendingNode), sinO = sin(lonAscendingNode);
            double cosi = cos(inclination),      sini = sin(inclination);
            double cosw = cos(argOfPeriapsis),   sinw = sin(argOfPeriapsis);

            double x = (cosO * cosw - sinO * sinw * cosi) * x_orb
                     + (-cosO * sinw - sinO * cosw * cosi) * y_orb;
            double y = (sinO * cosw + cosO * sinw * cosi) * x_orb
                     + (-sinO * sinw + cosO * cosw * cosi) * y_orb;
            double z = (sinw * sini) * x_orb
                     + (cosw * sini) * y_orb;

            DVec3 worldPt = {
                parent->position.x + x,
                parent->position.y + y,
                parent->position.z + z
            };

            Vector3 pt = worldPt.toVec3(camPos);
            if (i > 0) DrawLine3D(prev, pt, DARKGRAY);
            prev = pt;
        }
    }
};