#pragma once
#include "raylib.h"
#include "raymath.h"
#include <vector>
#include <string>
#include <cmath>

class Planet {
public:
    // Physical properties
    std::string name;
    double mass;
    float radius;
    double soiRadius;

    // Orbital parameters (Keplerian elements)
    double semiMajorAxis;
    double eccentricity;
    double inclination;
    double argOfPeriapsis;
    double lonAscendingNode;
    double meanAnomalyEpoch;
    double renderScale = 1.0; // multiply computed positions by this for rendering

	// Visuals
	Color color = WHITE;

    // Hierarchy
    Planet* parent;
    std::vector<Planet*> children;

    // Computed each frame
    Vector3 position;

    static constexpr double G = 6.674e-11;

    Planet(std::string name,
           double mass,
           float radius,
           double semiMajorAxis,
           double eccentricity,
           double inclination,
           double argOfPeriapsis,
           double lonAscendingNode,
           double meanAnomalyEpoch,
           Planet* parent = nullptr,
           double renderScale = 1.0)
        : name(name),
          mass(mass),
          radius(radius),
          semiMajorAxis(semiMajorAxis),
          eccentricity(eccentricity),
          inclination(inclination),
          argOfPeriapsis(argOfPeriapsis),
          lonAscendingNode(lonAscendingNode),
          meanAnomalyEpoch(meanAnomalyEpoch),
          renderScale(renderScale),
          parent(parent),
          position({0,0,0})
    {
        if (parent) {
            parent->children.push_back(this);
            // SOI = a * (m / M)^(2/5), stored in render units
            soiRadius = semiMajorAxis * renderScale * pow(mass / parent->mass, 2.0 / 5.0);
        } else {
            soiRadius = 1e30; // root body SOI is effectively infinite
        }
    }

    void update(double t) {
        if (parent)
            position = computePosition(t);
        for (auto* child : children)
            child->update(t);
    }

    Planet* getSOI(Vector3 worldPos) {
        for (auto* child : children) {
            float dist = Vector3Distance(worldPos, child->position);
            if (dist < (float)child->soiRadius)
                return child->getSOI(worldPos);
        }
        return this;
    }

    void draw() const {
        DrawSphereWires(position, radius, 12, 12, color);
        drawOrbit();
        for (auto* child : children)
            child->draw();
    }

private:
    // Solve Kepler's equation M = E - e*sin(E) for E via Newton-Raphson
    double solveKepler(double M, double e) const {
        // Normalize M to [0, 2π]
        M = fmod(M, 2.0 * M_PI);
        if (M < 0) M += 2.0 * M_PI;

        // Initial guess — works well for low/moderate eccentricity
        double E = (e < 0.8) ? M : M_PI;

        for (int i = 0; i < 100; i++) {
            double dE = (M - E + e * sin(E)) / (1.0 - e * cos(E));
            E += dE;
            if (fabs(dE) < 1e-10) break;
        }

        return E;
    }

    Vector3 computePosition(double t) const {
        // 1. Mean motion (radians per second)
        double n = sqrt(G * parent->mass / pow(semiMajorAxis, 3.0));

        // 2. Mean anomaly at time t
        double M = meanAnomalyEpoch + n * t;

        // 3. Eccentric anomaly via Newton-Raphson
        double E = solveKepler(M, eccentricity);

        // 4. True anomaly
        double nu = 2.0 * atan2(
            sqrt(1.0 + eccentricity) * sin(E / 2.0),
            sqrt(1.0 - eccentricity) * cos(E / 2.0)
        );

        // 5. Distance from parent
        double r = semiMajorAxis * (1.0 - eccentricity * cos(E));

        // 6. Position in orbital plane (perifocal frame)
        double x_orb = r * cos(nu);
        double y_orb = r * sin(nu);

        // 7. Rotate into 3D world space using orbital elements
        //    Rz(-Ω) * Rx(-i) * Rz(-ω)
        double cosO = cos(lonAscendingNode), sinO = sin(lonAscendingNode);
        double cosi = cos(inclination),      sini = sin(inclination);
        double cosw = cos(argOfPeriapsis),   sinw = sin(argOfPeriapsis);

        // Standard perifocal -> ecliptic (XZ) frame
		// X = right, Y = up (north), Z = toward viewer
		double x = (cosO * cosw - sinO * sinw * cosi) * x_orb
				 + (-cosO * sinw - sinO * cosw * cosi) * y_orb;

		double y = (sinO * cosw + cosO * sinw * cosi) * x_orb
				 + (-sinO * sinw + cosO * cosw * cosi) * y_orb;

		double z = (sinw * sini) * x_orb
				 + (cosw * sini) * y_orb;

        // 8. Offset by parent world position
        return {
            parent->position.x + (float)(x * renderScale),
            parent->position.y + (float)(y * renderScale),
            parent->position.z + (float)(z * renderScale)
        };
    }

    // Draw the orbit path as a series of line segments
    void drawOrbit() const {
        if (!parent) return;

        const int segments = 4096;
        Vector3 prev = {0,0,0};

        for (int i = 0; i <= segments; i++) {
            // Step through one full orbit by sweeping mean anomaly
            double M = (2.0 * M_PI / segments) * i;
            double E = solveKepler(M, eccentricity);

            double nu = 2.0 * atan2(
                sqrt(1.0 + eccentricity) * sin(E / 2.0),
                sqrt(1.0 - eccentricity) * cos(E / 2.0)
            );

            double r = semiMajorAxis * (1.0 - eccentricity * cos(E));
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

            Vector3 pt = {
                parent->position.x + (float)(x * renderScale),
                parent->position.y + (float)(y * renderScale),
                parent->position.z + (float)(z * renderScale),
            };

            if (i > 0) DrawLine3D(prev, pt, DARKGRAY);
            prev = pt;
        }
    }
};