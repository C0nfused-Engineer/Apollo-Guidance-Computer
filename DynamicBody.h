#pragma once
#include "Planet.h"
#include <cmath>
#include <string>
#include <algorithm>

// ---------------------------------------------------------------------------
// DynamicBody — a spacecraft or free-flying object that:
//   - Coasts analytically on osculating Keplerian elements (zero drift)
//   - Burns with RK4 Newtonian integration when thrusting
//   - Refits osculating elements from Cartesian state on burn-end
//
// State vectors (rx,ry,rz,vx,vy,vz) are in the parent body's
// inertial frame, SI units (meters, m/s).
//
// Coordinate convention matches Planet.h:
//   X = right, Y = forward, Z = ecliptic north (up)
// ---------------------------------------------------------------------------

class DynamicBody {
public:
    std::string name;
    double      mass;          // kg
    double      thrustAccel;   // m/s² (specific force = F/m)
    Color       color = GREEN;

    Planet* parent = nullptr;  // current SOI body

    // Cartesian state in parent-centered inertial frame (meters, m/s)
    double rx = 0, ry = 0, rz = 0;
    double vx = 0, vy = 0, vz = 0;

    // Osculating Keplerian elements — refit on burn-end, used for coasting
    struct Elements {
        double sma   = 0;  // semi-major axis (m)
        double ecc   = 0;  // eccentricity
        double inc   = 0;  // inclination (rad)
        double lan   = 0;  // longitude of ascending node (rad)
        double aop   = 0;  // argument of periapsis (rad)
        double mae   = 0;  // mean anomaly at epoch (rad)
        double epoch = 0;  // sim time when mae was measured (s)
        bool   valid = false;
    } elements;

    bool thrusting  = false;
    bool valid      = false;   // false until initialised

    // Thrust direction unit vector — set each frame by the controller
    // expressed in parent-centered inertial frame (same as state vectors)
    double tdx = 0, tdy = 0, tdz = 0;

    // World-space position (DVec3, meters) — set by update(), used for draw/camera
    DVec3 worldPos;

    static constexpr double G      = 6.674e-11;
    static constexpr double TWO_PI = 2.0 * M_PI;

    // -----------------------------------------------------------------------
    // Initialise from an explicit Cartesian state at sim time t0
    // -----------------------------------------------------------------------
    void initFromState(Planet* soi, double t0,
                       double px, double py, double pz,
                       double pvx, double pvy, double pvz)
    {
        parent = soi;
        rx = px;  ry = py;  rz = pz;
        vx = pvx; vy = pvy; vz = pvz;
        fitElements(t0);
        valid = true;
    }

    // -----------------------------------------------------------------------
    // Initialise in a prograde circular orbit at given altitude above parent
    // surface, starting at longitude lonDeg in the ecliptic (XY) plane.
    // -----------------------------------------------------------------------
    void initCircularOrbit(Planet* soi, double t0,
                           double altitudeM, double lonDeg = 0.0)
    {
        parent = soi;
        double r  = soi->radius + altitudeM;
        double v  = std::sqrt(G * soi->mass / r);
        double lo = lonDeg * M_PI / 180.0;

        // Position in XY plane
        rx = r * std::cos(lo);
        ry = r * std::sin(lo);
        rz = 0.0;

        // Velocity perpendicular to radius (prograde), in XY plane
        vx = -v * std::sin(lo);
        vy =  v * std::cos(lo);
        vz = 0.0;

        fitElements(t0);
        valid = true;
    }

    // -----------------------------------------------------------------------
    // Main update — call once per frame
    //   t  = current sim time (s)
    //   dt = sim-time step this frame (s)
    // -----------------------------------------------------------------------
    void update(double t, double dt)
    {
        if (!valid || !parent) return;

        if (thrusting) {
            // Sub-step RK4 — cap step size at 10 s for accuracy
            const double maxStep = 10.0;
            double remaining = dt;
            while (remaining > 0.0) {
                double step = std::min(remaining, maxStep);
                integrateRK4(step);
                remaining -= step;
            }
        } else {
            // Analytic Keplerian propagation — no drift
            propagateKepler(t);
        }

        // Update world position (parent position + local offset)
        worldPos = {
            parent->position.x + rx,
            parent->position.y + ry,
            parent->position.z + rz
        };
    }

    // Call when thrust stops to re-sync elements from current state
    void onThrustStop(double t) {
        fitElements(t);
    }

    // -----------------------------------------------------------------------
    // Draw — call between BeginMode3D / EndMode3D
    // camPos: world-space camera position (meters), same as Planet::draw
    // -----------------------------------------------------------------------
    void draw(const DVec3& camPos) const {
        if (!valid) return;

        Vector3 pos = worldPos.toVec3(camPos);

        // Body — small bright sphere
        DrawSphere(pos, (float)(parent->radius * 0.005), color);
        DrawSphereWires(pos, (float)(parent->radius * 0.005), 8, 8, color);

        // Velocity vector (scaled for visibility)
        double speed = std::sqrt(vx*vx + vy*vy + vz*vz);
        if (speed > 0.0) {
            double scale = parent->soiRadius * 0.02 / speed;
            DVec3 velTip = {
                worldPos.x + vx * scale,
                worldPos.y + vy * scale,
                worldPos.z + vz * scale
            };
            DrawLine3D(pos, velTip.toVec3(camPos), LIME);
        }

        // Thrust vector
        if (thrusting) {
            double tscale = parent->soiRadius * 0.03;
            DVec3 thrustTip = {
                worldPos.x + tdx * tscale,
                worldPos.y + tdy * tscale,
                worldPos.z + tdz * tscale
            };
            DrawLine3D(pos, thrustTip.toVec3(camPos), ORANGE);
        }

        // Predicted orbit path (current osculating elements)
        drawOrbit(camPos);
    }

    // -----------------------------------------------------------------------
    // Refit osculating Keplerian elements from current Cartesian state.
    // Uses angular momentum + eccentricity vector approach (Bate, Mueller, White).
    // -----------------------------------------------------------------------
    void fitElements(double t)
    {
        double mu = G * parent->mass;
        double r  = std::sqrt(rx*rx + ry*ry + rz*rz);
        double v2 = vx*vx + vy*vy + vz*vz;

        // Angular momentum h = r × v
        double hx = ry*vz - rz*vy;
        double hy = rz*vx - rx*vz;
        double hz = rx*vy - ry*vx;
        double h  = std::sqrt(hx*hx + hy*hy + hz*hz);

        // Node vector n = ẑ × h  (ẑ = ecliptic north = Z-axis)
        double nx   = -hy;
        double ny   =  hx;
        double nMag = std::sqrt(nx*nx + ny*ny);

        // Eccentricity vector e = (v × h)/μ − r̂
        double evx = (vy*hz - vz*hy)/mu - rx/r;
        double evy = (vz*hx - vx*hz)/mu - ry/r;
        double evz = (vx*hy - vy*hx)/mu - rz/r;
        double ecc = std::sqrt(evx*evx + evy*evy + evz*evz);

        // Semi-major axis (vis-viva: ε = v²/2 − μ/r = −μ/2a)
        double sma = 1.0 / (2.0/r - v2/mu);

        // Inclination: cos i = hz / |h|
        double inc = std::acos(std::clamp(hz / h, -1.0, 1.0));

        // Longitude of ascending node
        double lan = 0.0;
        if (nMag > 1e-10) {
            lan = std::acos(std::clamp(nx / nMag, -1.0, 1.0));
            if (ny < 0.0) lan = TWO_PI - lan;
        }

        // Argument of periapsis
        double aop = 0.0;
        if (nMag > 1e-10 && ecc > 1e-10) {
            double ndote = (nx*evx + ny*evy) / (nMag * ecc);
            aop = std::acos(std::clamp(ndote, -1.0, 1.0));
            if (evz < 0.0) aop = TWO_PI - aop;
        }

        // True anomaly
        double nu = 0.0;
        if (ecc > 1e-10) {
            double rdote = (rx*evx + ry*evy + rz*evz) / (r * ecc);
            nu = std::acos(std::clamp(rdote, -1.0, 1.0));
            double rdotv = rx*vx + ry*vy + rz*vz;
            if (rdotv < 0.0) nu = TWO_PI - nu;
        }

        // True → eccentric → mean anomaly
        double cosE = (ecc + std::cos(nu)) / (1.0 + ecc * std::cos(nu));
        double sinE = std::sqrt(1.0 - ecc*ecc) * std::sin(nu) / (1.0 + ecc * std::cos(nu));
        double E    = std::atan2(sinE, cosE);
        double mae  = E - ecc * std::sin(E);

        elements = { sma, ecc, inc, lan, aop, mae, t, true };
    }

    // Altitude above parent surface (meters)
    double altitude() const {
        return std::sqrt(rx*rx + ry*ry + rz*rz) - parent->radius;
    }

    // Speed (m/s)
    double speed() const {
        return std::sqrt(vx*vx + vy*vy + vz*vz);
    }

private:
    // -----------------------------------------------------------------------
    // Analytic Keplerian propagation — solves for position AND velocity at t
    // -----------------------------------------------------------------------
    void propagateKepler(double t)
    {
        if (!elements.valid) return;

        double mu  = G * parent->mass;
        double dt  = t - elements.epoch;
        double n   = std::sqrt(mu / std::pow(std::abs(elements.sma), 3.0));
        double M   = std::fmod(elements.mae + n * dt, TWO_PI);
        if (M < 0.0) M += TWO_PI;

        double E = solveKepler(M, elements.ecc);

        double nu = 2.0 * std::atan2(
            std::sqrt(1.0 + elements.ecc) * std::sin(E / 2.0),
            std::sqrt(1.0 - elements.ecc) * std::cos(E / 2.0)
        );

        double r    = elements.sma * (1.0 - elements.ecc * std::cos(E));
        double xOrb = r * std::cos(nu);
        double yOrb = r * std::sin(nu);

        // Perifocal velocity
        double p      = elements.sma * (1.0 - elements.ecc * elements.ecc);
        double sqMuP  = std::sqrt(mu / p);
        double vxOrb  = -sqMuP * std::sin(nu);
        double vyOrb  =  sqMuP * (elements.ecc + std::cos(nu));

        // Perifocal → inertial (matches Planet.h convention: XY ecliptic, Z north)
        double cosO = std::cos(elements.lan), sinO = std::sin(elements.lan);
        double cosi = std::cos(elements.inc), sini = std::sin(elements.inc);
        double cosw = std::cos(elements.aop), sinw = std::sin(elements.aop);

        double Rxx =  cosO*cosw - sinO*sinw*cosi;
        double Rxy = -cosO*sinw - sinO*cosw*cosi;
        double Ryx =  sinO*cosw + cosO*sinw*cosi;
        double Ryy = -sinO*sinw + cosO*cosw*cosi;
        double Rzx =  sinw*sini;
        double Rzy =  cosw*sini;

        rx = Rxx*xOrb + Rxy*yOrb;
        ry = Ryx*xOrb + Ryy*yOrb;
        rz = Rzx*xOrb + Rzy*yOrb;

        vx = Rxx*vxOrb + Rxy*vyOrb;
        vy = Ryx*vxOrb + Ryy*vyOrb;
        vz = Rzx*vxOrb + Rzy*vyOrb;
    }

    // Newton-Raphson Kepler solver
    double solveKepler(double M, double e) const {
        double E = (e < 0.8) ? M : M_PI;
        for (int i = 0; i < 100; i++) {
            double dE = (M - E + e * std::sin(E)) / (1.0 - e * std::cos(E));
            E += dE;
            if (std::fabs(dE) < 1e-12) break;
        }
        return E;
    }

    // -----------------------------------------------------------------------
    // RK4 integration — equations of motion: gravity + thrust
    // -----------------------------------------------------------------------
    struct State6 { double rx,ry,rz,vx,vy,vz; };

    State6 deriv(const State6& s) const {
        double mu = G * parent->mass;
        double r  = std::sqrt(s.rx*s.rx + s.ry*s.ry + s.rz*s.rz);
        double r3 = r * r * r;

        double ax = -mu * s.rx / r3 + (thrusting ? thrustAccel * tdx : 0.0);
        double ay = -mu * s.ry / r3 + (thrusting ? thrustAccel * tdy : 0.0);
        double az = -mu * s.rz / r3 + (thrusting ? thrustAccel * tdz : 0.0);

        return { s.vx, s.vy, s.vz, ax, ay, az };
    }

    static State6 add(const State6& a, const State6& b, double h) {
        return {
            a.rx + h*b.rx, a.ry + h*b.ry, a.rz + h*b.rz,
            a.vx + h*b.vx, a.vy + h*b.vy, a.vz + h*b.vz
        };
    }

    void integrateRK4(double dt)
    {
        State6 s = { rx, ry, rz, vx, vy, vz };
        State6 k1 = deriv(s);
        State6 k2 = deriv(add(s, k1, dt*0.5));
        State6 k3 = deriv(add(s, k2, dt*0.5));
        State6 k4 = deriv(add(s, k3, dt));

        rx = s.rx + (dt/6.0)*(k1.rx + 2*k2.rx + 2*k3.rx + k4.rx);
        ry = s.ry + (dt/6.0)*(k1.ry + 2*k2.ry + 2*k3.ry + k4.ry);
        rz = s.rz + (dt/6.0)*(k1.rz + 2*k2.rz + 2*k3.rz + k4.rz);
        vx = s.vx + (dt/6.0)*(k1.vx + 2*k2.vx + 2*k3.vx + k4.vx);
        vy = s.vy + (dt/6.0)*(k1.vy + 2*k2.vy + 2*k3.vy + k4.vy);
        vz = s.vz + (dt/6.0)*(k1.vz + 2*k2.vz + 2*k3.vz + k4.vz);
    }

    // -----------------------------------------------------------------------
    // Draw the predicted osculating orbit ellipse
    // -----------------------------------------------------------------------
    void drawOrbit(const DVec3& camPos) const {
        if (!elements.valid || elements.ecc >= 1.0) return;

        const int segments = 256;
        Vector3 prev = {};

        double cosO = std::cos(elements.lan), sinO = std::sin(elements.lan);
        double cosi = std::cos(elements.inc), sini = std::sin(elements.inc);
        double cosw = std::cos(elements.aop), sinw = std::sin(elements.aop);

        double Rxx =  cosO*cosw - sinO*sinw*cosi;
        double Rxy = -cosO*sinw - sinO*cosw*cosi;
        double Ryx =  sinO*cosw + cosO*sinw*cosi;
        double Ryy = -sinO*sinw + cosO*cosw*cosi;
        double Rzx =  sinw*sini;
        double Rzy =  cosw*sini;

        for (int i = 0; i <= segments; i++) {
            double M  = TWO_PI * i / segments;
            double E  = solveKepler(M, elements.ecc);
            double nu = 2.0 * std::atan2(
                std::sqrt(1.0 + elements.ecc) * std::sin(E / 2.0),
                std::sqrt(1.0 - elements.ecc) * std::cos(E / 2.0)
            );
            double r    = elements.sma * (1.0 - elements.ecc * std::cos(E));
            double xOrb = r * std::cos(nu);
            double yOrb = r * std::sin(nu);

            DVec3 wp = {
                parent->position.x + Rxx*xOrb + Rxy*yOrb,
                parent->position.y + Ryx*xOrb + Ryy*yOrb,
                parent->position.z + Rzx*xOrb + Rzy*yOrb
            };

            Vector3 pt = wp.toVec3(camPos);
            if (i > 0) DrawLine3D(prev, pt, { 0, 200, 100, 180 });
            prev = pt;
        }
    }
};