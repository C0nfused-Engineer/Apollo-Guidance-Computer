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
//   - During thrust: draws a live forward-integrated trajectory preview
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
    bool valid      = false;

    // Thrust direction unit vector in parent-centered inertial frame
    double tdx = 0, tdy = 0, tdz = 0;

    // World-space position (meters) — set by update()
    DVec3 worldPos;

    // ---------------------------------------------------------------------------
    // Trajectory preview during thrust.
    // previewSteps  : how many RK4 steps to project forward
    // previewDt     : sim-seconds per step (e.g. 60 = one step per minute)
    // These can be tuned from main without recompiling.
    // ---------------------------------------------------------------------------
    int    previewSteps = 500;
    double previewDt    = 60.0;   // seconds per step

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

        rx = r * std::cos(lo);
        ry = r * std::sin(lo);
        rz = 0.0;

        vx = -v * std::sin(lo);
        vy =  v * std::cos(lo);
        vz = 0.0;

        fitElements(t0);
        valid = true;
    }

    // -----------------------------------------------------------------------
    // Main update — call once per frame
    // -----------------------------------------------------------------------
    void update(double t, double dt)
    {
        if (!valid || !parent) return;

        if (thrusting) {
            const double maxStep = 10.0;
            double remaining = dt;
            while (remaining > 0.0) {
                double step = std::min(remaining, maxStep);
                integrateRK4(step, rx, ry, rz, vx, vy, vz);
                remaining -= step;
            }
        } else {
            propagateKepler(t);
        }

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
    // -----------------------------------------------------------------------
    void draw(const DVec3& camPos) const {
        if (!valid) return;

        Vector3 pos = worldPos.toVec3(camPos);

        // Body
        DrawSphere(pos, (float)(parent->radius * 0.005), color);
        DrawSphereWires(pos, (float)(parent->radius * 0.005), 8, 8, color);

        // Velocity vector
        double spd = speed();
        if (spd > 0.0) {
            double scale = parent->soiRadius * 0.02 / spd;
            DVec3 velTip = {
                worldPos.x + vx * scale,
                worldPos.y + vy * scale,
                worldPos.z + vz * scale
            };
            DrawLine3D(pos, velTip.toVec3(camPos), LIME);
        }

        // Thrust vector indicator
        if (thrusting) {
            double tscale = parent->soiRadius * 0.03;
            DVec3 thrustTip = {
                worldPos.x + tdx * tscale,
                worldPos.y + tdy * tscale,
                worldPos.z + tdz * tscale
            };
            DrawLine3D(pos, thrustTip.toVec3(camPos), ORANGE);
        }

        if (thrusting)
            drawThrustTrajectory(camPos);
        else
            drawOrbit(camPos);
    }

    // -----------------------------------------------------------------------
    // Refit osculating Keplerian elements from current Cartesian state
    // -----------------------------------------------------------------------
    void fitElements(double t)
    {
        double mu = G * parent->mass;
        double r  = std::sqrt(rx*rx + ry*ry + rz*rz);
        double v2 = vx*vx + vy*vy + vz*vz;

        double hx = ry*vz - rz*vy;
        double hy = rz*vx - rx*vz;
        double hz = rx*vy - ry*vx;
        double h  = std::sqrt(hx*hx + hy*hy + hz*hz);

        double nx   = -hy;
        double ny   =  hx;
        double nMag = std::sqrt(nx*nx + ny*ny);

        double evx = (vy*hz - vz*hy)/mu - rx/r;
        double evy = (vz*hx - vx*hz)/mu - ry/r;
        double evz = (vx*hy - vy*hx)/mu - rz/r;
        double ecc = std::sqrt(evx*evx + evy*evy + evz*evz);

        double sma = 1.0 / (2.0/r - v2/mu);
        double inc = std::acos(std::clamp(hz / h, -1.0, 1.0));

        double lan = 0.0;
        if (nMag > 1e-10) {
            lan = std::acos(std::clamp(nx / nMag, -1.0, 1.0));
            if (ny < 0.0) lan = TWO_PI - lan;
        }

        double aop = 0.0;
        if (nMag > 1e-10 && ecc > 1e-10) {
            double ndote = (nx*evx + ny*evy) / (nMag * ecc);
            aop = std::acos(std::clamp(ndote, -1.0, 1.0));
            if (evz < 0.0) aop = TWO_PI - aop;
        }

        double nu = 0.0;
        if (ecc > 1e-10) {
            double rdote = (rx*evx + ry*evy + rz*evz) / (r * ecc);
            nu = std::acos(std::clamp(rdote, -1.0, 1.0));
            if (rx*vx + ry*vy + rz*vz < 0.0) nu = TWO_PI - nu;
        }

        double cosE = (ecc + std::cos(nu)) / (1.0 + ecc * std::cos(nu));
        double sinE = std::sqrt(1.0 - ecc*ecc) * std::sin(nu) / (1.0 + ecc * std::cos(nu));
        double E    = std::atan2(sinE, cosE);
        double mae  = E - ecc * std::sin(E);

        elements = { sma, ecc, inc, lan, aop, mae, t, true };
    }

    double altitude() const {
        return std::sqrt(rx*rx + ry*ry + rz*rz) - parent->radius;
    }

    double speed() const {
        return std::sqrt(vx*vx + vy*vy + vz*vz);
    }

private:
    // -----------------------------------------------------------------------
    // Analytic Keplerian propagation
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

        double p     = elements.sma * (1.0 - elements.ecc * elements.ecc);
        double sqMuP = std::sqrt(mu / p);
        double vxOrb = -sqMuP * std::sin(nu);
        double vyOrb =  sqMuP * (elements.ecc + std::cos(nu));

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
    // RK4 step — operates on explicit state passed by reference so we can
    // use it both for live integration AND for trajectory preview without
    // touching the real state.
    // thrustOn: whether to include thrust acceleration in this integration
    // -----------------------------------------------------------------------
    void integrateRK4(double dt,
                      double& prx, double& pry, double& prz,
                      double& pvx, double& pvy, double& pvz,
                      bool thrustOn = true) const
    {
        auto deriv = [&](double irx, double iry, double irz,
                         double ivx, double ivy, double ivz,
                         double& oax, double& oay, double& oaz,
                         double& odx, double& ody, double& odz)
        {
            double mu  = G * parent->mass;
            double r   = std::sqrt(irx*irx + iry*iry + irz*irz);
            double r3  = r * r * r;
            oax = -mu * irx / r3 + (thrustOn && thrusting ? thrustAccel * tdx : 0.0);
            oay = -mu * iry / r3 + (thrustOn && thrusting ? thrustAccel * tdy : 0.0);
            oaz = -mu * irz / r3 + (thrustOn && thrusting ? thrustAccel * tdz : 0.0);
            odx = ivx; ody = ivy; odz = ivz;
        };

        double ax1,ay1,az1, dx1,dy1,dz1;
        deriv(prx,pry,prz, pvx,pvy,pvz, ax1,ay1,az1, dx1,dy1,dz1);

        double ax2,ay2,az2, dx2,dy2,dz2;
        deriv(prx+dx1*dt*.5, pry+dy1*dt*.5, prz+dz1*dt*.5,
              pvx+ax1*dt*.5, pvy+ay1*dt*.5, pvz+az1*dt*.5,
              ax2,ay2,az2, dx2,dy2,dz2);

        double ax3,ay3,az3, dx3,dy3,dz3;
        deriv(prx+dx2*dt*.5, pry+dy2*dt*.5, prz+dz2*dt*.5,
              pvx+ax2*dt*.5, pvy+ay2*dt*.5, pvz+az2*dt*.5,
              ax3,ay3,az3, dx3,dy3,dz3);

        double ax4,ay4,az4, dx4,dy4,dz4;
        deriv(prx+dx3*dt, pry+dy3*dt, prz+dz3*dt,
              pvx+ax3*dt, pvy+ay3*dt, pvz+az3*dt,
              ax4,ay4,az4, dx4,dy4,dz4);

        prx += (dt/6.0)*(dx1+2*dx2+2*dx3+dx4);
        pry += (dt/6.0)*(dy1+2*dy2+2*dy3+dy4);
        prz += (dt/6.0)*(dz1+2*dz2+2*dz3+dz4);
        pvx += (dt/6.0)*(ax1+2*ax2+2*ax3+ax4);
        pvy += (dt/6.0)*(ay1+2*ay2+2*ay3+ay4);
        pvz += (dt/6.0)*(az1+2*az2+2*az3+az4);
    }

    // -----------------------------------------------------------------------
    // Draw the predicted osculating orbit ellipse (coasting mode)
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

    // -----------------------------------------------------------------------
    // Draw a live forward-integrated trajectory during thrust (KSP-style).
    //
    // Phase 1 — BURN segment: integrate with thrust ON for previewSteps steps.
    //   Drawn in orange, fading toward yellow as it extends forward.
    //
    // Phase 2 — COAST segment: from the burn-end state, fit osculating elements
    //   and draw one full Keplerian orbit in green, so you can see exactly
    //   where the resulting orbit ends up.
    // -----------------------------------------------------------------------
    void drawThrustTrajectory(const DVec3& camPos) const {
        if (!valid || !parent) return;

        // --- Phase 1: burn arc ---
        double prx = rx, pry = ry, prz = rz;
        double pvx = vx, pvy = vy, pvz = vz;

        Vector3 prev = worldPos.toVec3(camPos);
        bool hitGround = false;

        // Colour gradient: orange → yellow over the burn arc
        for (int i = 0; i < previewSteps; i++) {
            integrateRK4(previewDt, prx, pry, prz, pvx, pvy, pvz, true);

            // Stop drawing if we hit the parent body
            double dist = std::sqrt(prx*prx + pry*pry + prz*prz);
            if (dist < parent->radius) { hitGround = true; break; }

            float t   = (float)i / (float)(previewSteps - 1);  // 0→1
            Color col = {
                255,
                (unsigned char)(165 + (int)(90 * t)),   // 165→255 (orange→yellow)
                0,
                (unsigned char)(220 - (int)(80 * t))    // fade alpha slightly
            };

            DVec3 wp = {
                parent->position.x + prx,
                parent->position.y + pry,
                parent->position.z + prz
            };
            Vector3 pt = wp.toVec3(camPos);
            DrawLine3D(prev, pt, col);
            prev = pt;
        }

        if (hitGround) return;  // impactor — no coast phase

        // --- Phase 2: coast arc (osculating orbit from burn-end state) ---
        // Fit elements from the projected burn-end state
        double mu  = G * parent->mass;
        double r   = std::sqrt(prx*prx + pry*pry + prz*prz);
        double v2  = pvx*pvx + pvy*pvy + pvz*pvz;

        double hx  = pry*pvz - prz*pvy;
        double hy  = prz*pvx - prx*pvz;
        double hz  = prx*pvy - pry*pvx;
        double h   = std::sqrt(hx*hx + hy*hy + hz*hz);

        double nx   = -hy, ny = hx;
        double nMag = std::sqrt(nx*nx + ny*ny);

        double evx = (pvy*hz - pvz*hy)/mu - prx/r;
        double evy = (pvz*hx - pvx*hz)/mu - pry/r;
        double evz = (pvx*hy - pvy*hx)/mu - prz/r;
        double ecc = std::sqrt(evx*evx + evy*evy + evz*evz);

        double sma = 1.0 / (2.0/r - v2/mu);

        // Skip hyperbolic or nearly-degenerate orbits
        if (ecc >= 1.0 || sma <= 0.0) return;

        double inc = std::acos(std::clamp(hz / h, -1.0, 1.0));

        double lan = 0.0;
        if (nMag > 1e-10) {
            lan = std::acos(std::clamp(nx / nMag, -1.0, 1.0));
            if (ny < 0.0) lan = TWO_PI - lan;
        }

        double aop = 0.0;
        if (nMag > 1e-10 && ecc > 1e-10) {
            double ndote = (nx*evx + ny*evy) / (nMag * ecc);
            aop = std::acos(std::clamp(ndote, -1.0, 1.0));
            if (evz < 0.0) aop = TWO_PI - aop;
        }

        // Build rotation matrix for the projected orbit
        double cosO = std::cos(lan),  sinO = std::sin(lan);
        double cosi = std::cos(inc),  sini = std::sin(inc);
        double cosw = std::cos(aop),  sinw = std::sin(aop);

        double Rxx =  cosO*cosw - sinO*sinw*cosi;
        double Rxy = -cosO*sinw - sinO*cosw*cosi;
        double Ryx =  sinO*cosw + cosO*sinw*cosi;
        double Ryy = -sinO*sinw + cosO*cosw*cosi;
        double Rzx =  sinw*sini;
        double Rzy =  cosw*sini;

        // Draw one full orbit of the projected coast ellipse
        const int segments = 256;
        Vector3 cprev = {};
        for (int i = 0; i <= segments; i++) {
            double M  = TWO_PI * i / segments;
            double E  = solveKepler(M, ecc);
            double nu = 2.0 * std::atan2(
                std::sqrt(1.0 + ecc) * std::sin(E / 2.0),
                std::sqrt(1.0 - ecc) * std::cos(E / 2.0)
            );
            double cr   = sma * (1.0 - ecc * std::cos(E));
            double xOrb = cr * std::cos(nu);
            double yOrb = cr * std::sin(nu);

            DVec3 wp = {
                parent->position.x + Rxx*xOrb + Rxy*yOrb,
                parent->position.y + Ryx*xOrb + Ryy*yOrb,
                parent->position.z + Rzx*xOrb + Rzy*yOrb
            };
            Vector3 pt = wp.toVec3(camPos);
            if (i > 0) DrawLine3D(cprev, pt, { 0, 220, 120, 160 });
            cprev = pt;
        }

        // Periapsis marker — lowest point of projected orbit
        {
            double E_pe  = 0.0;  // E=0 → periapsis
            double nu_pe = 0.0;
            double rpe   = sma * (1.0 - ecc);
            DVec3 peri = {
                parent->position.x + Rxx*rpe + Rxy*0,
                parent->position.y + Ryx*rpe + Ryy*0,
                parent->position.z + Rzx*rpe + Rzy*0
            };
            (void)E_pe; (void)nu_pe;
            DrawSphere(peri.toVec3(camPos), (float)(parent->radius * 0.003), YELLOW);
        }

        // Apoapsis marker
        {
            double rAp   = sma * (1.0 + ecc);
            DVec3 apo = {
                parent->position.x + Rxx*(-rAp) + Rxy*0,
                parent->position.y + Ryx*(-rAp) + Ryy*0,
                parent->position.z + Rzx*(-rAp) + Rzy*0
            };
            DrawSphere(apo.toVec3(camPos), (float)(parent->radius * 0.003), SKYBLUE);
        }
    }
};