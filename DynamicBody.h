#pragma once
#include "Planet.h"
#include <cmath>
#include <string>
#include <algorithm>

class DynamicBody {
public:
    std::string name;
    double      mass;
    double      thrustAccel;
    Color       color = GREEN;

    Planet* parent = nullptr;

    double rx = 0, ry = 0, rz = 0;
    double vx = 0, vy = 0, vz = 0;

    struct Elements {
        double sma   = 0;
        double ecc   = 0;
        double inc   = 0;
        double lan   = 0;
        double aop   = 0;
        double mae   = 0;
        double epoch = 0;
        bool   valid = false;
    } elements;

    bool thrusting = false;
    bool valid     = false;

    double tdx = 0, tdy = 0, tdz = 0;

    DVec3 worldPos;

    int    previewSteps = 500;
    double previewDt    = 60.0;

    static constexpr double G      = 6.674e-11;
    static constexpr double TWO_PI = 2.0 * M_PI;

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

    void initCircularOrbit(Planet* soi, double t0,
                           double altitudeM, double lonDeg = 0.0)
    {
        parent = soi;
        double r  = soi->radius + altitudeM;
        double v  = std::sqrt(G * soi->mass / r);
        double lo = lonDeg * M_PI / 180.0;
        rx = r * std::cos(lo);  ry = r * std::sin(lo);  rz = 0.0;
        vx = -v * std::sin(lo); vy = v * std::cos(lo);  vz = 0.0;
        fitElements(t0);
        valid = true;
    }

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

        // worldPos uses the CURRENT frame's parent position — consistent
        // with what draw() will also see this frame.
        worldPos = {
            parent->position.x + rx,
            parent->position.y + ry,
            parent->position.z + rz
        };
    }

    void onThrustStop(double t) { fitElements(t); }

    void draw(const DVec3& camPos) const {
        if (!valid) return;

        // Snapshot parent position ONCE for this entire draw call.
        // Every piece of geometry — ship body, vectors, orbit, arc — uses
        // this same origin so nothing is offset relative to anything else.
        const DVec3 orig = parent->position;

        // Ship body — derived from orig+state, NOT from worldPos,
        // so it's consistent with the arc start point below.
        DVec3 shipWorld = { orig.x + rx, orig.y + ry, orig.z + rz };
        Vector3 pos = shipWorld.toVec3(camPos);

        DrawSphere(pos, (float)(parent->radius * 0.005), color);
        DrawSphereWires(pos, (float)(parent->radius * 0.005), 8, 8, color);

        // Velocity vector
        double spd = speed();
        if (spd > 0.0) {
            double scale = parent->soiRadius * 0.02 / spd;
            DVec3 velTip = { shipWorld.x + vx*scale,
                             shipWorld.y + vy*scale,
                             shipWorld.z + vz*scale };
            DrawLine3D(pos, velTip.toVec3(camPos), LIME);
        }

        // Thrust vector indicator
        if (thrusting) {
            double tscale = parent->soiRadius * 0.03;
            DVec3 thrustTip = { shipWorld.x + tdx*tscale,
                                shipWorld.y + tdy*tscale,
                                shipWorld.z + tdz*tscale };
            DrawLine3D(pos, thrustTip.toVec3(camPos), ORANGE);
        }

        if (thrusting)
            drawThrustTrajectory(camPos, orig);
        else
            drawOrbit(camPos, orig);
    }

    void fitElements(double t)
    {
        double mu = G * parent->mass;
        double r  = std::sqrt(rx*rx + ry*ry + rz*rz);
        double v2 = vx*vx + vy*vy + vz*vz;

        double hx = ry*vz - rz*vy, hy = rz*vx - rx*vz, hz = rx*vy - ry*vx;
        double h  = std::sqrt(hx*hx + hy*hy + hz*hz);

        double nx = -hy, ny = hx;
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
        double cosE = (ecc + std::cos(nu)) / (1.0 + ecc*std::cos(nu));
        double sinE = std::sqrt(1.0 - ecc*ecc) * std::sin(nu) / (1.0 + ecc*std::cos(nu));
        double mae  = std::atan2(sinE, cosE) - ecc * std::sin(std::atan2(sinE, cosE));

        elements = { sma, ecc, inc, lan, aop, mae, t, true };
    }

    double altitude() const {
        return std::sqrt(rx*rx + ry*ry + rz*rz) - parent->radius;
    }
    double speed() const {
        return std::sqrt(vx*vx + vy*vy + vz*vz);
    }

private:

    // Builds the 3x2 perifocal→inertial rotation matrix
    static void buildRotMatrix(double lan, double inc, double aop,
                               double& Rxx, double& Rxy,
                               double& Ryx, double& Ryy,
                               double& Rzx, double& Rzy)
    {
        double cosO = std::cos(lan), sinO = std::sin(lan);
        double cosi = std::cos(inc), sini = std::sin(inc);
        double cosw = std::cos(aop), sinw = std::sin(aop);
        Rxx =  cosO*cosw - sinO*sinw*cosi;
        Rxy = -cosO*sinw - sinO*cosw*cosi;
        Ryx =  sinO*cosw + cosO*sinw*cosi;
        Ryy = -sinO*sinw + cosO*cosw*cosi;
        Rzx =  sinw*sini;
        Rzy =  cosw*sini;
    }

    void propagateKepler(double t)
    {
        if (!elements.valid) return;
        double mu = G * parent->mass;
        double dt = t - elements.epoch;
        double n  = std::sqrt(mu / std::pow(std::abs(elements.sma), 3.0));
        double M  = std::fmod(elements.mae + n * dt, TWO_PI);
        if (M < 0.0) M += TWO_PI;

        double E  = solveKepler(M, elements.ecc);
        double nu = 2.0 * std::atan2(
            std::sqrt(1.0 + elements.ecc) * std::sin(E / 2.0),
            std::sqrt(1.0 - elements.ecc) * std::cos(E / 2.0));

        double r    = elements.sma * (1.0 - elements.ecc * std::cos(E));
        double xOrb = r * std::cos(nu);
        double yOrb = r * std::sin(nu);

        double p     = elements.sma * (1.0 - elements.ecc * elements.ecc);
        double sqMuP = std::sqrt(mu / p);
        double vxOrb = -sqMuP * std::sin(nu);
        double vyOrb =  sqMuP * (elements.ecc + std::cos(nu));

        double Rxx,Rxy,Ryx,Ryy,Rzx,Rzy;
        buildRotMatrix(elements.lan, elements.inc, elements.aop,
                       Rxx,Rxy,Ryx,Ryy,Rzx,Rzy);

        rx = Rxx*xOrb + Rxy*yOrb;  ry = Ryx*xOrb + Ryy*yOrb;  rz = Rzx*xOrb + Rzy*yOrb;
        vx = Rxx*vxOrb + Rxy*vyOrb; vy = Ryx*vxOrb + Ryy*vyOrb; vz = Rzx*vxOrb + Rzy*vyOrb;
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
            double mu = G * parent->mass;
            double r  = std::sqrt(irx*irx + iry*iry + irz*irz);
            double r3 = r * r * r;
            oax = -mu*irx/r3 + (thrustOn && thrusting ? thrustAccel*tdx : 0.0);
            oay = -mu*iry/r3 + (thrustOn && thrusting ? thrustAccel*tdy : 0.0);
            oaz = -mu*irz/r3 + (thrustOn && thrusting ? thrustAccel*tdz : 0.0);
            odx = ivx; ody = ivy; odz = ivz;
        };

        double ax1,ay1,az1,dx1,dy1,dz1;
        deriv(prx,pry,prz,pvx,pvy,pvz,ax1,ay1,az1,dx1,dy1,dz1);
        double ax2,ay2,az2,dx2,dy2,dz2;
        deriv(prx+dx1*dt*.5,pry+dy1*dt*.5,prz+dz1*dt*.5,
              pvx+ax1*dt*.5,pvy+ay1*dt*.5,pvz+az1*dt*.5,ax2,ay2,az2,dx2,dy2,dz2);
        double ax3,ay3,az3,dx3,dy3,dz3;
        deriv(prx+dx2*dt*.5,pry+dy2*dt*.5,prz+dz2*dt*.5,
              pvx+ax2*dt*.5,pvy+ay2*dt*.5,pvz+az2*dt*.5,ax3,ay3,az3,dx3,dy3,dz3);
        double ax4,ay4,az4,dx4,dy4,dz4;
        deriv(prx+dx3*dt,pry+dy3*dt,prz+dz3*dt,
              pvx+ax3*dt,pvy+ay3*dt,pvz+az3*dt,ax4,ay4,az4,dx4,dy4,dz4);

        prx += (dt/6.0)*(dx1+2*dx2+2*dx3+dx4);
        pry += (dt/6.0)*(dy1+2*dy2+2*dy3+dy4);
        prz += (dt/6.0)*(dz1+2*dz2+2*dz3+dz4);
        pvx += (dt/6.0)*(ax1+2*ax2+2*ax3+ax4);
        pvy += (dt/6.0)*(ay1+2*ay2+2*ay3+ay4);
        pvz += (dt/6.0)*(az1+2*az2+2*az3+az4);
    }

    // Draws a full Keplerian ellipse given elements, using the provided origin
    void drawEllipse(const DVec3& orig, const DVec3& camPos,
                     double sma, double ecc,
                     double lan, double inc, double aop,
                     Color col) const
    {
        double Rxx,Rxy,Ryx,Ryy,Rzx,Rzy;
        buildRotMatrix(lan, inc, aop, Rxx,Rxy,Ryx,Ryy,Rzx,Rzy);

        const int segments = 256;
        Vector3 prev = {};
        for (int i = 0; i <= segments; i++) {
            double M  = TWO_PI * i / segments;
            double E  = solveKepler(M, ecc);
            double nu = 2.0 * std::atan2(
                std::sqrt(1.0 + ecc) * std::sin(E / 2.0),
                std::sqrt(1.0 - ecc) * std::cos(E / 2.0));
            double r    = sma * (1.0 - ecc * std::cos(E));
            double xOrb = r * std::cos(nu);
            double yOrb = r * std::sin(nu);
            DVec3 wp = {
                orig.x + Rxx*xOrb + Rxy*yOrb,
                orig.y + Ryx*xOrb + Ryy*yOrb,
                orig.z + Rzx*xOrb + Rzy*yOrb
            };
            Vector3 pt = wp.toVec3(camPos);
            if (i > 0) DrawLine3D(prev, pt, col);
            prev = pt;
        }

        // Pe marker (nu=0 → xOrb=rPe, yOrb=0)
        double rPe = sma * (1.0 - ecc);
        DVec3 peri = { orig.x + Rxx*rPe, orig.y + Ryx*rPe, orig.z + Rzx*rPe };
        DrawSphere(peri.toVec3(camPos), (float)(parent->radius * 0.003), YELLOW);

        // Ap marker (nu=π → xOrb=-rAp, yOrb=0)
        double rAp = sma * (1.0 + ecc);
        DVec3 apo  = { orig.x - Rxx*rAp, orig.y - Ryx*rAp, orig.z - Rzx*rAp };
        DrawSphere(apo.toVec3(camPos), (float)(parent->radius * 0.003), SKYBLUE);
    }

    // Coasting orbit — uses snapshotted orig passed from draw()
    void drawOrbit(const DVec3& camPos, const DVec3& orig) const {
        if (!elements.valid || elements.ecc >= 1.0) return;
        drawEllipse(orig, camPos,
                    elements.sma, elements.ecc,
                    elements.lan, elements.inc, elements.aop,
                    { 0, 200, 100, 180 });
    }

    // Thrust trajectory preview — uses the same snapshotted orig from draw()
    void drawThrustTrajectory(const DVec3& camPos, const DVec3& orig) const {
        if (!valid || !parent) return;

        // ---- Phase 1: burn arc ----
        // Start from current local state — ship position is orig+rx/ry/rz
        double prx = rx, pry = ry, prz = rz;
        double pvx = vx, pvy = vy, pvz = vz;

        // First point: ship's current position, derived from same orig
        Vector3 prev = DVec3{ orig.x+rx, orig.y+ry, orig.z+rz }.toVec3(camPos);
        bool hitGround = false;

        for (int i = 0; i < previewSteps; i++) {
            integrateRK4(previewDt, prx, pry, prz, pvx, pvy, pvz, true);

            if (std::sqrt(prx*prx + pry*pry + prz*prz) < parent->radius) {
                hitGround = true;
                break;
            }

            float  t   = (float)i / (float)(previewSteps - 1);
            Color  col = {
                255,
                (unsigned char)(165 + (int)(90.0f * t)),
                0,
                (unsigned char)(220 - (int)(80.0f * t))
            };

            // Arc point: orig + local position — same orig as ship body
            Vector3 pt = DVec3{ orig.x+prx, orig.y+pry, orig.z+prz }.toVec3(camPos);
            DrawLine3D(prev, pt, col);
            prev = pt;
        }

        if (hitGround) return;

        // ---- Fit elements from burn-end state ----
        double mu = G * parent->mass;
        double r  = std::sqrt(prx*prx + pry*pry + prz*prz);
        double v2 = pvx*pvx + pvy*pvy + pvz*pvz;

        double hx = pry*pvz - prz*pvy, hy = prz*pvx - prx*pvz, hz = prx*pvy - pry*pvx;
        double h  = std::sqrt(hx*hx + hy*hy + hz*hz);
        double nx = -hy, ny = hx;
        double nMag = std::sqrt(nx*nx + ny*ny);

        double evx = (pvy*hz - pvz*hy)/mu - prx/r;
        double evy = (pvz*hx - pvx*hz)/mu - pry/r;
        double evz = (pvx*hy - pvy*hx)/mu - prz/r;
        double ecc = std::sqrt(evx*evx + evy*evy + evz*evz);
        double sma = 1.0 / (2.0/r - v2/mu);

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

        // ---- Phase 2: coast ellipse — SAME orig as burn arc ----
        drawEllipse(orig, camPos, sma, ecc, lan, inc, aop, { 0, 220, 120, 160 });
    }
};