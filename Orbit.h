#pragma once
#include "raylib.h"
#include "raymath.h"
#include <cmath>
#include <vector>

// Keplerian orbital elements + drawing for any conic section.
// Works for elliptic (ecc < 1), parabolic (ecc == 1) and hyperbolic (ecc > 1) orbits.
struct Orbit {
    double sma  = 0.0;   // semi-major axis (meters); negative = hyperbolic
    double ecc  = 0.0;   // eccentricity
    double inc  = 0.0;   // inclination (radians)
    double lan  = 0.0;   // longitude of ascending node (radians)
    double aop  = 0.0;   // argument of periapsis (radians)
    double mae  = 0.0;   // mean anomaly at epoch (radians)
    double epoch = 0.0;  // epoch time (seconds)
    bool   valid = false;

    static constexpr double G      = 6.674e-11;
    static constexpr double TWO_PI = 2.0 * M_PI;

    bool isElliptic()    const { return valid && ecc <  1.0 && sma > 0.0; }
    bool isHyperbolic()  const { return valid && ecc >= 1.0; }

    // -----------------------------------------------------------------------
    // Build rotation matrix columns from orbital elements.
    // x_world = Rxx*x_orb + Rxy*y_orb   etc.
    // -----------------------------------------------------------------------
    void buildRotMatrix(double& Rxx, double& Rxy,
                        double& Ryx, double& Ryy,
                        double& Rzx, double& Rzy) const
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

    // -----------------------------------------------------------------------
    // Solve Kepler's equation for elliptic orbits (Newton iteration).
    // -----------------------------------------------------------------------
    double solveKepler(double M, double e) const {
        M = std::fmod(M, TWO_PI);
        if (M < 0.0) M += TWO_PI;
        double E = (e < 0.8) ? M : M_PI;
        for (int i = 0; i < 100; i++) {
            double dE = (M - E + e * std::sin(E)) / (1.0 - e * std::cos(E));
            E += dE;
            if (std::fabs(dE) < 1e-12) break;
        }
        return E;
    }

    // -----------------------------------------------------------------------
    // Solve the hyperbolic Kepler equation: M = e*sinh(H) - H
    // -----------------------------------------------------------------------
    double solveKeplerHyperbolic(double M, double e) const {
        double H = (std::fabs(M) < 1.0) ? M : std::copysign(std::log(2.0*std::fabs(M)/e + 1.8), M);
        for (int i = 0; i < 100; i++) {
            double dH = (M - e*std::sinh(H) + H) / (e*std::cosh(H) - 1.0);
            H += dH;
            if (std::fabs(dH) < 1e-12) break;
        }
        return H;
    }

    // -----------------------------------------------------------------------
    // Position on orbit at true anomaly nu → orbital-plane coords.
    // -----------------------------------------------------------------------
    void posAtNu(double nu, double& xOrb, double& yOrb) const {
        double absSma = std::fabs(sma);
        double r      = absSma * (ecc*ecc - 1.0);   // semi-latus rectum (hyperbola)
        if (isElliptic())
            r = sma * (1.0 - ecc * std::cos(solveKepler(
                    std::atan2(std::sqrt(1.0-ecc*ecc)*std::sin(nu),
                               ecc+std::cos(nu)), ecc)));
        // Use the conic section formula directly — works for both:
        r        = (isElliptic() ? sma*(1.0-ecc*ecc) : absSma*(ecc*ecc-1.0))
                   / (1.0 + ecc*std::cos(nu));
        xOrb = r * std::cos(nu);
        yOrb = r * std::sin(nu);
    }

    // -----------------------------------------------------------------------
    // Fit orbital elements from a state vector (position + velocity)
    // relative to body `mu = G*parentMass`.
    // Returns false if the state is degenerate.
    // -----------------------------------------------------------------------
    bool fitFromState(double prx, double pry, double prz,
                      double pvx, double pvy, double pvz,
                      double parentMass, double t)
    {
        double mu = G * parentMass;
        double r  = std::sqrt(prx*prx + pry*pry + prz*prz);
        double v2 = pvx*pvx + pvy*pvy + pvz*pvz;
        if (r < 1.0) return false;

        double hx = pry*pvz - prz*pvy;
        double hy = prz*pvx - prx*pvz;
        double hz = prx*pvy - pry*pvx;
        double h  = std::sqrt(hx*hx + hy*hy + hz*hz);
        if (h < 1.0) return false;

        double nx = -hy, ny = hx;
        double nMag = std::sqrt(nx*nx + ny*ny);

        double ex_ = (pvy*hz - pvz*hy)/mu - prx/r;
        double ey_ = (pvz*hx - pvx*hz)/mu - pry/r;
        double ez_ = (pvx*hy - pvy*hx)/mu - prz/r;
        double e   = std::sqrt(ex_*ex_ + ey_*ey_ + ez_*ez_);
        double a   = 1.0 / (2.0/r - v2/mu);
        if (a == 0.0) return false;

        double i_  = std::acos(std::clamp(hz / h, -1.0, 1.0));
        double lo  = 0.0;
        if (nMag > 1e-10) {
            lo = std::acos(std::clamp(nx / nMag, -1.0, 1.0));
            if (ny < 0.0) lo = TWO_PI - lo;
        }
        double wo = 0.0;
        if (nMag > 1e-10 && e > 1e-10) {
            double ndote = (nx*ex_ + ny*ey_) / (nMag * e);
            wo = std::acos(std::clamp(ndote, -1.0, 1.0));
            if (ez_ < 0.0) wo = TWO_PI - wo;
        }
        double nu = 0.0;
        if (e > 1e-10) {
            double rdote = (prx*ex_ + pry*ey_ + prz*ez_) / (r * e);
            nu = std::acos(std::clamp(rdote, -1.0, 1.0));
            if (prx*pvx + pry*pvy + prz*pvz < 0.0) nu = TWO_PI - nu;
        }

        // Mean anomaly at epoch
        double mae_;
        if (e < 1.0) {   // elliptic
            double cosE = (e + std::cos(nu)) / (1.0 + e*std::cos(nu));
            double sinE = std::sqrt(std::max(0.0, 1.0-e*e)) * std::sin(nu)
                          / (1.0 + e*std::cos(nu));
            double E_   = std::atan2(sinE, cosE);
            mae_ = E_ - e * std::sin(E_);
        } else {          // hyperbolic
            double cosH  = (e + std::cos(nu)) / (1.0 + e*std::cos(nu));
            double H_    = std::acosh(std::clamp(cosH, 1.0, 1e18));
            if (nu > M_PI) H_ = -H_;
            mae_ = e*std::sinh(H_) - H_;
        }

        sma   = a;
        ecc   = e;
        inc   = i_;
        lan   = lo;
        aop   = wo;
        mae   = mae_;
        epoch = t;
        valid = true;
        return true;
    }

    // -----------------------------------------------------------------------
    // Draw the orbit conic in 3-D.
    //   orig    — world-space position of the parent body
    //   camPos  — world-space camera position (for float-origin rendering)
    //   col     — line colour
    //   segments — point count (ellipse) or half-arc points (hyperbola)
    //   markApPe — draw Pe/Ap markers (only for elliptic)
    //   peRadius — parent body radius, used to size Pe/Ap dot (0 → skip dots)
    // -----------------------------------------------------------------------
    void draw(const DVec3& orig, const DVec3& camPos,
              Color col,
              int   segments  = 256,
              bool  markApPe  = false,
              double peRadius = 0.0) const
    {
        if (!valid) return;

        double Rxx,Rxy,Ryx,Ryy,Rzx,Rzy;
        buildRotMatrix(Rxx,Rxy,Ryx,Ryy,Rzx,Rzy);

        if (isElliptic()) {
            drawEllipse(orig, camPos, col, segments,
                        Rxx,Rxy,Ryx,Ryy,Rzx,Rzy,
                        markApPe, peRadius);
        } else {
            drawHyperbola(orig, camPos, col, segments,
                          Rxx,Rxy,Ryx,Ryy,Rzx,Rzy,
                          markApPe, peRadius);
        }
    }

private:

    // Inline helper: orbital coords → world DVec3
    static DVec3 orbToWorld(const DVec3& orig,
                            double xOrb, double yOrb,
                            double Rxx, double Rxy,
                            double Ryx, double Ryy,
                            double Rzx, double Rzy)
    {
        return { orig.x + Rxx*xOrb + Rxy*yOrb,
                 orig.y + Ryx*xOrb + Ryy*yOrb,
                 orig.z + Rzx*xOrb + Rzy*yOrb };
    }

    void drawEllipse(const DVec3& orig, const DVec3& camPos, Color col,
                     int segs,
                     double Rxx,double Rxy,double Ryx,double Ryy,double Rzx,double Rzy,
                     bool markApPe, double peRadius) const
    {
        Vector3 prev = {};
        for (int i = 0; i <= segs; i++) {
            double M   = TWO_PI * i / segs;
            double E   = solveKepler(M, ecc);
            double nu  = 2.0 * std::atan2(
                std::sqrt(1.0+ecc) * std::sin(E/2.0),
                std::sqrt(1.0-ecc) * std::cos(E/2.0));
            double r   = sma * (1.0 - ecc*std::cos(E));
            double xO  = r * std::cos(nu);
            double yO  = r * std::sin(nu);
            Vector3 pt = orbToWorld(orig,xO,yO,Rxx,Rxy,Ryx,Ryy,Rzx,Rzy).toVec3(camPos);
            if (i > 0) DrawLine3D(prev, pt, col);
            prev = pt;
        }

        if (markApPe && peRadius > 0.0) {
            double dotR = (float)(peRadius * 0.003);
            // Periapsis  (+x in orbital plane)
            double rPe = sma * (1.0 - ecc);
            DVec3  pe  = orbToWorld(orig, rPe, 0.0, Rxx,Rxy,Ryx,Ryy,Rzx,Rzy);
            DrawSphere(pe.toVec3(camPos), dotR, YELLOW);
            // Apoapsis  (-x in orbital plane)
            double rAp = sma * (1.0 + ecc);
            DVec3  ap  = orbToWorld(orig,-rAp, 0.0, Rxx,Rxy,Ryx,Ryy,Rzx,Rzy);
            DrawSphere(ap.toVec3(camPos), dotR, SKYBLUE);
        }
    }

    // Draw the physical (near) branch of a hyperbola.
    // The asymptotic half-angle is acos(1/e); we draw up to 95% of that.
    void drawHyperbola(const DVec3& orig, const DVec3& camPos, Color col,
                       int halfSegs,
                       double Rxx,double Rxy,double Ryx,double Ryy,double Rzx,double Rzy,
                       bool markApPe, double peRadius) const
    {
        // Max true anomaly just inside the asymptote
        double nuMax = std::acos(-1.0/ecc) * 0.97;
        int    segs  = halfSegs * 2;

        Vector3 prev = {};
        for (int i = 0; i <= segs; i++) {
            double nu  = -nuMax + (2.0*nuMax) * i / segs;
            double p   = std::fabs(sma) * (ecc*ecc - 1.0);   // semi-latus rectum
            double r   = p / (1.0 + ecc*std::cos(nu));
            double xO  = r * std::cos(nu);
            double yO  = r * std::sin(nu);
            Vector3 pt = orbToWorld(orig,xO,yO,Rxx,Rxy,Ryx,Ryy,Rzx,Rzy).toVec3(camPos);
            if (i > 0) DrawLine3D(prev, pt, col);
            prev = pt;
        }

        if (markApPe && peRadius > 0.0) {
            double dotR = (float)(peRadius * 0.003);
            double rPe  = std::fabs(sma) * (ecc - 1.0);
            DVec3  pe   = orbToWorld(orig, rPe, 0.0, Rxx,Rxy,Ryx,Ryy,Rzx,Rzy);
            DrawSphere(pe.toVec3(camPos), dotR, YELLOW);
        }
    }
};