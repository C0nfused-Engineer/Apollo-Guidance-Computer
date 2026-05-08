#pragma once
#include "Planet.h"
#include "Orbit.h"
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>

class DynamicBody {
public:
    std::string name;
    double      mass;
    double      thrustAccel;
    Color       color = GREEN;

    Planet* parent = nullptr;

    double rx = 0, ry = 0, rz = 0;
    double vx = 0, vy = 0, vz = 0;

    Orbit orbit;   // ← replaces the old Elements inner struct

    bool thrusting = false;
    bool valid     = false;

    double tdx = 0, tdy = 0, tdz = 0;

    DVec3 worldPos;

    int    previewSteps = 540;
    double previewDt    = 10.0;

    int    patchMaxPatches = 5;
    double patchStepDt     = 30.0;
    int    patchMaxSteps   = 2000;

    static constexpr double G      = 6.674e-11;
    static constexpr double TWO_PI = 2.0 * M_PI;

    // -----------------------------------------------------------------------
    void initFromState(Planet* soi, double t0,
                       double px, double py, double pz,
                       double pvx, double pvy, double pvz)
    {
        parent = soi;
        rx=px; ry=py; rz=pz;
        vx=pvx; vy=pvy; vz=pvz;
        orbit.fitFromState(rx,ry,rz,vx,vy,vz, parent->mass, t0);
        valid = true;
    }

    void initCircularOrbit(Planet* soi, double t0,
                           double altitudeM, double lonDeg = 0.0)
    {
        parent = soi;
        double r  = soi->radius + altitudeM;
        double v  = std::sqrt(G * soi->mass / r);
        double lo = lonDeg * M_PI / 180.0;
        rx = r*std::cos(lo);  ry = r*std::sin(lo);  rz = 0.0;
        vx = -v*std::sin(lo); vy = v*std::cos(lo);  vz = 0.0;
        orbit.fitFromState(rx,ry,rz,vx,vy,vz, parent->mass, t0);
        valid = true;
    }

    // -----------------------------------------------------------------------
    void update(double t, double dt)
    {
        if (!valid || !parent) return;

        const double maxStep = 10.0;
        double remaining = dt;
        while (remaining > 0.0) {
            double step = std::min(remaining, maxStep);
            integrateRK4(step, rx,ry,rz, vx,vy,vz);
            remaining -= step;
        }

        orbit.fitFromState(rx,ry,rz,vx,vy,vz, parent->mass, t);
        checkAndUpdateSOI();

        worldPos = {
            parent->position.x + rx,
            parent->position.y + ry,
            parent->position.z + rz
        };
    }

    void onThrustStop(double t) {
        orbit.fitFromState(rx,ry,rz,vx,vy,vz, parent->mass, t);
    }

    // -----------------------------------------------------------------------
    void draw(const DVec3& camPos) const {
        if (!valid) return;

        const DVec3 orig    = parent->position;
        DVec3       shipWorld = { orig.x+rx, orig.y+ry, orig.z+rz };
        Vector3     pos       = shipWorld.toVec3(camPos);

        DrawSphere(pos, (float)(parent->radius * 0.005), color);
        DrawSphereWires(pos, (float)(parent->radius * 0.005), 8, 8, color);

        double spd = speed();
        if (spd > 0.0) {
            double scale  = parent->soiRadius * 0.02 / spd;
            DVec3  velTip = { shipWorld.x+vx*scale,
                              shipWorld.y+vy*scale,
                              shipWorld.z+vz*scale };
            DrawLine3D(pos, velTip.toVec3(camPos), LIME);
        }

        if (thrusting) {
            double tscale    = parent->soiRadius * 0.03;
            DVec3  thrustTip = { shipWorld.x+tdx*tscale,
                                 shipWorld.y+tdy*tscale,
                                 shipWorld.z+tdz*tscale };
            DrawLine3D(pos, thrustTip.toVec3(camPos), ORANGE);
        }

        if (thrusting)
            drawThrustTrajectory(camPos, orig);
        else
            drawPatchedConics(camPos);
    }

    double altitude() const {
        return std::sqrt(rx*rx + ry*ry + rz*rz) - parent->radius;
    }
    double speed() const {
        return std::sqrt(vx*vx + vy*vy + vz*vz);
    }

private:

    // -----------------------------------------------------------------------
    // Patched-conics prediction
    // -----------------------------------------------------------------------
    struct PatchResult {
        Planet* soi    = nullptr;
        Orbit   orbit;
        struct Pt3 { double x,y,z; };
        std::vector<Pt3> arcPts3;   // used only for the hyperbolic arc
    };

    // -----------------------------------------------------------------------
    double solveKepler(double M, double e) const {
        return orbit.solveKepler(M, e);
    }

    // -----------------------------------------------------------------------
    void integrateRK4(double dt,
                      double& prx, double& pry, double& prz,
                      double& pvx_, double& pvy_, double& pvz_,
                      bool thrustOn = true) const
    {
        auto deriv = [&](double irx, double iry, double irz,
                         double ivx, double ivy, double ivz,
                         double& oax, double& oay, double& oaz,
                         double& odx, double& ody, double& odz)
        {
            double mu = G * parent->mass;
            double r  = std::sqrt(irx*irx+iry*iry+irz*irz);
            if (r < 1.0) r = 1.0;
            double r3 = r*r*r;
            oax = -mu*irx/r3 + (thrustOn && thrusting ? thrustAccel*tdx : 0.0);
            oay = -mu*iry/r3 + (thrustOn && thrusting ? thrustAccel*tdy : 0.0);
            oaz = -mu*irz/r3 + (thrustOn && thrusting ? thrustAccel*tdz : 0.0);
            odx = ivx; ody = ivy; odz = ivz;
        };

        double ax1,ay1,az1,dx1,dy1,dz1;
        deriv(prx,pry,prz,pvx_,pvy_,pvz_,ax1,ay1,az1,dx1,dy1,dz1);
        double ax2,ay2,az2,dx2,dy2,dz2;
        deriv(prx+dx1*dt*.5,pry+dy1*dt*.5,prz+dz1*dt*.5,
              pvx_+ax1*dt*.5,pvy_+ay1*dt*.5,pvz_+az1*dt*.5,
              ax2,ay2,az2,dx2,dy2,dz2);
        double ax3,ay3,az3,dx3,dy3,dz3;
        deriv(prx+dx2*dt*.5,pry+dy2*dt*.5,prz+dz2*dt*.5,
              pvx_+ax2*dt*.5,pvy_+ay2*dt*.5,pvz_+az2*dt*.5,
              ax3,ay3,az3,dx3,dy3,dz3);
        double ax4,ay4,az4,dx4,dy4,dz4;
        deriv(prx+dx3*dt,pry+dy3*dt,prz+dz3*dt,
              pvx_+ax3*dt,pvy_+ay3*dt,pvz_+az3*dt,
              ax4,ay4,az4,dx4,dy4,dz4);

        prx  += (dt/6.0)*(dx1+2*dx2+2*dx3+dx4);
        pry  += (dt/6.0)*(dy1+2*dy2+2*dy3+dy4);
        prz  += (dt/6.0)*(dz1+2*dz2+2*dz3+dz4);
        pvx_ += (dt/6.0)*(ax1+2*ax2+2*ax3+ax4);
        pvy_ += (dt/6.0)*(ay1+2*ay2+2*ay3+ay4);
        pvz_ += (dt/6.0)*(az1+2*az2+2*az3+az4);
    }

    static void integrateRK4Coast(double dt,
                                  double& prx, double& pry, double& prz,
                                  double& pvx_, double& pvy_, double& pvz_,
                                  Planet* soi)
    {
        auto deriv = [&](double irx, double iry, double irz,
                         double ivx, double ivy, double ivz,
                         double& oax, double& oay, double& oaz,
                         double& odx, double& ody, double& odz)
        {
            double mu = G * soi->mass;
            double r  = std::sqrt(irx*irx+iry*iry+irz*irz);
            if (r < 1.0) r = 1.0;
            double r3 = r*r*r;
            oax=-mu*irx/r3; oay=-mu*iry/r3; oaz=-mu*irz/r3;
            odx=ivx; ody=ivy; odz=ivz;
        };

        double ax1,ay1,az1,dx1,dy1,dz1;
        deriv(prx,pry,prz,pvx_,pvy_,pvz_,ax1,ay1,az1,dx1,dy1,dz1);
        double ax2,ay2,az2,dx2,dy2,dz2;
        deriv(prx+dx1*dt*.5,pry+dy1*dt*.5,prz+dz1*dt*.5,
              pvx_+ax1*dt*.5,pvy_+ay1*dt*.5,pvz_+az1*dt*.5,
              ax2,ay2,az2,dx2,dy2,dz2);
        double ax3,ay3,az3,dx3,dy3,dz3;
        deriv(prx+dx2*dt*.5,pry+dy2*dt*.5,prz+dz2*dt*.5,
              pvx_+ax2*dt*.5,pvy_+ay2*dt*.5,pvz_+az2*dt*.5,
              ax3,ay3,az3,dx3,dy3,dz3);
        double ax4,ay4,az4,dx4,dy4,dz4;
        deriv(prx+dx3*dt,pry+dy3*dt,prz+dz3*dt,
              pvx_+ax3*dt,pvy_+ay3*dt,pvz_+az3*dt,
              ax4,ay4,az4,dx4,dy4,dz4);

        prx  += (dt/6.0)*(dx1+2*dx2+2*dx3+dx4);
        pry  += (dt/6.0)*(dy1+2*dy2+2*dy3+dy4);
        prz  += (dt/6.0)*(dz1+2*dz2+2*dz3+dz4);
        pvx_ += (dt/6.0)*(ax1+2*ax2+2*ax3+ax4);
        pvy_ += (dt/6.0)*(ay1+2*ay2+2*ay3+ay4);
        pvz_ += (dt/6.0)*(az1+2*az2+2*az3+az4);
    }

    // -----------------------------------------------------------------------
    // SOI velocity of a body (orbital velocity relative to its parent)
    // -----------------------------------------------------------------------
    static void soiVelocity(Planet* body,
                            double& vbx, double& vby, double& vbz)
    {
        Planet* par = body->parent;
        if (!par) { vbx=vby=vbz=0; return; }
        double lx = body->position.x - par->position.x;
        double ly = body->position.y - par->position.y;
        double lz = body->position.z - par->position.z;
        double sr = std::sqrt(lx*lx+ly*ly+lz*lz);
        if (sr < 1.0) { vbx=vby=vbz=0; return; }
        double mu = G * par->mass;
        double sv = std::sqrt(std::abs(mu*(2.0/sr - 1.0/body->orbit.sma)));
        double rxn=lx/sr, ryn=ly/sr;
        double pvx_=-ryn, pvy_=rxn;
        double pl = std::sqrt(pvx_*pvx_+pvy_*pvy_);
        if (pl>1e-10) { pvx_/=pl; pvy_/=pl; }
        vbx=sv*pvx_; vby=sv*pvy_; vbz=0.0;
    }

    static void exitSOI(Planet* curSoi,
                        double lx,double ly,double lz,
                        double lvx,double lvy,double lvz,
                        double& ox,double& oy,double& oz,
                        double& ovx,double& ovy,double& ovz)
    {
        Planet* par = curSoi->parent;
        double soiLx = curSoi->position.x - par->position.x;
        double soiLy = curSoi->position.y - par->position.y;
        double soiLz = curSoi->position.z - par->position.z;
        ox=soiLx+lx; oy=soiLy+ly; oz=soiLz+lz;
        double vbx,vby,vbz;
        soiVelocity(curSoi,vbx,vby,vbz);
        ovx=lvx+vbx; ovy=lvy+vby; ovz=lvz+vbz;
    }

    static void enterSOI(Planet* child, Planet* par,
                         double lx,double ly,double lz,
                         double lvx,double lvy,double lvz,
                         double& ox,double& oy,double& oz,
                         double& ovx,double& ovy,double& ovz)
    {
        double cx=child->position.x-par->position.x;
        double cy=child->position.y-par->position.y;
        double cz=child->position.z-par->position.z;
        ox=lx-cx; oy=ly-cy; oz=lz-cz;
        double vbx,vby,vbz;
        soiVelocity(child,vbx,vby,vbz);
        ovx=lvx-vbx; ovy=lvy-vby; ovz=lvz-vbz;
    }

    static Planet* findChildSOI(Planet* soi,
                                double wx,double wy,double wz)
    {
        for (Planet* child : soi->children) {
            double dx=wx-child->position.x;
            double dy=wy-child->position.y;
            double dz=wz-child->position.z;
            if (std::sqrt(dx*dx+dy*dy+dz*dz) < child->soiRadius)
                return child;
        }
        return nullptr;
    }

    void checkAndUpdateSOI()
    {
        if (!parent) return;
        double wx=parent->position.x+rx;
        double wy=parent->position.y+ry;
        double wz=parent->position.z+rz;

        Planet* child=findChildSOI(parent,wx,wy,wz);
        if (child) {
            double ox,oy,oz,ovx,ovy,ovz;
            enterSOI(child,parent,rx,ry,rz,vx,vy,vz,ox,oy,oz,ovx,ovy,ovz);
            parent=child; rx=ox; ry=oy; rz=oz; vx=ovx; vy=ovy; vz=ovz;
            return;
        }
        if (parent->parent) {
            double dist=std::sqrt(rx*rx+ry*ry+rz*rz);
            if (dist>parent->soiRadius) {
                double ox,oy,oz,ovx,ovy,ovz;
                exitSOI(parent,rx,ry,rz,vx,vy,vz,ox,oy,oz,ovx,ovy,ovz);
                parent=parent->parent; rx=ox; ry=oy; rz=oz; vx=ovx; vy=ovy; vz=ovz;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Build patched-conic patches for the coast trajectory
    // -----------------------------------------------------------------------
    std::vector<PatchResult> buildPatches() const
    {
        std::vector<PatchResult> patches;
        Planet* curSoi = parent;
        double  crx=rx,cry=ry,crz=rz;
        double  cvx=vx,cvy=vy,cvz=vz;

        for (int patch=0; patch<patchMaxPatches; patch++) {
            PatchResult pr;
            pr.soi = curSoi;

            if (!pr.orbit.fitFromState(crx,cry,crz,cvx,cvy,cvz,
                                       curSoi->mass, 0.0))
                break;

            bool elliptic = pr.orbit.isElliptic();

            if (elliptic && curSoi->parent==nullptr) {
                patches.push_back(pr);
                break;
            }

            if (elliptic) {
                double period = TWO_PI * std::sqrt(
                    std::pow(std::abs(pr.orbit.sma), 3.0)
                    / (G * curSoi->mass));
                int steps = std::min(patchMaxSteps,
                                     (int)(period/patchStepDt)+1);

                double prx2=crx,pry2=cry,prz2=crz;
                double pvx2=cvx,pvy2=cvy,pvz2=cvz;
                bool    exited=false;
                Planet* nextSoi=nullptr;
                double  nrx=0,nry=0,nrz=0,nvx=0,nvy=0,nvz=0;

                for (int s=0; s<steps; s++) {
                    integrateRK4Coast(patchStepDt,
                                      prx2,pry2,prz2,pvx2,pvy2,pvz2, curSoi);
                    double wx2=curSoi->position.x+prx2;
                    double wy2=curSoi->position.y+pry2;
                    double wz2=curSoi->position.z+prz2;

                    Planet* ch=findChildSOI(curSoi,wx2,wy2,wz2);
                    if (ch) {
                        enterSOI(ch,curSoi,prx2,pry2,prz2,pvx2,pvy2,pvz2,
                                 nrx,nry,nrz,nvx,nvy,nvz);
                        nextSoi=ch; exited=true; break;
                    }
                    if (curSoi->parent) {
                        double d=std::sqrt(prx2*prx2+pry2*pry2+prz2*prz2);
                        if (d>curSoi->soiRadius) {
                            exitSOI(curSoi,prx2,pry2,prz2,pvx2,pvy2,pvz2,
                                    nrx,nry,nrz,nvx,nvy,nvz);
                            nextSoi=curSoi->parent; exited=true; break;
                        }
                    }
                }

                patches.push_back(pr);
                if (!exited) break;
                curSoi=nextSoi;
                crx=nrx; cry=nry; crz=nrz;
                cvx=nvx; cvy=nvy; cvz=nvz;

            } else {
                // Hyperbolic: integrate arc until SOI exit
                double prx2=crx,pry2=cry,prz2=crz;
                double pvx2=cvx,pvy2=cvy,pvz2=cvz;
                bool    exited=false;
                Planet* nextSoi=nullptr;
                double  nrx=0,nry=0,nrz=0,nvx=0,nvy=0,nvz=0;

                for (int s=0; s<patchMaxSteps; s++) {
                    integrateRK4Coast(patchStepDt,
                                      prx2,pry2,prz2,pvx2,pvy2,pvz2, curSoi);
                    // No need to store arc pts — Orbit::draw handles hyperbola
                    double wx2=curSoi->position.x+prx2;
                    double wy2=curSoi->position.y+pry2;
                    double wz2=curSoi->position.z+prz2;

                    Planet* ch=findChildSOI(curSoi,wx2,wy2,wz2);
                    if (ch) {
                        enterSOI(ch,curSoi,prx2,pry2,prz2,pvx2,pvy2,pvz2,
                                 nrx,nry,nrz,nvx,nvy,nvz);
                        nextSoi=ch; exited=true; break;
                    }
                    if (curSoi->parent) {
                        double d=std::sqrt(prx2*prx2+pry2*pry2+prz2*prz2);
                        if (d>curSoi->soiRadius) {
                            exitSOI(curSoi,prx2,pry2,prz2,pvx2,pvy2,pvz2,
                                    nrx,nry,nrz,nvx,nvy,nvz);
                            nextSoi=curSoi->parent; exited=true; break;
                        }
                    }
                }

                patches.push_back(pr);
                if (!exited || !nextSoi) break;
                curSoi=nextSoi;
                crx=nrx; cry=nry; crz=nrz;
                cvx=nvx; cvy=nvy; cvz=nvz;
            }
        }
        return patches;
    }

    // -----------------------------------------------------------------------
    void drawPatchedConics(const DVec3& camPos) const
    {
        static const Color patchCols[] = {
            {   0, 200, 100, 180 },
            {  80, 160, 255, 180 },
            { 255, 200,  50, 180 },
            { 200,  80, 255, 180 },
            { 255, 100, 100, 180 },
        };
        std::vector<PatchResult> patches = buildPatches();
        for (int pi=0; pi<(int)patches.size(); pi++) {
            const PatchResult& pr = patches[pi];
            Color col  = patchCols[pi % 5];
            DVec3 orig = pr.soi->position;
            // Orbit::draw handles both elliptic and hyperbolic automatically
            pr.orbit.draw(orig, camPos, col,
                          256,          // segments
                          true,         // mark Pe/Ap
                          pr.soi->radius);
        }
    }

    // -----------------------------------------------------------------------
    void drawThrustTrajectory(const DVec3& camPos, const DVec3& orig) const
    {
        if (!valid || !parent) return;

        // Projected coast orbit after burn ends
        Orbit projected;
        if (projected.fitFromState(rx,ry,rz,vx,vy,vz,parent->mass,0.0)) {
            projected.draw(orig, camPos,
                           { 0, 220, 120, 200 },
                           256, true, parent->radius);
        }

        // Short numerical burn arc (orange gradient)
        const int    arcSteps = 60;
        const double arcDt    = 10.0;
        double prx=rx,pry=ry,prz=rz,pvx_=vx,pvy_=vy,pvz_=vz;
        Vector3 prev = DVec3{orig.x+rx,orig.y+ry,orig.z+rz}.toVec3(camPos);

        for (int i=0; i<arcSteps; i++) {
            integrateRK4(arcDt, prx,pry,prz, pvx_,pvy_,pvz_, true);
            if (std::sqrt(prx*prx+pry*pry+prz*prz) < parent->radius) break;
            float  t_  = (float)i / (float)(arcSteps-1);
            Color  col = { 255, (unsigned char)(165+(int)(90.0f*t_)), 0, 200 };
            Vector3 pt = DVec3{orig.x+prx,orig.y+pry,orig.z+prz}.toVec3(camPos);
            DrawLine3D(prev, pt, col);
            prev = pt;
        }
    }
};