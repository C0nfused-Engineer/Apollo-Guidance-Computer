// Apollo Guidance Computer
// Compile: g++ main.cpp -o AGC.exe -lraylib -lopengl32 -lgdi32 -lwinmm
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include "Planet.h"
#include "SystemLoader.h"
#include "DynamicBody.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

static void normalize3(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-12) { x/=len; y/=len; z/=len; }
}

static void progradeDir(const DynamicBody& db, double& ox, double& oy, double& oz) {
    ox=db.vx; oy=db.vy; oz=db.vz; normalize3(ox,oy,oz);
}
static void retrogradeDir(const DynamicBody& db, double& ox, double& oy, double& oz) {
    progradeDir(db,ox,oy,oz); ox=-ox; oy=-oy; oz=-oz;
}
static void radialOutDir(const DynamicBody& db, double& ox, double& oy, double& oz) {
    ox=db.rx; oy=db.ry; oz=db.rz; normalize3(ox,oy,oz);
}
static void normalDir(const DynamicBody& db, double& ox, double& oy, double& oz) {
    ox=db.ry*db.vz-db.rz*db.vy;
    oy=db.rz*db.vx-db.rx*db.vz;
    oz=db.rx*db.vy-db.ry*db.vx;
    normalize3(ox,oy,oz);
}

int main(int argc, char* argv[]) {
    InitWindow(1920, 1080, "Orbital Sim");
    SetTargetFPS(60);

    const char* configPath = (argc > 1) ? argv[1] : "solar_system.cfg";
    SystemLoader loader;
    try {
        loader.load(configPath);
    } catch (const std::exception& e) {
        while (!WindowShouldClose()) {
            BeginDrawing();
                ClearBackground(BLACK);
                DrawText("Failed to load config:", 20, 20, 20, RED);
                DrawText(e.what(), 20, 50, 18, RAYWHITE);
                DrawText("Press ESC to quit.", 20, 80, 18, DARKGRAY);
            EndDrawing();
        }
        CloseWindow();
        return 1;
    }

    Planet* root                    = loader.root;
    std::vector<Planet*>& allBodies = loader.allBodies;
    loader.release();

    Planet* earth = nullptr;
    for (auto* p : allBodies)
        if (p->name == "Earth") { earth=p; break; }
    if (!earth) earth = root;

    DynamicBody ship;
    ship.name         = "Ship";
    ship.mass         = 10000.0;
    ship.thrustAccel  = 10.0;
    ship.color        = GREEN;
    ship.previewSteps = 500;
    ship.previewDt    = 60.0;
    ship.initCircularOrbit(earth, 0.0, 400e3, 0.0);

    float  camYaw   = 0.0f;
    float  camPitch = 80.0f;
    double camDist  = 1.5e11;

    int camTargetIndex = 0;
    int totalTargets   = (int)allBodies.size() + 1;

    Camera3D camera   = {};
    camera.fovy       = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    camera.up         = { 0.0f, 0.0f, 1.0f };

    const int   NUM_MODES   = 6;
    const char* modeNames[] = {
        "Prograde","Retrograde","Radial Out","Radial In","Normal","Anti-Normal"
    };
    int thrustMode = 0;

    double simTime   = 0.0;
    double timeScale = 1.0;
    bool   paused    = false;

    while (!WindowShouldClose()) {
        double dt = GetFrameTime();

        if (IsKeyPressed(KEY_UP))    timeScale *= 10.0;
        if (IsKeyPressed(KEY_DOWN))  timeScale /= 10.0;
        if (IsKeyPressed(KEY_SPACE)) paused = !paused;
        if (IsKeyPressed(KEY_TAB))   camTargetIndex = (camTargetIndex+1) % totalTargets;
        if (IsKeyPressed(KEY_M))     thrustMode = (thrustMode+1) % NUM_MODES;

        double effectiveTScale = paused ? 0.0 : timeScale;
        if (ship.thrusting && effectiveTScale > 1000.0) effectiveTScale = 1000.0;

        bool wantThrust   = IsKeyDown(KEY_Z);
        bool wasThrusting = ship.thrusting;
        ship.thrusting    = wantThrust;

        if (wantThrust) {
            double tdx=0,tdy=0,tdz=0;
            switch (thrustMode) {
                case 0: progradeDir  (ship,tdx,tdy,tdz); break;
                case 1: retrogradeDir(ship,tdx,tdy,tdz); break;
                case 2: radialOutDir (ship,tdx,tdy,tdz); break;
                case 3: radialOutDir (ship,tdx,tdy,tdz); tdx=-tdx; tdy=-tdy; tdz=-tdz; break;
                case 4: normalDir    (ship,tdx,tdy,tdz); break;
                case 5: normalDir    (ship,tdx,tdy,tdz); tdx=-tdx; tdy=-tdy; tdz=-tdz; break;
            }
            ship.tdx=tdx; ship.tdy=tdy; ship.tdz=tdz;
        }

        if (wasThrusting && !ship.thrusting)
            ship.onThrustStop(simTime);

        if (!paused) {
            double simDt = dt * effectiveTScale;
            simTime += simDt;
            root->update(simTime);
            ship.update(simTime, simDt);
        }

        DVec3 targetPos;
        if (camTargetIndex < (int)allBodies.size())
            targetPos = allBodies[camTargetIndex]->position;
        else
            targetPos = ship.worldPos;

        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            Vector2 delta = GetMouseDelta();
            camYaw   += delta.x * 0.3f;
            camPitch += delta.y * 0.3f;
            camPitch  = Clamp(camPitch, 5.0f, 89.9f);
        }
        camDist -= GetMouseWheelMove() * camDist * 0.1;
        if (camDist < 1.0) camDist = 1.0;

        float pitchRad = camPitch * DEG2RAD;
        float yawRad   = camYaw   * DEG2RAD;

        DVec3 camWorldPos = {
            targetPos.x + camDist * cos(pitchRad) * sin(yawRad),
            targetPos.y + camDist * cos(pitchRad) * cos(yawRad),
            targetPos.z + camDist * sin(pitchRad)
        };

        camera.position = { 0.0f, 0.0f, 0.0f };
        camera.target   = targetPos.toVec3(camWorldPos);
        camera.up       = { 0.0f, 0.0f, 1.0f };
        rlSetClipPlanes(camDist * 0.0001, camDist * 1000.0);

        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode3D(camera);
            root->draw(camWorldPos);
            ship.draw(camWorldPos);

            const double AU       = 1.496e11;
            const double gridExt  = 6.0 * AU;
            const double gridStep = 0.1 * AU;
            Color gridColor       = { 30, 30, 30, 255 };
            for (double i=-gridExt; i<=gridExt+gridStep*0.5; i+=gridStep) {
                DVec3 a1={-gridExt,i,0.0}, b1={gridExt,i,0.0};
                DVec3 a2={i,-gridExt,0.0}, b2={i,gridExt,0.0};
                DrawLine3D(a1.toVec3(camWorldPos), b1.toVec3(camWorldPos), gridColor);
                DrawLine3D(a2.toVec3(camWorldPos), b2.toVec3(camWorldPos), gridColor);
            }
        EndMode3D();

        // HUD
        int sh = GetScreenHeight();
        DrawText("UP/DOWN: timescale  |  TAB: target  |  M: thrust mode  |  Z: fire  |  SPACE: pause  |  RMB: rotate  |  Scroll: zoom",
                 10, 10, 16, RAYWHITE);

        char buf[256];
        const char* targetName = (camTargetIndex < (int)allBodies.size())
                               ? allBodies[camTargetIndex]->name.c_str()
                               : ship.name.c_str();
        sprintf(buf, "Target: %s", targetName);
        DrawText(buf, 10, 30, 18, YELLOW);

        if      (timeScale>=86400.0) sprintf(buf,"Time scale: %.1f days/s",timeScale/86400.0);
        else if (timeScale>=3600.0)  sprintf(buf,"Time scale: %.1f hrs/s", timeScale/3600.0);
        else if (timeScale>=60.0)    sprintf(buf,"Time scale: %.1f min/s", timeScale/60.0);
        else                         sprintf(buf,"Time scale: %.0fx",       timeScale);
        DrawText(buf, 10, 50, 18, RAYWHITE);

        long long totalSec = (long long)simTime;
        sprintf(buf,"Sim time: %lld d %02lld h %02lld m %02lld s",
            (long long)(totalSec/86400),
            (long long)((totalSec%86400)/3600),
            (long long)((totalSec%3600)/60),
            (long long)(totalSec%60));
        DrawText(buf, 10, 70, 18, RAYWHITE);

        if      (camDist>=1.496e11) sprintf(buf,"Cam dist: %.3f AU", camDist/1.496e11);
        else if (camDist>=1e6)      sprintf(buf,"Cam dist: %.0f km", camDist/1e3);
        else if (camDist>=1e3)      sprintf(buf,"Cam dist: %.1f km", camDist/1e3);
        else                        sprintf(buf,"Cam dist: %.1f m",  camDist);
        DrawText(buf, 10, 90, 18, RAYWHITE);

        if (paused) DrawText("[ PAUSED ]", 10, 115, 20, RED);

        // Ship panel
        int py = sh - 155;
        DrawRectangle(0, py-6, 420, 161, {0,0,0,180});

        DrawText("-- SHIP --", 10, py, 18, GREEN); py+=22;

        sprintf(buf,"Mode: %s  (M to cycle)", modeNames[thrustMode]);
        DrawText(buf, 10, py, 16, RAYWHITE); py+=20;

        sprintf(buf,"Thrust: %s  (hold Z)", ship.thrusting ? "FIRING" : "off");
        DrawText(buf, 10, py, 16, ship.thrusting ? ORANGE : GRAY); py+=20;

        double alt = ship.altitude();
        if      (alt>1e6) sprintf(buf,"Alt:   %.0f km", alt/1e3);
        else if (alt>0)   sprintf(buf,"Alt:   %.1f m",  alt);
        else              sprintf(buf,"Alt:   IMPACT");
        DrawText(buf, 10, py, 16, alt>0 ? RAYWHITE : RED); py+=20;

        sprintf(buf,"Speed: %.1f m/s  (%.3f km/s)", ship.speed(), ship.speed()/1e3);
        DrawText(buf, 10, py, 16, RAYWHITE); py+=20;

        if (ship.orbit.valid) {
            double sma = ship.orbit.sma;
            double ecc = ship.orbit.ecc;
            double pe  = (sma*(1.0-ecc) - ship.parent->radius) / 1e3;
            double ap  = (sma*(1.0+ecc) - ship.parent->radius) / 1e3;
            sprintf(buf,"Pe: %.0f km   Ap: %.0f km", pe, ap);
            DrawText(buf, 10, py, 16, RAYWHITE); py+=20;

            sprintf(buf,"SMA: %.0f km   Ecc: %.5f", sma/1e3, ecc);
            DrawText(buf, 10, py, 16, RAYWHITE);
        } else {
            DrawText("SMA: ---   Ecc: ---", 10, py, 16, RAYWHITE);
        }

        if (ship.thrusting) {
            DrawRectangle(10, py-44, 14, 14, ORANGE);
            DrawText("burn arc", 28, py-44, 15, RAYWHITE);
            DrawRectangle(10, py-26, 14, 14, {0,220,120,200});
            DrawText("projected coast orbit", 28, py-26, 15, RAYWHITE);
        }

        EndDrawing();
    }

    for (auto* p : allBodies) delete p;
    CloseWindow();
    return 0;
}