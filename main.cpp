#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include "Planet.h"
#include "SystemLoader.h"
#include "DynamicBody.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// Thrust direction helpers — all return unit vectors in parent-inertial frame
// ---------------------------------------------------------------------------

static void normalize3(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-12) { x /= len; y /= len; z /= len; }
}

// Prograde: velocity direction
static void progradeDir(const DynamicBody& db, double& ox, double& oy, double& oz) {
    ox = db.vx; oy = db.vy; oz = db.vz;
    normalize3(ox, oy, oz);
}

// Retrograde
static void retrogradeDir(const DynamicBody& db, double& ox, double& oy, double& oz) {
    progradeDir(db, ox, oy, oz);
    ox = -ox; oy = -oy; oz = -oz;
}

// Radial outward: position direction
static void radialOutDir(const DynamicBody& db, double& ox, double& oy, double& oz) {
    ox = db.rx; oy = db.ry; oz = db.rz;
    normalize3(ox, oy, oz);
}

// Normal to orbital plane: h = r × v
static void normalDir(const DynamicBody& db, double& ox, double& oy, double& oz) {
    ox = db.ry*db.vz - db.rz*db.vy;
    oy = db.rz*db.vx - db.rx*db.vz;
    oz = db.rx*db.vy - db.ry*db.vx;
    normalize3(ox, oy, oz);
}

int main(int argc, char* argv[]) {
    // --- Fullscreen at native resolution ---
    InitWindow(1920, 1080, "Orbital Sim");
    SetTargetFPS(60);

    // --- Load system ---
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

    // --- Spawn ship in circular orbit around Earth at 400 km altitude ---
    Planet* earth = nullptr;
    for (auto* p : allBodies)
        if (p->name == "Earth") { earth = p; break; }
    if (!earth) earth = root;

    DynamicBody ship;
    ship.name        = "Ship";
    ship.mass        = 10000.0;   // kg
    ship.thrustAccel = 10.0;      // m/s² (~1g, fun but not unrealistic for a rocket)
    ship.color       = GREEN;
    ship.initCircularOrbit(earth, 0.0, 400e3, 0.0);

    // --- Camera ---
    float  camYaw   = 0.0f;
    float  camPitch = 80.0f;
    double camDist  = 1.5e11;

    // Targets: 0..N-1 = planets, N = ship
    int camTargetIndex = 0;
    int totalTargets   = (int)allBodies.size() + 1;

    Camera3D camera   = {};
    camera.fovy       = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    camera.up         = { 0.0f, 0.0f, 1.0f };

    // --- Thrust mode ---
    // 0=prograde  1=retrograde  2=radial-out  3=radial-in  4=normal  5=antinormal
    const int   NUM_MODES   = 6;
    const char* modeNames[] = {
        "Prograde", "Retrograde", "Radial Out", "Radial In", "Normal", "Anti-Normal"
    };
    int thrustMode = 0;

    // --- Time ---
    double simTime   = 0.0;
    double timeScale = 1.0;
    bool   paused    = false;

    while (!WindowShouldClose()) {
        double dt = GetFrameTime();

        // --- Input ---
        if (IsKeyPressed(KEY_UP))    timeScale *= 10.0;
        if (IsKeyPressed(KEY_DOWN))  timeScale /= 10.0;
        if (IsKeyPressed(KEY_SPACE)) paused = !paused;
        if (IsKeyPressed(KEY_TAB))   camTargetIndex = (camTargetIndex + 1) % totalTargets;
        if (IsKeyPressed(KEY_M))     thrustMode = (thrustMode + 1) % NUM_MODES;

        // Clamp time scale during burns so RK4 sub-steps stay manageable
        double effectiveTScale = paused ? 0.0 : timeScale;
        if (ship.thrusting && effectiveTScale > 1000.0) effectiveTScale = 1000.0;

        // --- Thrust ---
        bool wantThrust = IsKeyDown(KEY_Z);
        bool wasThrusting = ship.thrusting;
        ship.thrusting = wantThrust;

        if (wantThrust) {
            double tdx = 0, tdy = 0, tdz = 0;
            switch (thrustMode) {
                case 0: progradeDir  (ship, tdx, tdy, tdz); break;
                case 1: retrogradeDir(ship, tdx, tdy, tdz); break;
                case 2: radialOutDir (ship, tdx, tdy, tdz); break;
                case 3: radialOutDir (ship, tdx, tdy, tdz); tdx=-tdx; tdy=-tdy; tdz=-tdz; break;
                case 4: normalDir    (ship, tdx, tdy, tdz); break;
                case 5: normalDir    (ship, tdx, tdy, tdz); tdx=-tdx; tdy=-tdy; tdz=-tdz; break;
            }
            ship.tdx = tdx; ship.tdy = tdy; ship.tdz = tdz;
        }

        // On burn-end: refit osculating elements
        if (wasThrusting && !ship.thrusting)
            ship.onThrustStop(simTime);

        // --- Simulate ---
        if (!paused) {
            double simDt = dt * effectiveTScale;
            simTime += simDt;
            root->update(simTime);
            ship.update(simTime, simDt);
        }

        // --- Camera target position ---
        DVec3 targetPos;
        if (camTargetIndex < (int)allBodies.size())
            targetPos = allBodies[camTargetIndex]->position;
        else
            targetPos = ship.worldPos;

        // Camera orbit
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

        // --- Draw ---
        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode3D(camera);
            root->draw(camWorldPos);
            ship.draw(camWorldPos);

            // Ecliptic reference grid (XY plane, Z=0)
            const double AU       = 1.496e11;
            const double gridExt  = 6.0 * AU;
            const double gridStep = 0.1 * AU;
            Color gridColor       = { 30, 30, 30, 255 };
            for (double i = -gridExt; i <= gridExt + gridStep * 0.5; i += gridStep) {
                DVec3 a1 = { -gridExt, i, 0.0 }, b1 = { gridExt, i, 0.0 };
                DVec3 a2 = { i, -gridExt, 0.0 }, b2 = { i, gridExt, 0.0 };
                DrawLine3D(a1.toVec3(camWorldPos), b1.toVec3(camWorldPos), gridColor);
                DrawLine3D(a2.toVec3(camWorldPos), b2.toVec3(camWorldPos), gridColor);
            }
        EndMode3D();

        // --- HUD ---
        int sw = GetScreenWidth();
        int sh = GetScreenHeight();
        (void)sw;

        DrawText("UP/DOWN: timescale  |  TAB: next target  |  M: thrust mode  |  Z: fire  |  SPACE: pause  |  RMB: rotate  |  Scroll: zoom",
                 10, 10, 16, RAYWHITE);

        char buf[256];
        const char* targetName = (camTargetIndex < (int)allBodies.size())
                               ? allBodies[camTargetIndex]->name.c_str()
                               : ship.name.c_str();
        sprintf(buf, "Target: %s", targetName);
        DrawText(buf, 10, 30, 18, YELLOW);

        if (timeScale >= 86400.0)      sprintf(buf, "Time scale: %.1f days/s", timeScale / 86400.0);
        else if (timeScale >= 3600.0)  sprintf(buf, "Time scale: %.1f hrs/s",  timeScale / 3600.0);
        else if (timeScale >= 60.0)    sprintf(buf, "Time scale: %.1f min/s",  timeScale / 60.0);
        else                           sprintf(buf, "Time scale: %.0fx",        timeScale);
        DrawText(buf, 10, 50, 18, RAYWHITE);

        long long totalSec = (long long)simTime;
        sprintf(buf, "Sim time: %lld d %02lld h %02lld m %02lld s",
            (long long)(totalSec / 86400),
            (long long)((totalSec % 86400) / 3600),
            (long long)((totalSec % 3600)  / 60),
            (long long)(totalSec % 60));
        DrawText(buf, 10, 70, 18, RAYWHITE);

        if (camDist >= 1.496e11)  sprintf(buf, "Cam dist: %.3f AU",  camDist / 1.496e11);
        else if (camDist >= 1e6)  sprintf(buf, "Cam dist: %.0f km",  camDist / 1e3);
        else if (camDist >= 1e3)  sprintf(buf, "Cam dist: %.1f km",  camDist / 1e3);
        else                      sprintf(buf, "Cam dist: %.1f m",   camDist);
        DrawText(buf, 10, 90, 18, RAYWHITE);

        if (paused) DrawText("[ PAUSED ]", 10, 115, 20, RED);

        // --- Ship status panel (bottom-left) ---
        int py = sh - 130;
        DrawRectangle(0, py - 6, 380, 136, { 0, 0, 0, 180 });

        DrawText("-- SHIP --", 10, py, 18, GREEN); py += 22;

        sprintf(buf, "Mode: %s  (M to cycle)", modeNames[thrustMode]);
        DrawText(buf, 10, py, 16, RAYWHITE); py += 20;

        sprintf(buf, "Thrust: %s  (hold Z)", ship.thrusting ? "FIRING" : "off");
        DrawText(buf, 10, py, 16, ship.thrusting ? ORANGE : GRAY); py += 20;

        double alt = ship.altitude();
        if (alt > 1e6)       sprintf(buf, "Alt:   %.0f km",   alt / 1e3);
        else if (alt > 0)    sprintf(buf, "Alt:   %.1f m",    alt);
        else                 sprintf(buf, "Alt:   SUBORBITAL");
        DrawText(buf, 10, py, 16, alt > 0 ? RAYWHITE : RED); py += 20;

        sprintf(buf, "Speed: %.1f m/s  (%.3f km/s)", ship.speed(), ship.speed() / 1e3);
        DrawText(buf, 10, py, 16, RAYWHITE); py += 20;

        if (ship.elements.valid)
            sprintf(buf, "SMA: %.0f km   Ecc: %.5f", ship.elements.sma / 1e3, ship.elements.ecc);
        else
            sprintf(buf, "SMA: ---   Ecc: ---");
        DrawText(buf, 10, py, 16, RAYWHITE);

        EndDrawing();
    }

    for (auto* p : allBodies) delete p;
    CloseWindow();
    return 0;
}