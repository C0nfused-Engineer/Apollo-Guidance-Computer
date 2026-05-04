#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include "Planet.h"
#include "SystemLoader.h"
#include <cstdio>

int main(int argc, char* argv[]) {
    const int screenWidth  = 1280;
    const int screenHeight = 720;

    InitWindow(screenWidth, screenHeight, "Orbital Sim");
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

    // --- Camera state (all distances in meters) ---
    float   camYaw         = 0.0f;
    float   camPitch       = 80.0f;
    double  camDist        = 1.5e11;  // ~1 AU, good starting view of inner solar system
    int     camTargetIndex = 0;
    Planet* camTarget      = allBodies[0];

    Camera3D camera   = {};
    camera.fovy       = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    camera.up         = { 0.0f, 0.0f, 1.0f };

    // --- Time ---
    double simTime   = 0.0;
    double timeScale = 1.0;

    while (!WindowShouldClose()) {
        // --- Update ---
        double dt = GetFrameTime();
        simTime += dt * timeScale;
        root->update(simTime);

        if (IsKeyPressed(KEY_UP))   timeScale *= 10.0;
        if (IsKeyPressed(KEY_DOWN)) timeScale /= 10.0;

        if (IsKeyPressed(KEY_TAB)) {
            camTargetIndex = (camTargetIndex + 1) % (int)allBodies.size();
            camTarget      = allBodies[camTargetIndex];
        }

        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            Vector2 delta = GetMouseDelta();
            camYaw   += delta.x * 0.3f;
            camPitch += delta.y * 0.3f;
            camPitch  = Clamp(camPitch, 5.0f, 89.9f);
        }

        // Proportional zoom in meters
        camDist -= GetMouseWheelMove() * camDist * 0.1;
        if (camDist < 1.0) camDist = 1.0; // never closer than 1 meter

        // --- Compute camera world position (meters) ---
        float pitchRad = camPitch * DEG2RAD;
        float yawRad   = camYaw   * DEG2RAD;

        DVec3 targetPos = camTarget->position;

        DVec3 camWorldPos = {
            targetPos.x + camDist * cos(pitchRad) * sin(yawRad),
            targetPos.y + camDist * cos(pitchRad) * cos(yawRad),
            targetPos.z + camDist * sin(pitchRad)
        };

        // Camera rendered at float origin — all geometry shifted relative to camWorldPos
        camera.position = { 0.0f, 0.0f, 0.0f };
        camera.target   = targetPos.toVec3(camWorldPos);
        camera.up       = { 0.0f, 0.0f, 1.0f };

        // Clip planes scale with distance: near = 0.01% of dist, far = 1000x dist
        // This keeps depth buffer precision good at any zoom level.
        rlSetClipPlanes(camDist * 0.0001, camDist * 1000.0);

        // --- Draw ---
        BeginDrawing();
            ClearBackground(BLACK);

            BeginMode3D(camera);
                root->draw(camWorldPos);

                // Ecliptic reference grid in XY plane (Z = 0), in meters
                // Grid extends ~6 AU, spacing ~0.1 AU
                const double AU       = 1.496e11;
                const double gridExt  = 6.0  * AU;
                const double gridStep = 0.1  * AU;
                Color gridColor       = { 30, 30, 30, 255 };

                for (double i = -gridExt; i <= gridExt + gridStep * 0.5; i += gridStep) {
                    DVec3 a1 = { -gridExt, i,       0.0 };
                    DVec3 b1 = {  gridExt, i,       0.0 };
                    DVec3 a2 = { i,       -gridExt, 0.0 };
                    DVec3 b2 = { i,        gridExt, 0.0 };
                    DrawLine3D(a1.toVec3(camWorldPos), b1.toVec3(camWorldPos), gridColor);
                    DrawLine3D(a2.toVec3(camWorldPos), b2.toVec3(camWorldPos), gridColor);
                }
            EndMode3D();

            // --- HUD ---
            DrawText("UP/DOWN: time scale | TAB: next body | RMB drag: rotate | Scroll: zoom", 10, 10, 16, RAYWHITE);

            char buf[128];
            sprintf(buf, "Target: %s", camTarget->name.c_str());
            DrawText(buf, 10, 30, 18, YELLOW);

            // Time scale
            if      (timeScale >= 86400.0) sprintf(buf, "Time scale: %.1f days/sec",  timeScale / 86400.0);
            else if (timeScale >= 3600.0)  sprintf(buf, "Time scale: %.1f hrs/sec",   timeScale / 3600.0);
            else if (timeScale >= 60.0)    sprintf(buf, "Time scale: %.1f min/sec",   timeScale / 60.0);
            else                           sprintf(buf, "Time scale: %.0fx",           timeScale);
            DrawText(buf, 10, 50, 18, RAYWHITE);

            // Sim time
            long long totalSec = (long long)simTime;
            sprintf(buf, "Sim time: %dd %02dh %02dm %02ds",
                (int)(totalSec / 86400),
                (int)((totalSec % 86400) / 3600),
                (int)((totalSec % 3600)  / 60),
                (int)(totalSec % 60));
            DrawText(buf, 10, 70, 18, RAYWHITE);

            // Camera distance — auto-unit
            if      (camDist >= 1.496e11) sprintf(buf, "Cam dist: %.3f AU",   camDist / 1.496e11);
            else if (camDist >= 1e9)      sprintf(buf, "Cam dist: %.0f km",   camDist / 1e3);
            else if (camDist >= 1e3)      sprintf(buf, "Cam dist: %.0f m",    camDist);
            else                          sprintf(buf, "Cam dist: %.2f m",    camDist);
            DrawText(buf, 10, 90, 18, RAYWHITE);

        EndDrawing();
    }

    for (auto* p : allBodies) delete p;
    CloseWindow();
    return 0;
}