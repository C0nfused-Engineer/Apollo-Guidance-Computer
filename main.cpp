#include "raylib.h"
#include "raymath.h"
#include "Planet.h"
#include "SystemLoader.h"
#include <cstdio>

int main(int argc, char* argv[]) {
    const int screenWidth  = 1280;
    const int screenHeight = 720;

    InitWindow(screenWidth, screenHeight, "Orbital Sim");
    SetTargetFPS(60);

    // --- Load system from config ---
    const char* configPath = (argc > 1) ? argv[1] : "solar_system.cfg";

    SystemLoader loader;
    try {
        loader.load(configPath);
    } catch (const std::exception& e) {
        // Show error in window then exit
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
    loader.release(); // we now own the Planet pointers

    // --- Camera ---
    Camera3D camera   = {};
    camera.fovy       = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    float   camYaw         = 0.0f;
    float   camPitch       = 80.0f;
    float   camDist        = 150.0f;
    int     camTargetIndex = 0;
    Planet* camTarget      = allBodies[0];

    // --- Time ---
    double simTime   = 0.0;
    double timeScale = 1.0; // 1 real second = 1 sim second

    while (!WindowShouldClose()) {
        // --- Update ---
        double dt = GetFrameTime();
        simTime += dt * timeScale;
        root->update(simTime);

        // Time scale
        if (IsKeyPressed(KEY_UP))   timeScale *= 10.0;
        if (IsKeyPressed(KEY_DOWN)) timeScale /= 10.0;

        // TAB: cycle camera target
        if (IsKeyPressed(KEY_TAB)) {
            camTargetIndex = (camTargetIndex + 1) % (int)allBodies.size();
            camTarget = allBodies[camTargetIndex];
        }

        // RMB drag: rotate camera
        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            Vector2 delta = GetMouseDelta();
            camYaw   += delta.x * 0.3f;
            camPitch += delta.y * 0.3f;
            camPitch  = Clamp(camPitch, 5.0f, 89.9f);
        }

        // Scroll: zoom (proportional to current distance)
        camDist -= GetMouseWheelMove() * camDist * 0.1f;

        // Compute camera offset from yaw/pitch/dist
        float pitchRad = camPitch * DEG2RAD;
        float yawRad   = camYaw   * DEG2RAD;

        Vector3 offset = {
            camDist * cosf(pitchRad) * sinf(yawRad),
            camDist * cosf(pitchRad) * cosf(yawRad),
            camDist * sinf(pitchRad)
        };

        camera.target   = camTarget->position;
        camera.position = Vector3Add(camTarget->position, offset);
        camera.up       = { 0.0f, 0.0f, 1.0f };

        // --- Draw ---
        BeginDrawing();
            ClearBackground(BLACK);

            BeginMode3D(camera);
                root->draw();

                // Ecliptic reference grid in XY plane (Z=0)
                int   gridSize    = 500;
                int   gridSpacing = 10;
                Color gridColor   = { 30, 30, 30, 255 };

                for (int i = -gridSize; i <= gridSize; i += gridSpacing) {
                    DrawLine3D(
                        { (float)-gridSize, (float)i, 0.0f },
                        { (float) gridSize, (float)i, 0.0f },
                        gridColor
                    );
                    DrawLine3D(
                        { (float)i, (float)-gridSize, 0.0f },
                        { (float)i, (float) gridSize, 0.0f },
                        gridColor
                    );
                }
            EndMode3D();

            // HUD
            DrawText("UP/DOWN: time scale | TAB: next body | RMB drag: rotate | Scroll: zoom", 10, 10, 16, RAYWHITE);

            char buf[128];
            sprintf(buf, "Target: %s", camTarget->name.c_str());
            DrawText(buf, 10, 30, 18, YELLOW);

            // Time scale — auto-unit display
            if (timeScale >= 86400.0)
                sprintf(buf, "Time scale: %.1f days/sec", timeScale / 86400.0);
            else if (timeScale >= 3600.0)
                sprintf(buf, "Time scale: %.1f hrs/sec",  timeScale / 3600.0);
            else if (timeScale >= 60.0)
                sprintf(buf, "Time scale: %.1f min/sec",  timeScale / 60.0);
            else
                sprintf(buf, "Time scale: %.0fx", timeScale);
            DrawText(buf, 10, 50, 18, RAYWHITE);

            // Sim time
            long long totalSec = (long long)simTime;
            int simDays = totalSec / 86400;
            int simHrs  = (totalSec % 86400) / 3600;
            int simMins = (totalSec % 3600)  / 60;
            int simSecs = totalSec % 60;
            sprintf(buf, "Sim time: %dd %02dh %02dm %02ds", simDays, simHrs, simMins, simSecs);
            DrawText(buf, 10, 70, 18, RAYWHITE);

            // Camera distance
            sprintf(buf, "Cam dist: %.2f units", camDist);
            DrawText(buf, 10, 90, 18, RAYWHITE);

        EndDrawing();
    }

    // Clean up
    for (auto* p : allBodies)
        delete p;

    CloseWindow();
    return 0;
}