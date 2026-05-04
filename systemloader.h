#pragma once
#include "Planet.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Config file format
// ---------------------------------------------------------------------------
// Lines starting with '#' or empty lines are ignored.
//
// Global setting (must appear before any body blocks):
//   scale <value>          render_scale applied to all bodies (meters -> render units)
//
// Each body is defined as a block:
//   body <name>
//     parent        <name | none>
//     mass          <kg>
//     radius        <render units>
//     sma           <meters>          semi-major axis
//     ecc           <0..1>            eccentricity
//     inc           <degrees>         inclination
//     aop           <degrees>         argument of periapsis
//     lan           <degrees>         longitude of ascending node
//     mae           <degrees>         mean anomaly at epoch
//     color         <r> <g> <b>       0-255 each
//   end
// ---------------------------------------------------------------------------

struct SystemLoader {

    // Loaded bodies — owned by this struct. Call release() to take ownership.
    std::vector<Planet*> bodies;

    // The root body (no parent)
    Planet* root = nullptr;

    // Flat list in file order, useful for TAB cycling
    std::vector<Planet*> allBodies;

    ~SystemLoader() {
        for (auto* p : bodies)
            delete p;
    }

    // Transfer ownership to caller — bodies will NOT be deleted by destructor
    void release() { bodies.clear(); }

    // Load from file. Throws std::runtime_error on parse failure.
    void load(const std::string& path) {
        std::ifstream f(path);
        if (!f.is_open())
            throw std::runtime_error("Cannot open config file: " + path);

        double globalScale = 1.0;

        // Raw data per body before we can resolve parent pointers
        struct BodyDef {
            std::string name;
            std::string parentName;
            double mass          = 1.0;
            float  radius        = 1.0f;
            double sma           = 0.0;
            double ecc           = 0.0;
            double inc           = 0.0; // stored in radians after parse
            double aop           = 0.0;
            double lan           = 0.0;
            double mae           = 0.0;
            Color  color         = WHITE;
        };

        std::vector<BodyDef> defs;
        BodyDef current;
        bool inBody = false;

        std::string line;
        int lineNum = 0;

        auto trim = [](std::string s) -> std::string {
            size_t a = s.find_first_not_of(" \t\r\n");
            size_t b = s.find_last_not_of(" \t\r\n");
            return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
        };

        while (std::getline(f, line)) {
            ++lineNum;
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;

            std::istringstream ss(line);
            std::string token;
            ss >> token;

            // --- global scale ---
            if (token == "scale") {
                if (inBody) throw std::runtime_error("'scale' inside body block at line " + std::to_string(lineNum));
                ss >> globalScale;
                continue;
            }

            // --- body start ---
            if (token == "body") {
                if (inBody) throw std::runtime_error("Nested 'body' at line " + std::to_string(lineNum));
                current = BodyDef{};
                ss >> current.name;
                if (current.name.empty()) throw std::runtime_error("'body' missing name at line " + std::to_string(lineNum));
                inBody = true;
                continue;
            }

            // --- body end ---
            if (token == "end") {
                if (!inBody) throw std::runtime_error("'end' without 'body' at line " + std::to_string(lineNum));
                defs.push_back(current);
                inBody = false;
                continue;
            }

            if (!inBody) throw std::runtime_error("Unexpected token '" + token + "' at line " + std::to_string(lineNum));

            // --- body fields ---
            if (token == "parent") {
                ss >> current.parentName;
            } else if (token == "mass") {
                ss >> current.mass;
            } else if (token == "radius") {
                ss >> current.radius;
            } else if (token == "sma") {
                ss >> current.sma;
            } else if (token == "ecc") {
                ss >> current.ecc;
            } else if (token == "inc") {
                ss >> current.inc;
                current.inc *= DEG2RAD;
            } else if (token == "aop") {
                ss >> current.aop;
                current.aop *= DEG2RAD;
            } else if (token == "lan") {
                ss >> current.lan;
                current.lan *= DEG2RAD;
            } else if (token == "mae") {
                ss >> current.mae;
                current.mae *= DEG2RAD;
            } else if (token == "color") {
                int r, g, b;
                ss >> r >> g >> b;
                current.color = { (unsigned char)r, (unsigned char)g, (unsigned char)b, 255 };
            } else {
                throw std::runtime_error("Unknown field '" + token + "' at line " + std::to_string(lineNum));
            }
        }

        if (inBody)
            throw std::runtime_error("Unclosed 'body' block at end of file");

        // --- Build Planet objects ---
        // Map name -> Planet* for parent resolution
        std::map<std::string, Planet*> byName;

        for (auto& d : defs) {
            Planet* parent = nullptr;
            if (d.parentName != "none" && !d.parentName.empty()) {
                auto it = byName.find(d.parentName);
                if (it == byName.end())
                    throw std::runtime_error("Body '" + d.name + "' references unknown parent '" + d.parentName + "' (parent must be defined first)");
                parent = it->second;
            }

            Planet* p = new Planet(
                d.name,
                d.mass,
                (float)(d.radius * globalScale), // radius in meters -> render units
                d.sma,
                d.ecc,
                d.inc,
                d.aop,
                d.lan,
                d.mae,
                parent,
                globalScale
            );
            p->color = d.color;

            if (!parent) {
                if (root) throw std::runtime_error("Multiple root bodies (no parent) found; only one is allowed");
                root = p;
                root->position = { 0.0f, 0.0f, 0.0f };
            }

            bodies.push_back(p);
            allBodies.push_back(p);
            byName[d.name] = p;
        }

        if (!root)
            throw std::runtime_error("No root body (body with parent=none) found in config");
    }
};