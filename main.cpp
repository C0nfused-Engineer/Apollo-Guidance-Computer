// Apollo Guidance Computer: Aerospace Engineering Simulation
#include <iostream>
#include "constants.h"
#include "vec3.h"
#include "solar_system.h"
#include "gravity.h"
#include "diagnostics.h"

#include <cmath>
#include <cstdio>

int main(int /*argc*/, char* /*argv*/[]) {
	// ---- Hard-coded Sun / Earth / Moon test system (SI units) ----
	// Simplified circular, coplanar initial conditions at t = 0 (NOT Horizons data):
	//   Sun at the origin, Earth 1 AU out on +x, Moon 384,400 km beyond Earth on +x,
	//   each moving in +y at its circular-orbit speed. Sun's velocity is set so total
	//   momentum is zero (barycentric frame).
	const double muSun   = 1.32712440018e20;   // m^3/s^2
	const double muEarth = 3.986004418e14;
	const double muMoon  = 4.9048695e12;

	SolarSystem sys;
	int sun   = sys.addBody("Sun",   10,  muSun,   6.957e8,  noParent);
	int earth = sys.addBody("Earth", 399, muEarth, 6.371e6,  sun);
	int moon  = sys.addBody("Moon",  301, muMoon,  1.7374e6, earth);

	const double rEarth = constants::au_meters;      // Sun-Earth distance
	const double rMoon  = 3.844e8;                   // Earth-Moon distance

	double vEarth = std::sqrt(muSun / rEarth);       // circular speed about the Sun
	double vMoon  = vEarth + std::sqrt(muEarth / rMoon);  // Earth's velocity + circular speed about Earth

	sys.setState(sun,   BodyState{Vec3(0.0, 0.0, 0.0),             Vec3(0.0, 0.0, 0.0)});
	sys.setState(earth, BodyState{Vec3(rEarth, 0.0, 0.0),          Vec3(0.0, vEarth, 0.0)});
	sys.setState(moon,  BodyState{Vec3(rEarth + rMoon, 0.0, 0.0),  Vec3(0.0, vMoon, 0.0)});

	// Zero the total momentum: give the Sun the opposite of everyone else's mu*v.
	Vec3 p = muEarth * sys.getState(earth).velocity + muMoon * sys.getState(moon).velocity;
	sys.setState(sun, BodyState{sys.getState(sun).position, Vec3(0.0, 0.0, 0.0) - p / muSun});

	// ---- Print the initial state and the accelerations at t = 0 ----
	std::vector<BodyState> deriv = gravitation(sys, sys.getStates());

	std::printf("epoch = %.1f s\n\n", sys.epoch);
	for(size_t i = 0; i < sys.size(); ++i) {
		const BodyInfo& info = sys.getInfo(static_cast<int>(i));
		const BodyState& st  = sys.getState(static_cast<int>(i));
		std::printf("%-6s r = (%.6e, %.6e, %.6e) m\n", info.name.c_str(),
		            st.position.x(), st.position.y(), st.position.z());
		std::printf("       v = (%.6e, %.6e, %.6e) m/s\n",
		            st.velocity.x(), st.velocity.y(), st.velocity.z());
		std::printf("       a = (%.6e, %.6e, %.6e) m/s^2   |a| = %.6e\n\n",
		            deriv[i].velocity.x(), deriv[i].velocity.y(), deriv[i].velocity.z(),
		            deriv[i].velocity.mag());
	}
	
	double E0 = totalEnergy(sys, sys.getStates());
	double P0 = totalLinearMomentum(sys, sys.getStates()).mag();
	double L0 = totalAngularMomentum(sys, sys.getStates()).mag();
	
	std::printf("Total Energy: %.6e J, Total Momentum: %.6e Ns(Linear), %.6e Ns(Angular)\n", E0, P0, L0);
	
	return 0;
}