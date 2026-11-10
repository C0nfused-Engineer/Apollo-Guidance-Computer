// constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace constants {
	// Core Constants
	constexpr double pi = 3.141592653589793;
	constexpr double pi2 = 2*pi;
	constexpr double au_meters = 149597870700.0; // m
	constexpr double G = 6.6743e-11; // m3/kgs2
	constexpr double c = 299792458.0; // m/s
	constexpr double sec_day = 24*3600.0; // s
	constexpr double sec_year = 31557600.0; // s
	constexpr double rad_to_deg = 180.0/pi; // radians -> degree
	constexpr double deg_to_rad = pi/180.0; // degree -> radians
	constexpr double J2000 = 2451545.0; // Days before Jan 1, 2000 12:00
	
	
	// Rocketry Constants
	constexpr double g0 = 9.80665; // m/s2
	constexpr double earth_radius = 6378137.0; // m
	constexpr double earth_omega = 7.2921159e-5; // radians/s
	
	// Atmosphere Constants
	constexpr double R_univ = 8.314462618; // J/molK
	constexpr double P0 = 101325.0; // Pa
	constexpr double T0 = 288.15; // K
	constexpr double rho0 = 1.225; // kg/m3
	
	// Environmental
	constexpr double stef_boltz = 5.670374419e-8; // W/m2K4
	constexpr double boltz = 1.380649e-23; // J/K
}

#endif