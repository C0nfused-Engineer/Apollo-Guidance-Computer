// diagnostics.h
#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include <stdexcept>
#include <vector>
#include <cmath>
#include "solar_system.h"
#include "constants.h"
#include "vec3.h"

double totalEnergy(const SolarSystem& s, const std::vector<BodyState>& states) {
	if(states.size() != s.size()) throw std::invalid_argument("state vector size mismatch");
	
	double KE = 0.0; // KE in Joule*G
	double PE = 0.0; // PE in Joule*G
	for(size_t i = 0; i < states.size(); i++) {
		double mu_i = s.getInfo(static_cast<int>(i)).mu;
		double v2 = (states[i].velocity).dot(states[i].velocity);
		KE += 0.5 * mu_i * v2;
		
		for(size_t j = i + 1; j < states.size(); j++) {
			double mu_j = s.getInfo(static_cast<int>(j)).mu;
			
			double r_ij = (states[j].position - states[i].position).mag();
			PE += -mu_i*mu_j/r_ij;
		}
	}
	
	return (KE + PE) / constants::G;
}

Vec3 totalLinearMomentum(const SolarSystem& s, const std::vector<BodyState>& states) {
	if(states.size() != s.size()) throw std::invalid_argument("state vector size mismatch");
	
	Vec3 P_total = Vec3(); // linear momentum in Ns*G
	for(size_t i = 0; i < states.size(); i++) {
		double mu = s.getInfo(static_cast<int>(i)).mu;
		P_total += mu * states[i].velocity;
	}
	
	return P_total / constants::G;
}

Vec3 totalAngularMomentum(const SolarSystem& s, const std::vector<BodyState>& states) {
	if(states.size() != s.size()) throw std::invalid_argument("state vector size mismatch");
	
	Vec3 L_total = Vec3(); // angular momentum in Ns*G
	for(size_t i = 0; i < states.size(); i++) {
		double mu = s.getInfo(static_cast<int>(i)).mu;
		L_total += mu * (states[i].position.cross(states[i].velocity));
	}
	
	return L_total / constants::G;
}

double energyError(double E, double E0) {
	return (E - E0) / abs(E0);
}

double angMomentError(const Vec3& L, const Vec3& L0) {
	return (L - L0).mag() / L0.mag();
}

double linMomentError(const Vec3& P, const Vec3& P0, const SolarSystem& s, const std::vector<BodyState>& states) {
	if(states.size() != s.size()) throw std::invalid_argument("state vector size mismatch");
	
	double numerator = (P - P0).mag();
	double denom = 0;
	for(size_t i = 0; i < states.size(); i++) {
		double mu = s.getInfo(static_cast<int>(i)).mu;
		denom += mu * states[i].velocity.mag();
	}
	
	return numerator / (denom / constants::G);
}

#endif