// gravity.h
#ifndef GRAVITY_H
#define GRAVITY_H

#include <cassert>
#include <stdexcept>
#include "solar_system.h"

std::vector<BodyState> gravitation(const SolarSystem& s, const std::vector<BodyState>& state) {
	if(state.size() != s.size()) throw std::invalid_argument("state array size mismatch to solar system");
	std::vector<BodyState> derivState;
	derivState.reserve(state.size());
	
	for(size_t i = 0; i < state.size(); i++) {
		Vec3 vel = state[i].velocity;
		Vec3 acc = Vec3();
		for(size_t j = 0; j < state.size(); j++) {
			if(j == i) continue;
			
			double mu_j = s.getInfo(static_cast<int>(j)).mu;
			Vec3 r = state[j].position - state[i].position;
			double r2 = r.dot(r);
			assert(r2 > 0 && "Zero Distance");
			acc += mu_j * r/(r2 * sqrt(r2));
		}
		
		derivState.emplace_back(BodyState{vel, acc});
	}
	
	return derivState;
}

#endif