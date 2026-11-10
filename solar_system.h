// solar_system.h
#ifndef SOLAR_SYSTEM_H
#define SOLAR_SYSTEM_H

#include <unordered_map>
#include <string>
#include <vector>
#include <limits>
#include <stdexcept>
#include <utility>
#include "vec3.h"

constexpr int noParent = -1;

struct BodyInfo {
	std::string name;
	int naifId = -1;
	double mu = std::numeric_limits<double>::quiet_NaN(); // m3/s2
	double radius = std::numeric_limits<double>::quiet_NaN(); // mean, m
	int parentIndex = noParent;
};

struct BodyState {
	Vec3 position; // holds position in state, and velocity in derivitive
	Vec3 velocity; // holds velocity in state, and acceleration in derivitive
	
	BodyState operator+(const BodyState& b) const {
		return BodyState{position + b.position, velocity + b.velocity};
	}
	
	BodyState& operator+=(const BodyState& b) {
		position += b.position;
		velocity += b.velocity;
		return *this;
	}
	
	BodyState operator-(const BodyState& b) const {
		return BodyState{position - b.position, velocity - b.velocity};
	}
	
	BodyState& operator-=(const BodyState& b) {
		position -= b.position;
		velocity -= b.velocity;
		return *this;
	}

	
	BodyState operator*(double s) const {
		return BodyState{s * position, s * velocity};
	}
	
	BodyState& operator*=(double s) {
		position *= s;
		velocity *= s;
		return *this;
	}
};

inline BodyState operator*(double lhs, const BodyState& rhs) {
	return rhs * lhs;
}

inline std::vector<BodyState> operator*(double lhs, std::vector<BodyState> rhs) {
	for(auto& s : rhs) {
		s *= lhs;
	}
	
	return rhs;
}

inline std::vector<BodyState> operator*(std::vector<BodyState> lhs, double rhs) {
	return rhs * lhs;
}

inline std::vector<BodyState> operator+(std::vector<BodyState> lhs, const std::vector<BodyState>& rhs) {
	if(lhs.size() != rhs.size()) throw std::invalid_argument("state array sizes not compatible");
	
	for(size_t i = 0; i < lhs.size(); i++) {
		lhs[i] += rhs[i];
	}
	
	return lhs;
}

inline std::vector<BodyState> operator-(std::vector<BodyState> lhs, const std::vector<BodyState>& rhs) {
	if(lhs.size() != rhs.size()) throw std::invalid_argument("state array sizes not compatible");
	
	for(size_t i = 0; i < lhs.size(); i++) {
		lhs[i] -= rhs[i];
	}
	
	return lhs;
}


struct SolarSystem {
	private:
		std::vector<BodyInfo> bodyInfos;
		std::vector<BodyState> bodyStates;
		
		std::unordered_map<int, int> lookup;
	
	public:
		double epoch = 0.0; // seconds since J2000(TDB), s
		
		int addBody(std::string n, int nId, double mu, double r, int pId) {
			int idx = static_cast<int>(bodyInfos.size());
			if(lookup.find(nId) != lookup.end()) {
				throw std::invalid_argument("body with same naifId already exists");
			}
			
			if(pId != noParent && !(pId >= 0 && pId < idx)) {
				throw std::invalid_argument("parent body does not yet exist");
			}
			
			lookup.emplace(nId, idx);
			
			bodyInfos.emplace_back(BodyInfo{std::move(n), nId, mu, r, pId});
			bodyStates.emplace_back(BodyState{Vec3(), Vec3()});
			return idx;
		}
		
		int indexOf(int naifId) const {
			auto it = lookup.find(naifId);
			if(it == lookup.end()) throw std::invalid_argument("planet naifId: " + std::to_string(naifId) + " does not exist in lookup");
			
			return it->second;
		}
		
		size_t size() const {
			return bodyStates.size();
		}
		
		bool contains(int naifId) const {
			auto it = lookup.find(naifId);
			return it != lookup.end();
		}
		
		const BodyInfo& getInfo(int idx) const {
			return bodyInfos.at(idx);
		}
		
		const std::vector<BodyState>& getStates() const {
			return bodyStates;
		}
		
		void setStates(const std::vector<BodyState>& s, double deltaTime) {
			if(s.size() != bodyStates.size()) throw std::invalid_argument("size mismatch");
			bodyStates = s;
			epoch += deltaTime;
		}
		
		const BodyState& getState(int idx) const {
			return bodyStates.at(idx);
		}
		
		void setState(int idx, const BodyState& state) {
			bodyStates.at(idx) = state;
		}
};

#endif