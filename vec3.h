// vector3
#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <string>

struct Vec3 {
	public:
		// Constructor
		Vec3(double x, double y, double z) : _x(x), _y(y), _z(z) {}
		Vec3() : _x(0.0f), _y(0.0f), _z(0.0f) {}
		
		// Direct Access
		double x() const { return _x; }
		double y() const { return _y; }
		double z() const { return _z; }
		
		// Vector Functions
		double mag() const {
			return sqrt((*this).dot(*this));
		}
		
		Vec3 normal() const {
			return (*this)/mag();
		}
		
		double dot(const Vec3& b) const {
			return _x*b._x + _y*b._y + _z*b._z;
		}
		
		Vec3 cross(const Vec3& b) const {
			return Vec3(_y*b._z - _z*b._y, _z*b._x - _x*b._z, _x*b._y - _y*b._x);
		}
	
		// Operator Overloading
		Vec3 operator+(const Vec3& rhs) const {
			return Vec3(_x + rhs._x, _y + rhs._y, _z + rhs._z);
		}
		
		Vec3& operator+=(const Vec3& rhs) {
			_x += rhs._x;
			_y += rhs._y;
			_z += rhs._z;
			return *this;
		}
		
		Vec3 operator-() const {
			return Vec3(-_x, -_y, -_z);
		}
		
		Vec3 operator-(const Vec3& rhs) const {
			return Vec3(_x - rhs._x, _y - rhs._y, _z - rhs._z);
		}
		
		Vec3& operator-=(const Vec3& rhs) {
			_x -= rhs._x;
			_y -= rhs._y;
			_z -= rhs._z;
			return *this;
		}
		
		Vec3 operator*(double rhs) const {
			return Vec3(_x*rhs, _y*rhs, _z*rhs);
		}
		
		Vec3& operator*=(double rhs) {
			_x *= rhs;
			_y *= rhs;
			_z *= rhs;
			return *this;
		}
		
		Vec3 operator/(double rhs) const {
			return Vec3(_x/rhs, _y/rhs, _z/rhs);
		}
		
		Vec3& operator/=(double rhs) {
			_x /= rhs;
			_y /= rhs;
			_z /= rhs;
			return *this;
		}
	
		// Utility
		std::string toString() const {
			return "[ " + std::to_string(_x) + ", " + std::to_string(_y) + ", " + std::to_string(_z) + "]";
		}
	
	private:
		double _x, _y, _z;
};

inline Vec3 operator*(double lhs, const Vec3& rhs) {
	return rhs * lhs;
}

#endif