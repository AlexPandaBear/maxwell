#pragma once

#include <cmath>

class Vec3D
{
private:
	double m_x;
	double m_y;
	double m_z;

public:
	Vec3D();
	Vec3D(double x, double y, double z);
	~Vec3D();

	double get_x() const;
	double get_y() const;
	double get_z() const;

	void set_x(double x);
	void set_y(double y);
	void set_z(double z);
	void set(double x, double y, double z);

	double compute_squared_norm() const;
	double compute_norm() const;
	void normalize();

	bool operator==(Vec3D const& v);
	bool operator!=(Vec3D const& v);

	Vec3D operator+(Vec3D const& v) const;
	Vec3D operator-(Vec3D const& v) const;
	Vec3D operator*(double a) const;
	Vec3D operator/(double a) const;

	void operator+=(Vec3D const& v);
	void operator-=(Vec3D const& v);
	void operator*=(double a);
	void operator/=(double a);

	static double dot_product(Vec3D const& v1, Vec3D const& v2);
	static Vec3D cross_product(Vec3D const& v1, Vec3D const& v2);
};