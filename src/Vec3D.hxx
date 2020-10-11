#pragma once

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
};