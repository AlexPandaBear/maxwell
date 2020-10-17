#include "Vec3D.hxx"

Vec3D::Vec3D() {}

Vec3D::Vec3D(double x, double y, double z) :
	m_x(x),
	m_y(y),
	m_z(z) {}

Vec3D::~Vec3D() {}

double Vec3D::get_x() const
{
	return m_x;
}

double Vec3D::get_y() const
{
	return m_y;
}

double Vec3D::get_z() const
{
	return m_z;
}

void Vec3D::set_x(double x)
{
	m_x = x;
}

void Vec3D::set_y(double y)
{
	m_y = y;
}

void Vec3D::set_z(double z)
{
	m_z = z;
}

void Vec3D::set(double x, double y, double z)
{
	m_x = x;
	m_y = y;
	m_z = z;
}

bool Vec3D::operator==(Vec3D const& v)
	{
		if (v.get_x() != m_x || v.get_y() != m_y || v.get_z() != m_z)
		{
			return false;
		}

		return true;
	}

bool Vec3D::operator!=(Vec3D const& v)
{
	return !(*this == v);
}