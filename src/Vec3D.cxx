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

double Vec3D::compute_squared_norm() const
{
	return m_x*m_x + m_y*m_y + m_z*m_z;
}

double Vec3D::compute_norm() const
{
	return sqrt(compute_squared_norm());
}

void Vec3D::normalize()
{
	double norm(compute_norm());
	m_x /= norm;
	m_y /= norm;
	m_z /= norm;
}

Vec3D Vec3D::operator+(Vec3D const& v) const
{
	return Vec3D(m_x + v.get_x(), m_y + v.get_y(), m_z + v.get_z());
}

Vec3D Vec3D::operator-(Vec3D const& v) const
{
	return Vec3D(m_x - v.get_x(), m_y - v.get_y(), m_z - v.get_z());
}

Vec3D Vec3D::operator*(double a) const
{
	return Vec3D(a*m_x, a*m_y, a*m_z);
}

Vec3D Vec3D::operator/(double a) const
{
	return Vec3D(m_x/a, m_y/a, m_y/a);
}

void Vec3D::operator+=(Vec3D const& v)
{
	m_x += v.get_x();
	m_y += v.get_y();
	m_z += v.get_z();
}

void Vec3D::operator-=(Vec3D const& v)
{
	m_x -= v.get_x();
	m_y -= v.get_y();
	m_z -= v.get_z();
}

void Vec3D::operator*=(double a)
{
	m_x *= a;
	m_y *= a;
	m_z *= a;
}

void Vec3D::operator/=(double a)
{
	m_x /= a;
	m_y /= a;
	m_z /= a;
}

double Vec3D::dot_product(Vec3D const& v1, Vec3D const& v2)
{
	return v1.get_x()*v2.get_x() + v1.get_y()*v2.get_y() + v1.get_z()*v2.get_z();
}

Vec3D Vec3D::cross_product(Vec3D const& v1, Vec3D const& v2)
{
	return Vec3D(v1.get_y()*v2.get_z() - v1.get_z()*v2.get_y(), v1.get_z()*v2.get_x() - v1.get_x()*v2.get_z(), v1.get_x()*v2.get_y() - v1.get_y()*v2.get_x());
}