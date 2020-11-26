#include "Vec3D.hxx"

Vec3D::Vec3D() :
	m_view(false),
	m_x(0.),
	m_y(0.),
	m_z(0.),
	m_ptr_x(NULL),
	m_ptr_y(NULL),
	m_ptr_z(NULL) {}

Vec3D::Vec3D(double x, double y, double z) :
	m_view(false),
	m_x(x),
	m_y(y),
	m_z(z),
	m_ptr_x(NULL),
	m_ptr_y(NULL),
	m_ptr_z(NULL) {}

Vec3D::Vec3D(double* ptr_x, double* ptr_y, double* ptr_z) :
	m_view(true),
	m_x(0.),
	m_y(0.),
	m_z(0.),
	m_ptr_x(ptr_x),
	m_ptr_y(ptr_y),
	m_ptr_z(ptr_z) {}

Vec3D::Vec3D(Vec3D const& v) :
	m_view(false),
	m_x(v.get_x()),
	m_y(v.get_y()),
	m_z(v.get_z()),
	m_ptr_x(NULL),
	m_ptr_y(NULL),
	m_ptr_z(NULL) {}

Vec3D::~Vec3D() {}

double Vec3D::get_x() const
{
	if (m_view)
	{
		return *m_ptr_x;
	}

	return m_x;
}

double Vec3D::get_y() const
{
	if (m_view)
	{
		return *m_ptr_y;
	}

	return m_y;
}

double Vec3D::get_z() const
{
	if (m_view)
	{
		return *m_ptr_z;
	}

	return m_z;
}

void Vec3D::set_x(double x)
{
	if (m_view)
	{
		throw std::logic_error("Vec3D are read-only objects in view mode");
	}
	
	m_x = x;
}

void Vec3D::set_y(double y)
{
	if (m_view)
	{
		throw std::logic_error("Vec3D are read-only objects in view mode");
	}
		
	m_y = y;
}

void Vec3D::set_z(double z)
{
	if (m_view)
	{
		throw std::logic_error("Vec3D are read-only objects in view mode");
	}
	
	m_z = z;
}

void Vec3D::set(double x, double y, double z)
{
	set_x(x);
	set_y(y);
	set_z(z);
}

double Vec3D::compute_squared_norm() const
{
	return pow(get_x(), 2) + pow(get_y(), 2) + pow(get_z(), 2);
}

double Vec3D::compute_norm() const
{
	return sqrt(compute_squared_norm());
}

void Vec3D::normalize()
{
	(*this) /= compute_norm();
}

void Vec3D::display() const
{
	std::cout << "[ " << m_x << " " << m_y << " " << m_z << " ]^T" << std::endl;
}

bool Vec3D::operator==(Vec3D const& v) const
{
	if (v.get_x() != get_x() || v.get_y() != get_y() || v.get_z() != get_z())
	{
		return false;
	}

	return true;
}

bool Vec3D::operator!=(Vec3D const& v) const
{
	return !(*this == v);
}

void Vec3D::operator+=(Vec3D const& v)
{
	if (m_view)
	{
		throw std::logic_error("Vec3D are read-only objects in view mode");
	}

	m_x += v.get_x();
	m_y += v.get_y();
	m_z += v.get_z();
}

void Vec3D::operator-=(Vec3D const& v)
{
	if (m_view)
	{
		throw std::logic_error("Vec3D are read-only objects in view mode");
	}

	m_x -= v.get_x();
	m_y -= v.get_y();
	m_z -= v.get_z();
}

void Vec3D::operator*=(double a)
{
	if (m_view)
	{
		throw std::logic_error("Vec3D are read-only objects in view mode");
	}

	m_x *= a;
	m_y *= a;
	m_z *= a;
}

void Vec3D::operator/=(double a)
{
	if (m_view)
	{
		throw std::logic_error("Vec3D are read-only objects in view mode");
	}

	if (a == 0.)
	{
		throw std::invalid_argument("Cannot divide vector by zero");
	}

	m_x /= a;
	m_y /= a;
	m_z /= a;
}

Vec3D Vec3D::operator+(Vec3D const& v) const
{
	Vec3D result(*this);
	result += v;
	return result;
}

Vec3D Vec3D::operator-(Vec3D const& v) const
{
	Vec3D result(*this);
	result -= v;
	return result;
}

Vec3D Vec3D::operator*(double a) const
{
	Vec3D result(*this);
	result *= a;
	return result;
}

Vec3D Vec3D::operator/(double a) const
{
	Vec3D result(*this);
	result /= a;
	return result;
}

double Vec3D::dot_product(Vec3D const& v1, Vec3D const& v2)
{
	return v1.get_x()*v2.get_x() + v1.get_y()*v2.get_y() + v1.get_z()*v2.get_z();
}

Vec3D Vec3D::cross_product(Vec3D const& v1, Vec3D const& v2)
{
	return Vec3D(v1.get_y()*v2.get_z() - v1.get_z()*v2.get_y(), v1.get_z()*v2.get_x() - v1.get_x()*v2.get_z(), v1.get_x()*v2.get_y() - v1.get_y()*v2.get_x());
}