#include "Node.hxx"

Node::Node() {}

Node::Node(size_t id, double x, double y, double z, bool boundary) :
	m_id(id),
	m_x(x),
	m_y(y),
	m_z(z),
	m_boundary(boundary) {}

Node::~Node() {}

size_t Node::get_id() const
{
	return m_id;
}

double Node::get_x() const
{
	return m_x;
}

double Node::get_y() const
{
	return m_y;
}

double Node::get_z() const
{
	return m_z;
}

bool Node::get_boudary() const
{
	return m_boundary;
}

void Node::set_x(double x)
{
	m_x = x;
}

void Node::set_y(double y)
{
	m_y = y;
}

void Node::set_z(double z)
{
	m_z = z;
}

void Node::set_boundary(bool boundary)
{
	m_boundary = boundary;
}