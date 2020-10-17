#include "BoundaryCondition.hxx"

BoundaryCondition::BoundaryCondition(size_t node_nb, Vec3D E, Vec3D B) :
	m_node_nb(node_nb),
	m_E(E),
	m_B(B) {}

BoundaryCondition::~BoundaryCondition() {}

size_t BoundaryCondition::get_node_nb() const
{
	return m_node_nb;
}

Vec3D BoundaryCondition::get_E() const
{
	return m_E;
}

Vec3D BoundaryCondition::get_B() const
{
	return m_B;
}