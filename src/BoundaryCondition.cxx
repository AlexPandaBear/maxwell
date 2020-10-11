#include "BoundaryCondition.hxx"

BoundaryCondition::BoundaryCondition(size_t node_nb, double Ex, double Ey, double Bx, double By) :
	m_node_nb(node_nb),
	m_Ex(Ex),
	m_Ey(Ey),
	m_Bx(Bx),
	m_By(By) {}

BoundaryCondition::~BoundaryCondition() {}

size_t BoundaryCondition::get_node_nb() const
{
	return m_node_nb;
}

double BoundaryCondition::get_Ex() const
{
	return m_Ex;
}

double BoundaryCondition::get_Ey() const
{
	return m_Ey;
}

double BoundaryCondition::get_Bx() const
{
	return m_Bx;
}

double BoundaryCondition::get_By() const
{
	return m_By;
}