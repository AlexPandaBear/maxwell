#pragma once

class BoundaryCondition
{
private:
	const size_t m_node_nb;
	const double m_Ex, m_Ey;
	const double m_Bx, m_By;

public:
	BoundaryCondition(size_t node_nb, double Ex, double Ey, double Bx, double By);
	~BoundaryCondition();
	
	size_t get_node_nb() const;
	double get_Ex() const;
	double get_Ey() const;
	double get_Bx() const;
	double get_By() const;
};