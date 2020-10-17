#pragma once

#include <cstddef>
#include "Vec3D.hxx"

class BoundaryCondition
{
private:
	const size_t m_node_nb;
	const Vec3D m_E, m_B;

public:
	BoundaryCondition(size_t node_nb, Vec3D E, Vec3D B);
	~BoundaryCondition();
	
	size_t get_node_nb() const;
	Vec3D get_E() const;
	Vec3D get_B() const;
};