#pragma once

#include <cstddef>
#include <stdexcept>
#include "Node.hxx"
#include "Vec3D.hxx"

class Cell
{
private:
	Node m_N0, m_N1, m_N2, m_N3;
	
	double m_volume;
	double m_ddx_lamb0, m_ddy_lamb0, m_ddz_lamb0;
	double m_ddx_lamb1, m_ddy_lamb1, m_ddz_lamb1;
	double m_ddx_lamb2, m_ddy_lamb2, m_ddz_lamb2;
	double m_ddx_lamb3, m_ddy_lamb3, m_ddz_lamb3;

public:
	Cell();
	Cell(Node const& N0, Node const& N1, Node const& N2, Node const& N3);
	Cell(Cell const& cell);
	~Cell();

	void compute_nodes();

	Node const& get_node(size_t node_nb) const;
	bool contains(size_t node_nb) const;

	double get_volume() const;

	double get_ddx_lamb0() const;
	double get_ddx_lamb1() const;
	double get_ddx_lamb2() const;
	double get_ddx_lamb3() const;

	double get_ddy_lamb0() const;
	double get_ddy_lamb1() const;
	double get_ddy_lamb2() const;
	double get_ddy_lamb3() const;

	double get_ddz_lamb0() const;
	double get_ddz_lamb1() const;
	double get_ddz_lamb2() const;
	double get_ddz_lamb3() const;

	void display() const;
	
	Cell operator=(Cell const& cell);
};