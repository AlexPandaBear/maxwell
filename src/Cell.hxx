#pragma once

#include <cstddef>
#include <stdexcept>
#include "VectorField.hxx"

class Cell
{
private:
	VectorField const& m_nodes_ref;
	size_t m_N0_id, m_N1_id, m_N2_id, m_N3_id;
	
	double m_volume;
	double m_ddx_lamb0, m_ddy_lamb0, m_ddz_lamb0;
	double m_ddx_lamb1, m_ddy_lamb1, m_ddz_lamb1;
	double m_ddx_lamb2, m_ddy_lamb2, m_ddz_lamb2;
	double m_ddx_lamb3, m_ddy_lamb3, m_ddz_lamb3;

public:
	Cell(VectorField const& nodes);
	Cell(VectorField const& nodes, size_t N0_id, size_t N1_id, size_t N2_id, size_t N3_id);
	Cell(Cell const& cell);
	~Cell();

	void compute_nodes();

	void set_global_node_id(size_t local_node_id, size_t global_node_id);
	void set_nodes_id(size_t N0_id, size_t N1_id, size_t N2_id, size_t N3_id);
	size_t get_global_node_id(size_t local_node_id) const;
	Vec3D get_node_xyz(size_t local_node_id) const;

	//void set_nodes_ref(VectorField const& nodes);
	VectorField const& get_nodes_ref() const;
	
	bool contains(size_t global_node_id) const;

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

	//void display() const;
	
	Cell operator=(Cell const& cell);
};