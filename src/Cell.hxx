#pragma once

#include <cstddef>
#include <stdexcept>
#include "VectorField.hxx"

class Cell
{
private:
	const double m_a = 0.13819660112;
	const double m_b = 0.58541019662;
	const double m_w = 0.04166666666;

	VectorField const& m_nodes_ref;
	size_t m_N0_id, m_N1_id, m_N2_id, m_N3_id;
	
	double m_volume;
	Vec3D m_S_012, m_S_013, m_S_023, m_S_123;
	bool m_boundary_0, m_boundary_1, m_boundary_2, m_boundary_3;

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

	void set_boundary_face(size_t opposite_node_local_id);
	bool is_boundary_cell() const;

	VectorField const& get_nodes_ref() const;
	
	bool contains(size_t global_node_id) const;

	double get_volume() const;
	Vec3D get_surface(size_t opposite_node_local_id) const;

	double compute_int3d_phi_psi(size_t i, size_t j) const;
	double compute_int3d_grad_phi_grad_psi(size_t i, size_t j) const;
	Vec3D compute_int3d_phi_grad_psi(size_t i) const;
	Vec3D compute_int2d_phi_psi_n(size_t i, size_t j) const;
	
	Cell operator=(Cell const& cell);
};