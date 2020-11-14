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
	Vec3D m_S_012, m_S_013, m_S_023, m_S_123;

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
	Vec3D get_surface(size_t opposite_node_local_id) const;
	
	Cell operator=(Cell const& cell);
};