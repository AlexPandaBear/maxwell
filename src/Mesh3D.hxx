#pragma once

#include <list>
#include <string>
#include <vector>
#include <iostream>
#include "Node.hxx"
#include "Cell.hxx"

class Mesh3D
{
private:
	size_t m_nb_nodes;
	size_t m_nb_cells;

	std::vector<Node> m_nodes;
	std::vector<Cell> m_cells;

public:
	Mesh3D();
	Mesh3D(std::string file);
	~Mesh3D();

	void generate_grid_mesh(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, size_t nx, size_t ny, size_t nz);

	size_t get_nb_nodes() const;
	size_t get_nb_cells() const;

	Node get_node(size_t node_nb) const;
	Cell const& get_cell(size_t cell_nb) const;
	const Cell* get_cell_ptr(size_t cell_nb) const;
	std::vector<size_t> get_neighbor_cells_id(size_t node_nb) const;
	
	void save(std::string file) const;

	Vec3D get_node_xyz(size_t node_nb) const;
};