#pragma once

#include "Matrix.hxx"

class MatrixFEM : public Matrix
{
private:
	size_t m_nb_nodes;
	size_t m_nb_scalar_by_node;

public:
	MatrixFEM(size_t nb_nodes, size_t nb_scalar_by_node);
	~MatrixFEM();

	void set_interaction_bloc(size_t node_i, size_t node_j, Matrix const& bloc);
	void add_to_interaction_bloc(size_t node_i, size_t node_j, Matrix const& bloc);
};