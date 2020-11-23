#include "MatrixFEM.hxx"

MatrixFEM::MatrixFEM(size_t nb_nodes, size_t nb_scalar_by_node) :
	Matrix(nb_nodes*nb_scalar_by_node),
	m_nb_nodes(nb_nodes),
	m_nb_scalar_by_node(nb_scalar_by_node) {}

MatrixFEM::~MatrixFEM() {}

void MatrixFEM::set_interaction_bloc(size_t node_i, size_t node_j, Matrix const& bloc)
{
	std::pair<size_t, size_t> bloc_shape(bloc.get_shape());

	if (bloc_shape.first != m_nb_scalar_by_node or bloc_shape.second != m_nb_scalar_by_node)
	{
		throw std::invalid_argument("Node to node interaction bloc shape does not agree with expected shape");
	}

	for (size_t i = 0; i < m_nb_scalar_by_node; i++)
	{
		for (size_t j = 0; j < m_nb_scalar_by_node; j++)
		{
			set_coef(m_nb_scalar_by_node*node_i + i, m_nb_scalar_by_node*node_j + j, bloc.get_coef(i, j));
		}
	}
}

void MatrixFEM::add_to_interaction_bloc(size_t node_i, size_t node_j, Matrix const& bloc)
{
	std::pair<size_t, size_t> bloc_shape(bloc.get_shape());

	if (bloc_shape.first != m_nb_scalar_by_node or bloc_shape.second != m_nb_scalar_by_node)
	{
		throw std::invalid_argument("Node to node interaction bloc shape does not agree with expected shape");
	}

	for (size_t i = 0; i < m_nb_scalar_by_node; i++)
	{
		for (size_t j = 0; j < m_nb_scalar_by_node; j++)
		{
			add_to_coef(m_nb_scalar_by_node*node_i + i, m_nb_scalar_by_node*node_j + j, bloc.get_coef(i, j));
		}
	}
}