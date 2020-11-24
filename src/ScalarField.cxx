#include "ScalarField.hxx"

ScalarField::ScalarField(size_t nb_nodes) :
	m_nb_nodes(nb_nodes),
	m_values(std::vector<double>(nb_nodes, 0.)) {}

ScalarField::ScalarField(size_t nb_nodes, double value) :
	m_nb_nodes(nb_nodes),
	m_values(std::vector<double>(nb_nodes, value)) {}

ScalarField::~ScalarField() {}

void ScalarField::reset_nb_nodes(size_t nb_nodes)
	{
		m_nb_nodes = nb_nodes;
		m_values = std::vector<double>(nb_nodes);
	}

size_t ScalarField::get_nb_nodes() const
	{
		return m_nb_nodes;
	}

void ScalarField::set_value(size_t node_nb, double value)
{
	m_values[node_nb] = value;
}

double ScalarField::get_value(size_t node_nb) const
{
	return m_values[node_nb];
}

double* ScalarField::get_ptr()
{
	return &(m_values[0]);
}

size_t ScalarField::get_element_size() const
{
	return sizeof(double);
}