#include "VectorField.hxx"

VectorField::VectorField(size_t nb_nodes) :
	m_nb_nodes(nb_nodes),
	m_values(std::vector<double>(3*nb_nodes)) {}

VectorField::~VectorField() {}

void VectorField::reset_nb_nodes(size_t nb_nodes)
	{
		m_nb_nodes = nb_nodes;
		m_values = std::vector<double>(3*nb_nodes);
	}

size_t VectorField::get_nb_nodes() const
	{
		return m_nb_nodes;
	}

void VectorField::set_value(size_t node_nb, double x_value, double y_value, double z_value)
{
	m_values[3*node_nb] = x_value;
	m_values[3*node_nb + 1] = y_value;
	m_values[3*node_nb + 2] = z_value;
}

void VectorField::set_value(size_t node_nb, Vec3D value)
{
	set_value(node_nb, value.get_x(), value.get_y(), value.get_z());
}

double VectorField::get_x_value(size_t node_nb) const
{
	return m_values[3*node_nb];
}

double VectorField::get_y_value(size_t node_nb) const
{
	return m_values[3*node_nb + 1];
}

double VectorField::get_z_value(size_t node_nb) const
{
	return m_values[3*node_nb + 2];
}

Vec3D VectorField::get_value(size_t node_nb) const
{
	return Vec3D(m_values[3*node_nb], m_values[3*node_nb + 1], m_values[3*node_nb + 2]);
}

double* VectorField::get_ptr()
{
	return &(m_values[0]);
}

size_t VectorField::get_element_size() const
{
	return sizeof(double);
}